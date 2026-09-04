#!/usr/bin/env python3
"""Benchmark krepp read-level distance estimation (`dist`) and phylogenetic
placement (`place`) on artificial reads with known ground truth.

Artificial reads are simulated from toy reference genomes at multiple
divergence levels (i.e., read-level distances) by applying a known per-base
substitution rate. The ground truth (read -> source genome, true divergence)
is stored alongside the reads, so the benchmark can reuse committed data
instead of regenerating it on every run.

Each configuration is named after its actual command-line options (e.g.
"dist --hdist-th 4 -q reads.fq.gz"). The full query configuration matrix runs
on the baseline index; a set of core configurations (default `dist` and
`place --tabular`) is additionally run on every index configuration to explore
indexing parameters (k-mer length, window length, LSH positions, sketching).

Usage:
  python3 test/benchmark_distances.py                 # ensure data, run, summarize
  python3 test/benchmark_distances.py generate        # (re)generate the data
  python3 test/benchmark_distances.py generate --force
  python3 test/benchmark_distances.py run [--sub dist|place|all] [--out results.tsv]

Data (committed): test/benchmark/reads.fq.gz (FASTQ+gzip, all reads),
test/benchmark/reads.fq (FASTQ plain), test/benchmark/reads.fa (FASTA plain),
test/benchmark/reads.fa.gz (FASTA+gzip), test/benchmark/ground_truth.tsv.
Indexes (not committed): built on demand under test/index_benchmark/.
"""

import argparse
import csv
import gzip
import json
import math
import random
import re
import subprocess
import sys
import tarfile
import time
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
DATA_DIR = TEST_DIR / "benchmark"
INDEX_ROOT = TEST_DIR / "index_benchmark"
KREPP_DEFAULT = TEST_DIR.parent / "krepp"

SEED = 42
READ_LEN = 150
READS_PER_LEVEL = 40
LEVELS = [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25]
SOURCES = ["G000341695", "G001917855"]  # tight-clade member vs. isolated basal

QUERY_FILES = {"fq-gz": DATA_DIR / "reads.fq.gz", "fa-gz": DATA_DIR / "reads.fa.gz"}

# Index configurations. The first entry is the baseline index on which the
# full query matrix is run; every entry also receives the core configurations
# (default `dist` and `place --tabular`) for the index-parameter comparison.
# Constraints: 3 <= h <= 15, k-h <= 16, w >= k; h > 13 explodes memory
# (2^(2h) 64-bit offsets), so h is capped at 13 here.
INDEX_CONFIGS = [
    ["-k", "27", "-w", "35", "-h", "11"],  # baseline (README quickstart)
    [],  # all defaults (k=29, w=k+6, h=k-16)
    ["-k", "19", "-w", "25", "-h", "7"],
    ["-k", "23", "-w", "29", "-h", "9"],
    ["-k", "25", "-w", "31", "-h", "9"],
    ["-k", "27", "-w", "35", "-h", "12"],
    ["-k", "27", "-w", "35", "-h", "13"],
    ["-k", "27", "-w", "28", "-h", "11"],
    ["-k", "27", "-w", "60", "-h", "11"],
    ["-k", "27", "-w", "35", "-h", "11", "--no-frac"],
    ["-k", "27", "-w", "35", "-h", "11", "-m", "8"],
    ["-k", "27", "-w", "35", "-h", "11", "-m", "1", "-r", "0"],  # no LSH filtering
]

# Query configurations: (subcommand, global options, query options, query file key).
QUERY_CONFIGS = [
    ("dist", [], [], "fq-gz"),
    ("dist", ["--num-threads", "4"], [], "fq-gz"),
    ("dist", ["--num-threads", "8"], [], "fq-gz"),
    ("dist", ["--seed", "42"], [], "fq-gz"),
    ("dist", ["--verbose"], [], "fq-gz"),
    ("dist", [], ["--hdist-th", "0"], "fq-gz"),
    ("dist", [], ["--hdist-th", "2"], "fq-gz"),
    ("dist", [], ["--chisq", "1.0"], "fq-gz"),
    ("dist", [], ["--chisq", "5.0"], "fq-gz"),
    ("dist", [], ["--dist-max", "0.1"], "fq-gz"),
    ("dist", [], ["--no-multi"], "fq-gz"),
    ("dist", [], ["--filter"], "fq-gz"),
    ("dist", [], ["--summarize"], "fq-gz"),
    ("dist", [], [], "fa-gz"),
    ("place", [], ["--tabular"], "fq-gz"),
    ("place", [], [], "fq-gz"),
    ("place", ["--num-threads", "4"], ["--tabular"], "fq-gz"),
    ("place", ["--num-threads", "8"], ["--tabular"], "fq-gz"),
    ("place", ["--seed", "42"], ["--tabular"], "fq-gz"),
    ("place", ["--verbose"], ["--tabular"], "fq-gz"),
    ("place", [], ["-t", "tree_toy.nwk", "--tabular"], "fq-gz"),
    ("place", [], ["-l", "lineages_toy.txt", "--tabular"], "fq-gz"),
    ("place", [], ["--tau", "0", "--tabular"], "fq-gz"),
    ("place", [], ["--tau", "4", "--tabular"], "fq-gz"),
    ("place", [], ["--no-multi", "--tabular"], "fq-gz"),
    ("place", [], ["--no-filter", "--tabular"], "fq-gz"),
    ("place", [], ["--summarize"], "fq-gz"),
    ("place", [], ["--hdist-th", "2", "--tabular"], "fq-gz"),
    ("place", [], ["--tabular"], "fa-gz"),
]

# Run on every index configuration for the index-parameter comparison.
CORE_CONFIGS = [("dist", [], [], "fq-gz"), ("place", [], ["--tabular"], "fq-gz")]


def config_name(subcmd, global_opts, query_opts, query_key):
    return " ".join([subcmd] + global_opts + query_opts + ["-q", QUERY_FILES[query_key].name])


def index_label(index_args):
    return " ".join(index_args) if index_args else "(all defaults)"


def index_short(index_args):
    """Compact label for table headers, e.g. '-k 27 -w 35 -h 11' -> 'k27w35h11'."""
    if not index_args:
        return "defaults"
    parts, i = [], 0
    while i < len(index_args):
        arg = index_args[i].lstrip("-")
        if i + 1 < len(index_args) and not index_args[i + 1].startswith("-"):
            parts.append(arg + index_args[i + 1])
            i += 2
        else:
            parts.append(arg)
            i += 1
    return "-".join(parts)


def index_slug(index_args):
    return "-".join(re.sub(r"[^A-Za-z0-9]+", "", a) for a in index_args) or "defaults"


def index_dir_for(index_args):
    return INDEX_ROOT / index_slug(index_args)


def run_krepp(krepp, args, cwd=TEST_DIR):
    return subprocess.run([str(krepp)] + args, cwd=cwd, check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)


def ensure_references():
    refs = TEST_DIR / "references_toy"
    if not (refs / f"{SOURCES[0]}.fna").exists():
        print("Extracting test/references_toy.tar.gz ...")
        with tarfile.open(TEST_DIR / "references_toy.tar.gz") as tf:
            tf.extractall(TEST_DIR)
        for xz in refs.glob("*.xz"):
            subprocess.run(["xz", "-d", "-k", str(xz)], check=True)


def dir_size_mb(path):
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e6


def ensure_index(krepp, index_dir, index_args):
    """Builds the index if missing; returns (build time in s, size in MB)."""
    if index_dir.exists() and list(index_dir.glob("cmer-*")):
        return None, dir_size_mb(index_dir)
    ensure_references()
    print(f"Building index at {index_dir} ({index_label(index_args) or 'defaults'}) ...")
    index_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    run_krepp(krepp, ["index"] + index_args + ["-o", str(index_dir), "-i", "input_map.tsv", "-t", "tree_toy.nwk", "--num-threads", "8"])
    build_s = time.perf_counter() - t0
    # A stale/incomplete index dir cannot be loaded by queries.
    if not list(index_dir.glob("cmer-*")):
        sys.exit(f"Index build failed at {index_dir}")
    return build_s, dir_size_mb(index_dir)


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def mutate(seq, rate, rng):
    bases = "ACGT"
    out = []
    for b in seq:
        if rng.random() < rate:
            b = rng.choice([x for x in bases if x != b])
        out.append(b)
    return "".join(out)


def generate(force=False):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    marker = QUERY_FILES["fq-gz"]
    if marker.exists() and not force:
        print(f"Data already exists at {marker} (use --force to regenerate).")
        return
    ensure_references()
    rng = random.Random(SEED)
    genomes = {
        src: "".join(c for c in "".join(read_fasta(TEST_DIR / "references_toy" / f"{src}.fna").values()).upper() if c in "ACGT") for src in SOURCES
    }
    n_reads = READS_PER_LEVEL * len(LEVELS) * len(SOURCES)
    records = []  # (read_id, source, level, seq, qual)
    for src in SOURCES:
        genome = genomes[src]
        for level in LEVELS:
            for rep in range(READS_PER_LEVEL):
                start = rng.randrange(0, len(genome) - READ_LEN)
                seq = mutate(genome[start : start + READ_LEN], level, rng)
                read_id = f"b_{src}_d{level:g}_{rep}"
                records.append((read_id, src, level, seq, "I" * READ_LEN))
    assert len(records) == n_reads

    def write_fasta(fh, recs):
        for rid, _, _, seq, _ in recs:
            fh.write(f">{rid}\n{seq}\n")

    def write_fastq(fh, recs):
        for rid, _, _, seq, qual in recs:
            fh.write(f"@{rid}\n{seq}\n+\n{qual}\n")

    with gzip.open(QUERY_FILES["fq-gz"], "wt") as fh:
        write_fastq(fh, records)
    with gzip.open(QUERY_FILES["fa-gz"], "wt") as fh:
        write_fasta(fh, records)
    with open(DATA_DIR / "ground_truth.tsv", "w") as fh:
        fh.write("READ_ID\tSOURCE\tTRUE_DIV\n")
        for rid, src, level, _, _ in records:
            fh.write(f"{rid}\t{src}\t{level}\n")
    print(f"Generated {n_reads} reads " f"({len(SOURCES)} sources x {len(LEVELS)} levels x {READS_PER_LEVEL} reps) " f"in {DATA_DIR}")


def load_ground_truth():
    gt = {}
    with open(DATA_DIR / "ground_truth.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gt[row["READ_ID"]] = (row["SOURCE"], float(row["TRUE_DIV"]))
    return gt


# --- tree handling (for placement accuracy) -------------------------------


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name):
        self.name, self.length, self.children, self.parent = name, 0.0, [], None


def parse_newick(text):
    s, i = text.strip().rstrip(";"), 0

    def parse_node():
        nonlocal i
        node = Node("")
        if s[i] == "(":
            i += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if s[i] == ",":
                    i += 1
                elif s[i] == ")":
                    i += 1
                    break
        name = ""
        while i < len(s) and s[i] not in ",():":
            name += s[i]
            i += 1
        if i < len(s) and s[i] == ":":
            i += 1
            num = ""
            while i < len(s) and s[i] not in ",():":
                num += s[i]
                i += 1
            node.length = float(num)
        node.name = name
        return node

    return parse_node()


def leaf_set(node, memo):
    if node.name not in memo:
        if node.children:
            leaves = set()
            for ch in node.children:
                leaves |= leaf_set(ch, memo)
            memo[node.name] = leaves
        else:
            memo[node.name] = {node.name}
    return memo[node.name]


class Tree:
    def __init__(self, nwk_path):
        self.root = parse_newick(nwk_path.read_text())
        self.memo = {}
        self._cache = {}

    def contains(self, node_name, taxon):
        if node_name not in self._cache:
            self._cache[node_name] = self.by_name(node_name)
        node = self._cache[node_name]
        return node is not None and taxon in leaf_set(node, self.memo)

    def by_name(self, name):
        if name not in self._cache:
            self._cache[name] = self._find(self.root, name)
        return self._cache[name]

    def _find(self, node, name):
        if node.name == name:
            return node
        for ch in node.children:
            found = self._find(ch, name)
            if found:
                return found
        return None

    def path_length(self, name_a, name_b):
        """Total branch length of the path between two nodes (via their LCA).

        Returns None if either name is not in the tree."""
        if name_a not in self._cache:
            self._cache[name_a] = self.by_name(name_a)
        if name_b not in self._cache:
            self._cache[name_b] = self.by_name(name_b)
        na, nb = self._cache[name_a], self._cache[name_b]
        if na is None or nb is None:
            return None
        ancestors, dist = {}, 0.0
        cur = na
        while cur is not None:
            ancestors[cur] = dist
            dist += cur.length
            cur = cur.parent
        total, cur = 0.0, nb
        while cur is not None:
            if cur in ancestors:
                return ancestors[cur] + total
            total += cur.length
            cur = cur.parent
        return None


# --- output parsers --------------------------------------------------------


def parse_dist_tsv(out):
    hits = {}
    for line in out.splitlines():
        if line.startswith("#") or line.startswith("SEQ_ID"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        hits.setdefault(parts[0], []).append((parts[1], float(parts[2])))
    return hits


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float("nan")


def mean(vals):
    return sum(vals) / len(vals) if vals else None


def eval_dist(name, out, gt, wall_s):
    """Accuracy of estimated read distances against the simulated divergence.

    For every read, the error is the absolute difference between the estimated
    distance to the source genome and the true divergence used for simulation,
    in raw substitutions per site (no model correction)."""
    hits = parse_dist_tsv(out)
    n = len(gt)
    errs, ests, trues = [], [], []
    per_level = {}  # level -> {"ests": [...], "reads": int}
    for rid, (src, level) in gt.items():
        entry = per_level.setdefault(level, {"ests": [], "reads": 0})
        entry["reads"] += 1
        src_d = [d for r, d in hits.get(rid, []) if r == src]
        if not src_d:
            continue
        est = min(src_d)
        entry["ests"].append(est)
        errs.append(abs(est - level))
        ests.append(est)
        trues.append(level)
    levels = sorted(per_level)
    return {
        "invocation": name,
        "subcommand": "dist",
        "seconds": round(wall_s, 3),
        "reads_per_s": round(n / wall_s, 1),
        "reads": n,
        "reported": round(len(hits) / n, 4),
        "distance_mae": round(mean(errs), 4) if errs else None,
        "distance_bias": round(mean(ests) - mean(trues), 4) if ests else None,
        "correlation": round(pearson(trues, ests), 4) if len(ests) > 1 else None,
        "levels": levels,
        "level_coverage": [round(len(per_level[l]["ests"]) / per_level[l]["reads"], 4) for l in levels],
        "level_mean_est": [round(mean(per_level[l]["ests"]), 4) if per_level[l]["ests"] else None for l in levels],
        "level_bias": [round(mean(per_level[l]["ests"]) - l, 4) if per_level[l]["ests"] else None for l in levels],
        "level_mae": [round(mean([abs(e - l) for e in per_level[l]["ests"]]), 4) if per_level[l]["ests"] else None for l in levels],
    }


def eval_place_tabular(name, out, tree, gt, wall_s):
    """Placement accuracy as the tree (edge) distance between the best
    placement (highest like-weight ratio) and the true source genome,
    in substitutions per site along tree branches."""
    by_read = {}
    for line in out.splitlines():
        if line.startswith("#") or line.startswith("SEQ_ID"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        by_read.setdefault(parts[0], []).append((parts[1], int(parts[2]), float(parts[3]), float(parts[4])))
    return place_metrics(name, by_read, tree, gt, wall_s, n_placed_reads=len(by_read))


def place_metrics(name, by_read, tree, gt, wall_s, n_placed_reads):
    n = len(gt)
    lwr_true = 0
    edge_ds = []
    per_level = {}  # level -> {"edges": [...], "lwr": [...], "reads": int}
    for rid, (src, level) in gt.items():
        entry = per_level.setdefault(level, {"edges": [], "reads": 0})
        entry["reads"] += 1
        rows = by_read.get(rid, [])
        if not rows:
            continue
        best = max(rows, key=lambda r: r[2])  # highest like-weight ratio
        src_rows = [r for r in rows if tree.contains(r[0], src)]
        if src_rows:
            lwr_true += sum(r[2] for r in src_rows)
        d_path = tree.path_length(best[0], src)
        if d_path is not None:
            edge_ds.append(d_path)
            entry["edges"].append(d_path)
    levels = sorted(per_level)
    return {
        "invocation": name,
        "subcommand": "place",
        "seconds": round(wall_s, 3),
        "reads_per_s": round(n / wall_s, 1),
        "reads": n,
        "reported": round(n_placed_reads / n, 4),
        "source_lwr": round(lwr_true / n, 4),
        "edge_distance_mean": round(mean(edge_ds), 4) if edge_ds else None,
        "edge_distance_median": round(sorted(edge_ds)[len(edge_ds) // 2], 4) if edge_ds else None,
        "edge_distance_measured": round(len(edge_ds) / n, 4),
        "levels": levels,
        "level_coverage": [round(len(per_level[l]["edges"]) / per_level[l]["reads"], 4) for l in levels],
        "level_edge_mean": [round(mean(per_level[l]["edges"]), 4) if per_level[l]["edges"] else None for l in levels],
    }


def eval_place_summarize(name, out, tree, gt, wall_s):
    abundance = {}
    for line in out.splitlines():
        if line.startswith("#") or line.startswith("DISTAL_NODE"):
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            abundance[parts[0]] = float(parts[3])
    srcs = {src for src, _ in gt.values()}
    src_ab = sum(ab for node, ab in abundance.items() if any(tree.contains(node, s) for s in srcs))
    return {"invocation": name, "subcommand": "place", "seconds": round(wall_s, 3), "reads": len(gt), "source_abundance": round(src_ab, 4)}


def eval_dist_summarize(name, out, gt, wall_s):
    abundance = {}
    for line in out.splitlines():
        if line.startswith("#") or line.startswith("REFERENCE_NAME"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            abundance[parts[0]] = float(parts[2])
    srcs = {src for src, _ in gt.values()}
    return {
        "invocation": name,
        "subcommand": "dist",
        "seconds": round(wall_s, 3),
        "reads": len(gt),
        "source_abundance": round(sum(abundance.get(s, 0.0) for s in srcs), 4),
    }


def eval_place_jplace(name, out, tree, gt, wall_s):
    jplace = json.loads(out)
    fields = jplace["fields"]
    lwr_ix = fields.index("like_weight_ratio")
    # Map edge numbers to node names using the decorated tree emitted in-file
    # (node names precede the branch length: e.g. "G000341695:0.0022{39}").
    edge_to_node = {int(e): nm for nm, e in re.findall(r"([A-Za-z0-9_.|-]+):[^,():]*\{(\d+)\}", jplace["tree"])}
    by_read = {}
    for pl in jplace["placements"]:
        rid = pl["n"][0]
        # Same 4-field row shape as tabular output: (node, edge, LWR, distance).
        by_read[rid] = [(edge_to_node.get(p[0], "\0"), p[0], p[lwr_ix], p[fields.index("distance")]) for p in pl["p"]]
    return place_metrics(name, by_read, tree, gt, wall_s, n_placed_reads=len(by_read))


NUMERIC_KEYS = [
    "seconds",
    "reads_per_s",
    "source_lwr",
    "edge_distance_mean",
    "edge_distance_median",
    "distance_mae",
    "distance_bias",
    "correlation",
    "source_abundance",
]
SD_KEYS = ["source_lwr", "edge_distance_mean", "distance_mae"]


def aggregate(per_config):
    """Average metrics across repeats of the same (index, configuration) pair."""
    results = []
    for (_, _), runs in per_config.items():
        res = dict(runs[-1])
        for key in NUMERIC_KEYS:
            vals = [r[key] for r in runs if r.get(key) is not None]
            if vals:
                res[key] = round(sum(vals) / len(vals), 4)
                if len(vals) > 1 and key in SD_KEYS:
                    mean = res[key]
                    res[key + "_sd"] = round(math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)), 4)
        results.append(res)
    return results


def evaluate(name, subcmd, out, gt, tree, wall_s):
    if "--summarize" in name:
        return eval_dist_summarize(name, out, gt, wall_s) if subcmd == "dist" else eval_place_summarize(name, out, tree, gt, wall_s)
    if subcmd == "dist":
        return eval_dist(name, out, gt, wall_s)
    if "--tabular" in name:
        return eval_place_tabular(name, out, tree, gt, wall_s)
    return eval_place_jplace(name, out, tree, gt, wall_s)


def brief(res):
    if "source_abundance" in res:
        return f"source abundance={res['source_abundance']}"
    if res["subcommand"] == "dist":
        return f"mean absolute error={res.get('distance_mae')}, " f"bias={res.get('distance_bias')}, correlation={res.get('correlation')}"
    return f"edge distance={res.get('edge_distance_mean')}, " f"source LWR={res.get('source_lwr')}"


def run_configs(krepp, sub, out_path, reps=1):
    gt = load_ground_truth()
    tree = Tree(TEST_DIR / "tree_toy.nwk")
    per_config = {}
    baseline_label = index_label(INDEX_CONFIGS[0])
    total = sum(len([c for c in (QUERY_CONFIGS if i == 0 else CORE_CONFIGS) if sub in ("all", c[0])]) for i in range(len(INDEX_CONFIGS))) * reps
    done = 0
    for ix, index_args in enumerate(INDEX_CONFIGS):
        index_dir = index_dir_for(index_args)
        build_s, size_mb = ensure_index(krepp, index_dir, index_args)
        if build_s is not None:
            print(f"  built in {build_s:.1f}s, size {size_mb:.1f}MB")
        configs = [c for c in (QUERY_CONFIGS if ix == 0 else CORE_CONFIGS) if sub in ("all", c[0])]
        label = index_label(index_args)
        for _ in range(reps):
            for subcmd, global_opts, query_opts, qkey in configs:
                name = config_name(subcmd, global_opts, query_opts, qkey)
                args = global_opts + [subcmd, "-i", str(index_dir), "-q", str(QUERY_FILES[qkey])] + query_opts
                t0 = time.perf_counter()
                proc = run_krepp(krepp, args)
                wall_s = time.perf_counter() - t0
                res = evaluate(name, subcmd, proc.stdout, gt, tree, wall_s)
                res["subcommand"] = subcmd
                res["index"] = label
                res["index_short"] = index_short(index_args)
                res["query_options"] = " ".join(query_opts)
                res["global_options"] = " ".join(global_opts)
                res["index_build_s"] = build_s if build_s is not None else None
                res["index_size_mb"] = round(size_mb, 2) if size_mb is not None else None
                per_config.setdefault((label, name), []).append(res)
                done += 1
                print(f"  [{done:3d}/{total}] {name:55s} {wall_s:6.2f}s  {brief(res)}")
    results = aggregate(per_config)
    if out_path:
        write_results(results, out_path)
    print_summary(results, sub, baseline_label)
    return results


def fmt(value, digits=4):
    return "-" if value is None else f"{value:.{digits}f}"


def print_config_table(rows, title):
    summarize_rows = [r for r in rows if "source_abundance" in r]
    rows = [r for r in rows if "distance_mae" in r or "edge_distance_mean" in r]
    print(f"\n{title}")
    print("-" * 118)
    if rows:
        width = max(len("configuration"), max(len(r["invocation"]) for r in rows))
        if rows[0]["subcommand"] == "dist":
            header = f"{'configuration':<{width}} {'time (s)':>8} {'reads/s':>9} {'mean abs. error':>16} {'mean bias':>10} {'correlation':>12}"
            print(header)
            print("-" * 118)
            for r in rows:
                sd = r.get("distance_mae_sd")
                mae = fmt(r.get("distance_mae")) + (f" ±{sd:.3f}" if sd is not None else "")
                print(
                    f"{r['invocation']:<{width}} {fmt(r.get('seconds'), 2):>8} {fmt(r.get('reads_per_s'), 0):>9} {mae:>16} {fmt(r.get('distance_bias')):>10} {fmt(r.get('correlation')):>12}"
                )
            best = min((r for r in rows if r.get("distance_mae") is not None), key=lambda r: r["distance_mae"], default=None)
            fast = min(rows, key=lambda r: r.get("seconds") or 1e9)
            print("-" * 118)
            if best:
                print(f"lowest error:  {best['invocation']} (mean absolute error={fmt(best.get('distance_mae'))})")
            print(f"fastest:       {fast['invocation']} ({fmt(fast.get('seconds'), 2)}s, {fmt(fast.get('reads_per_s'), 0)} reads/s)")
        else:
            header = f"{'configuration':<{width}} {'time (s)':>8} {'reads/s':>9} {'edge distance':>14} {'edge dist. median':>18} {'source LWR':>12} {'measured':>10}"
            print(header)
            print("-" * 118)
            for r in rows:
                sd = r.get("edge_distance_mean_sd")
                edge = fmt(r.get("edge_distance_mean")) + (f" ±{sd:.3f}" if r.get("edge_distance_mean_sd") is not None else "")
                print(
                    f"{r['invocation']:<{width}} {fmt(r.get('seconds'), 2):>8} {fmt(r.get('reads_per_s'), 0):>9} {edge:>14} {fmt(r.get('edge_distance_median')):>18} {fmt(r.get('source_lwr')):>12} {fmt(r.get('edge_distance_measured')):>10}"
                )
            best = min((r for r in rows if r.get("edge_distance_mean") is not None), key=lambda r: r["edge_distance_mean"], default=None)
            fast = min(rows, key=lambda r: r.get("seconds") or 1e9)
            print("-" * 118)
            if best:
                print(f"closest:       {best['invocation']} (mean edge distance={fmt(best.get('edge_distance_mean'))})")
            print(f"fastest:       {fast['invocation']} ({fmt(fast.get('seconds'), 2)}s, {fmt(fast.get('reads_per_s'), 0)} reads/s)")
    if summarize_rows:
        width = max(len("configuration"), max(len(r["invocation"]) for r in summarize_rows))
        print(f"\n{title} (aggregated with --summarize)")
        print(f"{'configuration':<{width}} {'source abundance':>17}")
        print("-" * (width + 18))
        for r in summarize_rows:
            print(f"{r['invocation']:<{width}} {fmt(r.get('source_abundance')):>17}")


def print_index_table(results):
    core_names = {"dist -q reads.fq.gz", "place --tabular -q reads.fq.gz"}
    rows = [r for r in results if r["invocation"] in core_names]
    by_index = {}
    for r in rows:
        by_index.setdefault(r["index"], {})[r["subcommand"]] = r
    if not by_index:
        return
    width = max(len("index parameters"), max(len(ix) for ix in by_index))
    print("\nINDEX CONFIGURATIONS (default dist and place --tabular queries on each index)")
    print("-" * 130)
    print(
        f"{'index parameters':<{width}} {'build (s)':>10} {'size (MB)':>10} | {'mean abs. error':>16} {'mean bias':>10} | {'edge distance':>14} {'source LWR':>12}"
    )
    print(f"{'':<{width}} {'':>10} {'':>10} | {'(dist)':>16} {'(dist)':>11} | {'(place)':>14} {'(place)':>12}")
    print("-" * 130)
    for ix, subres in by_index.items():
        d = subres.get("dist", {})
        p = subres.get("place", {})
        print(
            f"{ix:<{width}} {fmt(d.get('index_build_s'), 1):>10} {fmt(d.get('index_size_mb'), 1):>10} | {fmt(d.get('distance_mae')):>16} {fmt(d.get('distance_bias')):>11} | {fmt(p.get('edge_distance_mean')):>14} {fmt(p.get('source_lwr')):>12}"
        )
    print("-" * 130)


def print_level_details(results, baseline_label):
    for subcmd in ("dist", "place"):
        rows = [r for r in results if r["subcommand"] == subcmd and r.get("levels") and r["index"] == baseline_label and "source_abundance" not in r]
        if not rows:
            continue
        r = rows[0]
        print(f"\n{subcmd.upper()} accuracy across divergence levels (baseline index: {baseline_label})")
        print("-" * 104)
        if subcmd == "dist":
            print(f"{'true divergence':>16} {'coverage':>10} {'mean estimated':>16} {'mean bias':>12} {'mean abs. error':>16}")
            print("-" * 104)
            for l, c, m, b, e in zip(r["levels"], r["level_coverage"], r["level_mean_est"], r["level_bias"], r["level_mae"]):
                print(f"{l:>16.3f} {fmt(c):>10} {fmt(m):>16} {fmt(b):>12} {fmt(e):>16}")
        else:
            print(f"{'true divergence':>16} {'coverage':>10} {'mean edge distance':>20}")
            print("-" * 104)
            for l, c, m in zip(r["levels"], r["level_coverage"], r["level_edge_mean"]):
                print(f"{l:>16.3f} {fmt(c):>10} {fmt(m):>20}")
        print("-" * 104)


def print_level_matrix(results):
    rows = [r for r in results if r["subcommand"] == "dist" and r.get("levels") and r["invocation"] == "dist -q reads.fq.gz"]
    if len(rows) < 2:
        return
    levels = sorted({l for r in rows for l in r["levels"]})
    col_w = max(max(len(r["index_short"]) for r in rows), 12) + 1
    print("\nDISTANCE MEAN ABSOLUTE ERROR BY TRUE DIVERGENCE ACROSS INDEX CONFIGURATIONS (default dist queries, substitutions per site)")
    print("-" * (16 + col_w * len(rows)))
    print(f"{'true divergence':<16}" + "".join(f"{r['index_short']:>{col_w}}" for r in rows))
    print("-" * (16 + col_w * len(rows)))
    for level in levels:
        line = f"{level:<16.3f}"
        for r in rows:
            if level in r["levels"]:
                line += f"{r['level_mae'][r['levels'].index(level)]:>{col_w}.4f}"
            else:
                line += f"{'-':>{col_w}}"
        print(line)
    print("-" * (16 + col_w * len(rows)))


def write_results(results, out_path):
    keys = [
        "subcommand",
        "index",
        "global_options",
        "query_options",
        "invocation",
        "seconds",
        "reads_per_s",
        "reads",
        "reported",
        "distance_mae",
        "distance_mae_sd",
        "distance_bias",
        "correlation",
        "source_lwr",
        "source_lwr_sd",
        "edge_distance_mean",
        "edge_distance_mean_sd",
        "edge_distance_median",
        "edge_distance_measured",
        "source_abundance",
        "index_build_s",
        "index_size_mb",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    print(f"Wrote results to {out_path}")


def print_summary(results, sub, baseline_label):
    print("\n" + "=" * 130)
    print("BENCHMARK SUMMARY")
    print("=" * 130)
    for subcmd in ("dist", "place"):
        rows = [r for r in results if r["subcommand"] == subcmd and r["index"] == baseline_label]
        if rows and sub in ("all", subcmd):
            print_config_table(rows, f"{subcmd.upper()} configurations (baseline index: {baseline_label})")
    if sub in ("all", "dist", "place"):
        print_index_table(results)
        print_level_details(results, baseline_label)
        print_level_matrix(results)
    print("=" * 130)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", nargs="?", choices=["generate", "run", "all"], default="all")
    ap.add_argument("--force", action="store_true", help="regenerate data even if it exists")
    ap.add_argument("--sub", choices=["dist", "place", "all"], default="all", help="run only one subcommand's configurations")
    ap.add_argument("--krepp", default=str(KREPP_DEFAULT), help="path to the krepp binary")
    ap.add_argument("--out", default=None, help="write per-configuration results to a CSV file")
    ap.add_argument(
        "--reps",
        type=int,
        default=1,
        help="repeats per configuration; metrics are averaged and "
        "accuracy metrics get a standard deviation (krepp output is "
        "not deterministic run-to-run)",
    )
    args = ap.parse_args()
    if args.action in ("generate", "all"):
        generate(force=args.force)
    if args.action in ("run", "all"):
        krepp = Path(args.krepp).resolve()
        if not krepp.exists():
            sys.exit(f"krepp binary not found at {krepp}")
        run_configs(krepp, args.sub, args.out, reps=max(1, args.reps))


if __name__ == "__main__":
    main()
