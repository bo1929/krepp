#!/usr/bin/env bash
# Reproduces the benchmark corpus and the synthetic index used for the
# performance measurements (long contigs and short reads).
#
#   bash test/benchmark/make_bench.sh               # writes under tmp/bench (gitignored)
#
# Produces:
#   tmp/bench/bigidx/            synthetic 4000-taxon index (100 Mbp of mutated 25 kb fragments)
#   tmp/bench/bigidx/tree.nwk    random backbone tree with the same leaf labels
#   tmp/bench/query_1mbp.fa      one 1 Mbp query (first Mbp of G000341695)
#   tmp/bench/query_20x1mbp.fa   20 x 1 Mbp queries (20 Mbp total)
#   tmp/bench/reads_200k_150bp.fq 200k x 150 bp reads (30 Mbp total)
set -euo pipefail
cd "$(dirname "$0")/../.."   # repository root

python3 - <<'PY'
import os, random
random.seed(42)
os.makedirs("tmp/bench/bigidx/genomes", exist_ok=True)

refs = {}
for fn in sorted(os.listdir("test/references_toy")):
    if fn.endswith(".fna"):
        with open(f"test/references_toy/{fn}") as f:
            refs[fn[:-4]] = "".join(l.strip() for l in f if not l.startswith(">")).upper()
names = list(refs)

comp = {"A": "C", "C": "G", "G": "T", "T": "A"}
def mutate(s, rate):
    out = list(s)
    for i, ch in enumerate(out):
        if random.random() < rate:
            out[i] = comp.get(ch, "A")
    return "".join(out)

N, FL = 4000, 25000
lines = []
for i in range(N):
    g = refs[names[i % len(names)]]
    start = random.randrange(0, max(1, len(g) - FL))
    frag = mutate(g[start:start + FL], random.uniform(0.01, 0.12))
    p = f"tmp/bench/bigidx/genomes/S{i:05d}.fna"
    with open(p, "w") as f:
        f.write(f">S{i:05d}\n")
        for j in range(0, len(frag), 80):
            f.write(frag[j:j + 80] + "\n")
    lines.append(f"S{i:05d}\t{p}")
open("tmp/bench/bigidx/input_map.tsv", "w").write("\n".join(lines) + "\n")

# random binary backbone tree over the same labels
nodes = [(n, 1) for n in [f"S{i:05d}" for i in range(N)]]
while len(nodes) > 1:
    a = nodes.pop(random.randrange(len(nodes)))
    b = nodes.pop(random.randrange(len(nodes)))
    nodes.append((f"({a[0]}:{random.uniform(0.01, 0.3):.4f},{b[0]}:{random.uniform(0.01, 0.3):.4f})", a[1] + b[1]))
open("tmp/bench/bigidx/tree.nwk", "w").write(nodes[0][0] + ";\n")

src = "".join(l.strip() for l in open("test/references_toy/G000341695.fna") if not l.startswith(">"))
def write_fa(path, seq, name):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")
write_fa("tmp/bench/query_1mbp.fa", src[:1000000], "query_1Mbp")
with open("tmp/bench/query_20x1mbp.fa", "w") as f:
    for k in range(20):
        f.write(f">q{k}\n")
        for i in range(0, 1000000, 80):
            f.write(src[i:i + 80] + "\n")
with open("tmp/bench/reads_200k_150bp.fq", "w") as f:
    for i in range(200000):
        p = random.randrange(0, len(src) - 150)
        f.write(f"@r{i}\n{src[p:p + 150]}\n+\n{'I' * 150}\n")
print("benchmark inputs written")
PY

./krepp index -i tmp/bench/bigidx/input_map.tsv -o tmp/bench/bigidx/index -k 27 -w 35 -h 11

cat <<'EOF'
Now run, e.g.:
  ./krepp dist  -i tmp/bench/bigidx/index -q tmp/bench/query_20x1mbp.fa -o /dev/null
  ./krepp place -i tmp/bench/bigidx/index -q tmp/bench/query_20x1mbp.fa -t tmp/bench/bigidx/tree.nwk -o /dev/null
  ./krepp dist  -i tmp/bench/bigidx/index -q tmp/bench/reads_200k_150bp.fq -o /dev/null
EOF
