#!/usr/bin/env python3
"""Parity tests for the krepp Python bindings against the CLI.

Builds the toy index exactly like test/regression/run_regression.sh, runs
`krepp dist`/`krepp place` variants, and asserts the binding results match
the CLI output value-for-value (within the CLI's printed precision of 5
decimal places). Also sanity-checks matches() and the error paths.

Run from the repository root:  make test-python   (or:)
    PYTHONPATH=python python3 test/python/test_bindings.py
"""

import gzip
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
KREPP = os.path.join(ROOT, "krepp")
DATA = os.path.join(ROOT, "test", "regression", "data")

PASS = 0
FAIL = 0


def note(ok, label, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
        print(f"  ok    {label}")
    else:
        FAIL += 1
        print(f"  FAIL  {label}" + (f"\n          {detail}" if detail else ""))


def close(a, b, tol=1e-5):
    an, bn = math.isnan(a), math.isnan(b)
    if an or bn:
        return an and bn
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def rel_close(a, b, tol=1e-4):
    an, bn = math.isnan(a), math.isnan(b)
    if an or bn:
        return an and bn
    return abs(a - b) <= tol * max(abs(a), abs(b), 1e-300)


# ------------------------------------------------------------------ fixtures

def unpack_toy_references():
    toy_dir = os.path.join(ROOT, "test", "references_toy")
    probe = os.path.join(toy_dir, "G000341695.fna")
    if os.path.exists(probe):
        return
    tarball = os.path.join(ROOT, "test", "references_toy.tar.gz")
    if not os.path.exists(tarball):
        sys.exit("error: neither test/references_toy nor test/references_toy.tar.gz present")
    print("unpacking test/references_toy.tar.gz ...")
    subprocess.run(["tar", "-xzf", tarball, "-C", os.path.join(ROOT, "test")], check=True)
    import glob
    for member in glob.glob(os.path.join(toy_dir, "*.fna.xz")):
        subprocess.run(["xz", "-d", member], check=True)


def build_index(work):
    os.symlink(os.path.join(ROOT, "test", "references_toy"), os.path.join(work, "references_toy"))
    os.symlink(os.path.join(ROOT, "test", "tree_toy.nwk"), os.path.join(work, "tree_toy.nwk"))
    shutil.copy(os.path.join(ROOT, "test", "input_map.tsv"), os.path.join(work, "input_map.tsv"))
    os.symlink(os.path.join(DATA, "query_reads.fq"), os.path.join(work, "query_reads.fq"))
    cmd = [KREPP, "--num-threads", "1", "index",
           "-i", "input_map.tsv", "-o", "idx", "-t", "tree_toy.nwk",
           "-k", "27", "-w", "35", "-h", "11"]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        sys.exit("error: index build failed")
    return os.path.join(work, "idx")


def cli_dist(work, *extra):
    out = os.path.join(work, "out.tsv")
    cmd = [KREPP, "--num-threads", "1", "dist", "-i", "idx", "-q", "query_reads.fq", "-o", "out.tsv", *extra]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    assert r.returncode == 0, f"krepp dist failed: {r.stderr}"
    rows = []
    with open(out) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("SEQ_ID"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return rows


def cli_place_tabular(work, *extra):
    out = os.path.join(work, "out_place.tsv")
    cmd = [KREPP, "--num-threads", "1", "place", "-i", "idx", "-q", "query_reads.fq",
           "-o", "out_place.tsv", "--tabular", *extra]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    assert r.returncode == 0, f"krepp place failed: {r.stderr}"
    rows = []
    with open(out) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("SEQ_ID"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return rows


def cli_place_jplace(work, *extra):
    out = os.path.join(work, "out.jplace")
    cmd = [KREPP, "--num-threads", "1", "place", "-i", "idx", "-q", "query_reads.fq", "-o", "out.jplace", *extra]
    r = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    assert r.returncode == 0, f"krepp place jplace failed: {r.stderr}"
    with open(out) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------- main

def main():
    if not os.path.exists(KREPP):
        sys.exit("error: build the krepp binary first (make)")
    try:
        sys.path.insert(0, os.path.join(ROOT, "python"))
        import krepp
    except ImportError as e:
        sys.exit(f"error: import krepp failed ({e}); run `make python` first")
    print(f"krepp bindings version: {krepp.__version__}")

    unpack_toy_references()
    work = tempfile.mkdtemp(prefix="krepp-bindings.")
    try:
        idx = build_index(work)
        query = os.path.join(work, "query_reads.fq")

        s1 = krepp.Session(idx, krepp.Options(num_threads=1))

        # ---- dist: default (multi, no filter) ----------------------------
        cli_rows = cli_dist(work)
        res = s1.dist(query)
        bind_rows = []
        for qd in res:
            if qd.na:
                bind_rows.append([qd.query, "NA", float("nan")])
            else:
                for r in qd.rows:
                    bind_rows.append([qd.query, r.ref_name, r.dist])
        ok = len(cli_rows) == len(bind_rows)
        detail = f"row count CLI={len(cli_rows)} bindings={len(bind_rows)}"
        if ok:
            for c, b in zip(cli_rows, bind_rows):
                if c[0] != b[0] or c[1] != b[1] or not close(float(c[2]), float(b[2])):
                    ok, detail = False, f"CLI {c} != bindings {b}"
                    break
        note(ok, "dist default parity", detail)

        # ---- dist: --p-value ---------------------------------------------
        s_pv = krepp.Session(idx, krepp.Options(num_threads=1, p_value=True))
        cli_rows = cli_dist(work, "--p-value")
        res = s_pv.dist(query)
        bind_rows = [[qd.query, r.ref_name, r.dist, r.p_value] for qd in res for r in qd.rows]
        ok = len(cli_rows) == len(bind_rows)
        detail = f"row count CLI={len(cli_rows)} bindings={len(bind_rows)}"
        if ok:
            for c, b in zip(cli_rows, bind_rows):
                if (c[0] != b[0] or c[1] != b[1]
                        or not close(float(c[2]), b[2])
                        or not rel_close(float(c[3]), b[3])):
                    ok, detail = False, f"CLI {c} != bindings {b}"
                    break
        note(ok, "dist --p-value parity", detail)

        # ---- dist: --no-multi ---------------------------------------------
        s_nm = krepp.Session(idx, krepp.Options(num_threads=1, multi=False))
        cli_rows = cli_dist(work, "--no-multi")
        res = s_nm.dist(query)
        bind_rows = []
        for qd in res:
            if qd.na:
                bind_rows.append([qd.query, "NA", float("nan")])
            else:
                for r in qd.rows:
                    bind_rows.append([qd.query, r.ref_name, r.dist])
        ok = len(cli_rows) == len(bind_rows)
        detail = f"count CLI={len(cli_rows)} bindings={len(bind_rows)}"
        if ok:
            # Tie-aware: on exactly-tied distances the closest pick depends on
            # hash-map order (per process), so accept any name that the multi
            # run reports for this query at the same distance.
            multi_rows = {}
            for c in cli_dist(work):
                multi_rows.setdefault(c[0], []).append((c[1], float(c[2])))
            for i, (c, b) in enumerate(zip(cli_rows, bind_rows)):
                vals_ok = close(float(c[2]), float(b[2]))
                name_ok = c[1] == b[1] or (b[1], float(c[2])) in [
                    (n, d) for n, d in multi_rows.get(b[0], []) if close(d, float(c[2]))
                ]
                if c[0] != b[0] or not vals_ok or not name_ok:
                    ok, detail = False, f"row {i}: CLI {c} != bindings {b}"
                    break
        note(ok, "dist --no-multi parity", detail if not ok else "")

        # ---- dist: --filter ------------------------------------------------
        s_f = krepp.Session(idx, krepp.Options(num_threads=1, dist_filter=True))
        cli_rows = cli_dist(work, "--filter")
        res = s_f.dist(query)
        bind_rows = [[qd.query, r.ref_name, r.dist] for qd in res if not qd.na for r in qd.rows]
        bind_rows += [[qd.query, "NA", float("nan")] for qd in res if qd.na]
        # --filter can drop every row for a query; the CLI then prints NA/NaN.
        cli_by_query = {}
        for c in cli_rows:
            cli_by_query.setdefault(c[0], []).append(c)
        bind_by_query = {}
        for b in bind_rows:
            bind_by_query.setdefault(b[0], []).append(b)
        ok = set(cli_by_query) == set(bind_by_query)
        if ok:
            for q in cli_by_query:
                c_rows, b_rows = cli_by_query[q], bind_by_query[q]
                if len(c_rows) != len(b_rows) or any(
                        c[1] != b[1] or not close(float(c[2]), float(b[2]) if not isinstance(b[2], float) else b[2])
                        for c, b in zip(c_rows, b_rows)):
                    ok = False
                    break
        note(ok, "dist --filter parity", "" if ok else f"CLI {cli_rows[:4]} vs {bind_rows[:4]}")

        # ---- dist: --dist-max 0.1 ------------------------------------------
        s_dm = krepp.Session(idx, krepp.Options(num_threads=1, dist_max=0.1))
        cli_rows = cli_dist(work, "--dist-max", "0.1")
        res = s_dm.dist(query)
        bind_rows = []
        for qd in res:
            if qd.na:
                bind_rows.append([qd.query, "NA", float("nan")])
            else:
                for r in qd.rows:
                    bind_rows.append([qd.query, r.ref_name, r.dist])
        ok = len(cli_rows) == len(bind_rows) and all(
            c[0] == b[0] and c[1] == b[1] and close(float(c[2]), float(b[2]))
            for c, b in zip(cli_rows, bind_rows))
        note(ok, "dist --dist-max parity",
             "" if ok else f"CLI n={len(cli_rows)} vs bindings n={len(bind_rows)}")

        # ---- place: jplace parity -------------------------------------------
        # Normalisation follows test/regression/normalize_jplace.py: LWR (a
        # ratio of sums accumulated in hash-map order) and the placement order
        # are not stable between runs; edge numbers, lengths, likelihoods and
        # distances are (compared at 4 decimals, as the normaliser rounds).
        jp = cli_place_jplace(work)
        s_pl = krepp.Session(idx, krepp.Options(num_threads=1))
        res = s_pl.place(query)
        cli_by_query = {}
        for entry in jp["placements"]:
            # jplace: every "p" array belongs to every name in "n".
            for name in entry["n"]:
                for p in entry["p"]:
                    cli_by_query.setdefault(name, []).append(
                        (p[0], round(p[1], 4), round(p[2], 4), round(p[3], 4), round(p[5], 4)))
        bind_by_query = {}
        for qp in res:
            if not qp.placed:
                continue
            for p in qp.placements:
                bind_by_query.setdefault(qp.query, []).append(
                    (p.edge_num, round(p.pendant_length, 4), round(p.distal_length, 4),
                     round(p.likelihood, 4), round(p.distance, 4)))
        ok = set(cli_by_query) == set(bind_by_query)
        detail = ""
        if ok:
            for q in cli_by_query:
                c_list = sorted(cli_by_query[q])
                b_list = sorted(bind_by_query[q])
                if len(c_list) != len(b_list):
                    ok, detail = False, f"{q}: {len(c_list)} vs {len(b_list)} placements"
                    break
                for c, b in zip(c_list, b_list):
                    if c[0] != b[0] or any(not close(float(x), float(y), tol=2e-4) for x, y in zip(c[1:], b[1:])):
                        ok, detail = False, f"{q}: CLI {c} != bindings {b}"
                        break
                if not ok:
                    break
        else:
            detail = f"query sets differ: {set(cli_by_query) ^ set(bind_by_query)}"
        note(ok, "place jplace parity", detail)

        # ---- place: tabular parity (same normalisation; LWR skipped) ---------
        cli_rows = cli_place_tabular(work)
        cli_by_query = {}
        for c in cli_rows:  # SEQ_ID DISTAL_NODE EDGE_NUM LWR DIST
            cli_by_query.setdefault(c[0], []).append((c[1], int(c[2]), round(float(c[4]), 4)))
        bind_by_query = {}
        for qp in res:
            if not qp.placed:
                continue
            for p in qp.placements:
                bind_by_query.setdefault(qp.query, []).append(
                    (p.distal_node, p.edge_num, round(p.distance, 4)))
        ok = set(cli_by_query) == set(bind_by_query)
        detail = "" if ok else f"query sets differ: {set(cli_by_query) ^ set(bind_by_query)}"
        if ok:
            for q in cli_by_query:
                c_s = sorted(cli_by_query[q])
                b_s = sorted(bind_by_query[q])
                if len(c_s) != len(b_s):
                    ok, detail = False, f"{q}: {len(c_s)} tabular rows vs {len(b_s)} placements"
                    break
                for c, b in zip(c_s, b_s):
                    if c[0] != b[0] or c[1] != b[1] or not close(c[2], b[2], tol=2e-4):
                        ok, detail = False, f"{q}: CLI {c} != bindings {b}"
                        break
                if not ok:
                    break
        note(ok, "place tabular parity", detail)

        # ---- matches: consistency with dist ----------------------------------
        s_mt = krepp.Session(idx, krepp.Options(num_threads=1))
        mres = s_mt.matches(query)
        dres = s1.dist(query)  # default: multi + no filter -> same node set
        ok = len(mres) == len(dres)
        detail = f"query count matches={len(mres)} dist={len(dres)}"
        if ok:
            for qm, qd in zip(mres, dres):
                if qm.query != qd.query:
                    ok, detail = False, f"query order mismatch {qm.query} vs {qd.query}"
                    break
                if qm.enmers == 0 or qm.onmers > qm.enmers:
                    ok, detail = False, f"{qm.query}: onmers={qm.onmers} enmers={qm.enmers}"
                    break
                m_set = {(r.ref_name, round(r.dist, 5)) for r in qm.rows}
                d_set = {(r.ref_name, round(r.dist, 5)) for r in qd.rows} if not qd.na else set()
                if not qd.na and m_set != d_set:
                    ok, detail = False, f"{qm.query}: match/dist row sets differ: {m_set ^ d_set}"
                    break
                for r in qm.rows:
                    if len(r.hdisthist) != 5 or r.hdist_min > 4 or r.match_count > qm.onmers:
                        ok, detail = False, f"{qm.query}/{r.ref_name}: bad summary fields"
                        break
                if not ok:
                    break
                agg = [sum(r.hdisthist[d] for r in qm.rows) for d in range(5)]
                if len(qm.hdist_hist) != 5 or any(abs(a - b) > 1e-9 for a, b in zip(agg, qm.hdist_hist)):
                    ok, detail = False, f"{qm.query}: hdist_hist {qm.hdist_hist} != row sums {agg}"
                    break
        note(ok, "matches consistency + sanity (incl. aggregate hdist_hist)", detail)

        # ---- threads=2 gives the same set (order may differ) -------------------
        s2 = krepp.Session(idx, krepp.Options(num_threads=2))
        res2 = s2.dist(query)
        set1 = sorted((qd.query, r.ref_name, round(r.dist, 5)) for qd in dres for r in qd.rows)
        set2 = sorted((qd.query, r.ref_name, round(r.dist, 5)) for qd in res2 for r in qd.rows)
        note(set1 == set2, "num_threads=2 dist parity (as sets)",
             "" if set1 == set2 else f"{len(set1)} vs {len(set2)} rows")

        # ---- error paths -------------------------------------------------------
        def raises(fn, *a, **kw):
            try:
                fn(*a, **kw)
                return False
            except ValueError:
                return True

        ok = raises(krepp.Session, "/nonexistent/dir")
        ok &= raises(krepp.Session, idx, krepp.Options(hdist_th=99))
        ok &= raises(krepp.Session, idx, krepp.Options(dist_max=5.0))
        s_err = krepp.Session(idx, krepp.Options(num_threads=1, tau=9, hdist_th=4))
        ok &= raises(s_err.place, query)
        ok &= raises(s1.dist, "/nonexistent/query.fq")
        note(ok, "error paths raise ValueError")

        # ---- corrupt index raises RuntimeError (ErrorRelay), not exit ---------
        bad_idx = os.path.join(work, "bad_idx")
        shutil.copytree(idx, bad_idx)
        cmers = sorted(f for f in os.listdir(bad_idx) if f.startswith("cmer-"))
        with open(os.path.join(bad_idx, cmers[0]), "r+b") as fh:
            fh.truncate(24)  # keeps the magic, truncates the arrays
        try:
            krepp.Session(bad_idx, krepp.Options(num_threads=2))
            ok = False
        except RuntimeError:
            ok = True
        except BaseException as e:  # noqa: BLE001 - any other outcome is a failure
            ok = False
            print(f"          unexpected exception type: {type(e).__name__}: {e}")
        note(ok, "corrupt index raises RuntimeError (process survives)")

        # and the session object is still usable afterwards (handler restored)
        res_after = s1.dist(query)
        note(len(res_after) == len(dres), "existing Session still works after a failed load")

        # ---- total_qseq bookkeeping (s1 ran dist 3x: res + dres + res_after) --
        note(s1.total_qseq == 3 * len(dres),
             "total_qseq accumulates across calls", f"total_qseq={s1.total_qseq}")

        # ================= extended coverage =================
        def run_cli(args, outname):
            cmd = [KREPP, "--num-threads", "1", *args, "-o", outname]
            r = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
            assert r.returncode == 0, f"CLI failed ({args}): {r.stderr[-300:]}"
            rows = []
            with open(os.path.join(work, outname), errors="replace") as fh:
                for line in fh:
                    if line.startswith(("#", "SEQ_ID", "REFERENCE_NAME", "DISTAL_NODE")):
                        continue
                    rows.append(line.rstrip("\n").split("\t"))
            return rows

        def dist_key(rows):  # CLI rows -> comparable (query, ref, rounded dist)
            return sorted((r[0], r[1], round(float(r[2]), 5) if r[2] != "NaN" else "NaN") for r in rows)

        def bind_key(res):
            out = []
            for qd in res:
                if qd.na:
                    out.append((qd.query, "NA", "NaN"))
                else:
                    out += [(qd.query, r.ref_name, round(r.dist, 5)) for r in qd.rows]
            return sorted(out)

        # ---- place via lineage file (tabular parity) --------------------------
        os.symlink(os.path.join(ROOT, "test", "lineages_toy.txt"), os.path.join(work, "lineages_toy.txt"))
        s_lin = krepp.Session(idx, krepp.Options(num_threads=1, lineage_path=os.path.join(work, "lineages_toy.txt")))
        cli_t = {}
        for c in run_cli(["place", "-i", "idx", "-q", "query_reads.fq", "-l", "lineages_toy.txt", "--tabular"], "lin.tsv"):
            cli_t.setdefault(c[0], []).append((c[1], int(c[2]), round(float(c[4]), 4)))
        bind_t = {}
        for qp in s_lin.place(query):
            for p in qp.placements:
                bind_t.setdefault(qp.query, []).append((p.distal_node, p.edge_num, round(p.distance, 4)))
        ok = set(cli_t) == set(bind_t) and all(sorted(cli_t[q]) == sorted(bind_t[q]) for q in cli_t)
        note(ok, "place with lineage tree parity",
             "" if ok else f"{set(cli_t) ^ set(bind_t)}")

        # ---- place via external nwk (jplace parity, normalised) ----------------
        run_cli(["place", "-i", "idx", "-q", "query_reads.fq", "-t", "tree_toy.nwk"], "nwk.jplace")
        with open(os.path.join(work, "nwk.jplace")) as fh:
            jp = json.load(fh)
        cli_j = {}
        for entry in jp["placements"]:
            for nm in entry["n"]:
                for p in entry["p"]:
                    cli_j.setdefault(nm, []).append((p[0], round(p[1], 4), round(p[2], 4), round(p[3], 4), round(p[5], 4)))
        s_nwk = krepp.Session(idx, krepp.Options(num_threads=1, nwk_path=os.path.join(work, "tree_toy.nwk")))
        bind_j = {}
        for qp in s_nwk.place(query):
            for p in qp.placements:
                bind_j.setdefault(qp.query, []).append(
                    (p.edge_num, round(p.pendant_length, 4), round(p.distal_length, 4),
                     round(p.likelihood, 4), round(p.distance, 4)))
        ok = set(cli_j) == set(bind_j)
        if ok:
            for q in cli_j:
                c_s, b_s = sorted(cli_j[q]), sorted(bind_j[q])
                if len(c_s) != len(b_s) or any(
                        c[0] != b[0] or any(not close(float(x), float(y), tol=2e-4) for x, y in zip(c[1:], b[1:]))
                        for c, b in zip(c_s, b_s)):
                    ok = False
                    break
        note(ok, "place with external nwk parity (jplace)")

        # ---- dist --hdist-th 2 parity -------------------------------------------
        cli_rows = run_cli(["dist", "-i", "idx", "-q", "query_reads.fq", "--hdist-th", "2"], "ht2.tsv")
        s_h2 = krepp.Session(idx, krepp.Options(num_threads=1, hdist_th=2))
        note(dist_key(cli_rows) == bind_key(s_h2.dist(query)), "dist --hdist-th 2 parity")

        # ---- gzipped query parity ------------------------------------------------
        with open(query, "rb") as fin, gzip.open(os.path.join(work, "query_reads.fq.gz"), "wb") as fout:
            shutil.copyfileobj(fin, fout)
        cli_rows = run_cli(["dist", "-i", "idx", "-q", "query_reads.fq.gz"], "gz.tsv")
        note(dist_key(cli_rows) == bind_key(s1.dist(os.path.join(work, "query_reads.fq.gz"))),
             "dist on gzipped query parity")

        # ---- multi-batch parity (1000 reads > 76.8 kbp batch limit) --------------
        shutil.copy(os.path.join(ROOT, "test", "query_toy.fq"), os.path.join(work, "query_toy.fq"))
        with open(os.path.join(work, "query_toy.fq")) as fin, open(os.path.join(work, "query_x10.fq"), "w") as fout:
            body = fin.read()
            for rep in range(10):
                fout.write(body.replace("@61435-", f"@x{rep}-61435-"))
        cli_rows = run_cli(["dist", "-i", "idx", "-q", "query_x10.fq"], "x10.tsv")
        res_x10 = s1.dist(os.path.join(work, "query_x10.fq"))
        ok = len(cli_rows) == sum(1 if qd.na else len(qd.rows) for qd in res_x10) \
            and dist_key(cli_rows) == bind_key(res_x10)
        note(ok, f"dist multi-batch parity ({len(res_x10)} queries, {len(cli_rows)} rows)")

        # ---- corrupt index battery (RuntimeError, process survives) ---------------
        # (the truncated-cmer case is covered in the error-path section above)
        def corrupt(label, mutate):
            bad = os.path.join(work, "bad_" + label)
            shutil.copytree(idx, bad)
            mutate(bad)
            try:
                krepp.Session(bad, krepp.Options(num_threads=2))
                note(False, f"corrupt {label} raises RuntimeError", "no exception")
            except RuntimeError:
                note(True, f"corrupt {label} raises RuntimeError")
            except BaseException as e:
                note(False, f"corrupt {label} raises RuntimeError", f"got {type(e).__name__}: {e}")

        def trunc(prefix, n):
            def m(d):
                for f in sorted(os.listdir(d)):
                    if f.startswith(prefix):
                        with open(os.path.join(d, f), "r+b") as fh:
                            fh.truncate(n)
                        return
            return m

        def drop_cmer(d):
            for f in sorted(os.listdir(d)):
                if f.startswith("cmer-"):
                    os.remove(os.path.join(d, f))
                    return

        corrupt("crecord", trunc("crecord-", 8))
        corrupt("metadata", trunc("metadata-", 4))
        corrupt("cmer-missing", drop_cmer)
        empty_dir = os.path.join(work, "empty_dir")
        os.makedirs(empty_dir, exist_ok=True)
        try:
            krepp.Session(empty_dir)
            note(False, "empty index dir raises RuntimeError", "no exception")
        except RuntimeError:
            note(True, "empty index dir raises RuntimeError")

        # ---- in-process determinism (NaN-safe bitwise) -----------------------------
        def bitdist(res):
            return [(q.query, q.na, [(r.ref_name, struct.pack("<d", r.dist), struct.pack("<d", r.p_value))
                                     for r in q.rows]) for q in res]
        note(bitdist(dres) == bitdist(res_after), "repeated dist bit-identical in-process")

    finally:
        shutil.rmtree(work, ignore_errors=True)

    print(f"\nbindings tests: {PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
