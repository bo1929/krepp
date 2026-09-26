#!/usr/bin/env bash
#
# End-to-end regression tests for the krepp binary.
#
#   bash test/regression/run_regression.sh [path-to-krepp]
#   UPDATE_GOLDEN=1 bash test/regression/run_regression.sh    # refresh goldens
#
# Builds a small index from the toy references, runs every subcommand against a
# fixed set of queries and compares the results with the golden files in
# test/regression/golden/. Everything runs single threaded so the numbers do
# not depend on scheduling.
#
# What is compared:
#   * TSV output  - comment lines dropped and the remaining rows sorted, because
#                   the row order currently follows hash-map iteration order;
#   * jplace      - whitespace stripped, placements split per query and sorted,
#                   with the invocation string masked out;
#   * index files - the reference list and the text metadata (minus the date),
#                   plus size and checksum of the binary parts, which pins the
#                   on-disk format;
#   * exit codes of the failure paths.
set -u

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
krepp="${1:-$root/krepp}"
data="$here/data"
golden="$here/golden"
update="${UPDATE_GOLDEN:-0}"

if [ ! -x "$krepp" ]; then
  echo "error: krepp binary not found at $krepp (build it with 'make')" >&2
  exit 2
fi
krepp="$(cd "$(dirname "$krepp")" && pwd)/$(basename "$krepp")"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/krepp-regression.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

pass=0
fail=0
skip=0
failures=()

note_ok() { pass=$((pass + 1)); }
note_skip() {
  skip=$((skip + 1))
  echo "skip: $1" >&2
}
note_bad() {
  fail=$((fail + 1))
  failures+=("$1")
  echo "FAIL: $1" >&2
}

if [ "$update" = "1" ]; then
  mkdir -p "$golden"
fi

# ---------------------------------------------------------------- normalize

norm_tsv() { grep -v '^#' "$1" | grep -v '^$' | LC_ALL=C sort; }
# dist gained a p-value column; the goldens were produced before it existed and
# pin the query, the reference and the distance, so that column is dropped here
# and checked on its own below.
norm_dist_tsv() { norm_tsv "$1" | awk -F'\t' 'NF>=4 {print $1"\t"$2"\t"$3; next} {print}'; }
# The release version and the build date are metadata, not behaviour: masking
# them keeps the goldens valid across a version bump. Everything else in these
# dumps - the configuration, the counts, the paths - is compared as it is.
norm_text() {
  grep -v '^date:' "$1" |
    sed -E 's/^krepp version: .*/krepp version: VERSION/' |
    grep -v '^$' |
    LC_ALL=C sort
}

# The tabular placement reports a like-weight-ratio that is a ratio of sums
# computed in hash-map order, so its last digits move between runs. The edge
# number and the distance are exact and are compared as they are.
norm_place_tsv() {
  grep -v '^#' "$1" | grep -v '^$' |
    awk -F'\t' 'NF==5 {printf "%s\t%s\t%s\t%.3f\t%s\n", $1, $2, $3, $4, $5; next} {print}' |
    LC_ALL=C sort
}

# jplace: canonical JSON with the unstable parts (placement order, invocation
# and like-weight-ratio) normalised away. Needs python3; when it is missing the
# jplace check is skipped rather than reported as a failure.
HAVE_PYTHON=0
if command -v python3 > /dev/null 2>&1; then HAVE_PYTHON=1; fi
norm_jplace() { python3 "$here/normalize_jplace.py" "$1"; }

norm_sizehash() { printf '%s\n' "$(wc -c < "$1" | tr -d ' ')" "$(cksum < "$1")"; }

# check <label> <normalizer> <actual-file> <golden-name>
check() {
  local label="$1" normalizer="$2" actual="$3" name="$4"
  if [ ! -f "$actual" ]; then
    note_bad "$label (no output file)"
    return
  fi
  "$normalizer" "$actual" > "$tmp/norm.actual"
  if [ "$update" = "1" ]; then
    cp "$tmp/norm.actual" "$golden/$name"
    note_ok
    return
  fi
  if [ ! -f "$golden/$name" ]; then
    note_bad "$label (missing golden $name)"
    return
  fi
  # Goldens are stored already normalized.
  if diff -u "$golden/$name" "$tmp/norm.actual" > "$tmp/diff.txt"; then
    note_ok
  else
    note_bad "$label"
    head -40 "$tmp/diff.txt" >&2
  fi
}

# --------------------------------------------------------------- the corpus

# The toy references are checked in as a tarball of .xz members; the extracted
# .fna files are gitignored, so a fresh clone has to unpack them first.
toy_dir="$root/test/references_toy"
if [ ! -f "$toy_dir/G000341695.fna" ]; then
  if [ ! -f "$root/test/references_toy.tar.gz" ]; then
    echo "error: neither $toy_dir nor test/references_toy.tar.gz is present" >&2
    exit 2
  fi
  echo "unpacking test/references_toy.tar.gz ..." >&2
  tar -xzf "$root/test/references_toy.tar.gz" -C "$root/test" || exit 2
  for member in "$toy_dir"/*.fna.xz; do
    [ -f "$member" ] || continue
    if command -v xz > /dev/null 2>&1; then
      xz -d "$member" || exit 2
    else
      echo "error: xz is needed to unpack $member" >&2
      exit 2
    fi
  done
fi

mkdir -p "$tmp/work"
ln -s "$root/test/references_toy" "$tmp/work/references_toy"
ln -s "$root/test/tree_toy.nwk" "$tmp/work/tree_toy.nwk"
ln -s "$root/test/lineages_toy.txt" "$tmp/work/lineages_toy.txt"
ln -s "$data/query_reads.fq" "$tmp/work/query_reads.fq"
ln -s "$data/query_contigs.fa" "$tmp/work/query_contigs.fa"
cp "$root/test/input_map.tsv" "$tmp/work/input_map.tsv"
cd "$tmp/work"

KREPP=("$krepp" --num-threads 1)

# ------------------------------------------------------------ index build

if "${KREPP[@]}" index -i input_map.tsv -o idx -t tree_toy.nwk -k 27 -w 35 -h 11 > index.log 2>&1; then
  note_ok
else
  note_bad "index exits successfully"
  cat index.log >&2
fi

suffix="-m4r1-frac"
for f in "cmer$suffix" "inc$suffix" "crecord$suffix" "metadata$suffix" "metadata$suffix.txt" "reflist$suffix" "tree$suffix"; do
  if [ -f "idx/$f" ]; then note_ok; else note_bad "index produces $f"; fi
done

check "index metadata" norm_text "idx/metadata$suffix.txt" "metadata.txt"
check "index reflist" cat "idx/reflist$suffix" "reflist.txt"
check "index cmer bytes" norm_sizehash "idx/cmer$suffix" "cmer.sizehash"
check "index inc bytes" norm_sizehash "idx/inc$suffix" "inc.sizehash"
check "index crecord bytes" norm_sizehash "idx/crecord$suffix" "crecord.sizehash"

# ---------------------------------------------------------------- inspect

if "${KREPP[@]}" inspect -i idx > inspect.txt 2>/dev/null; then note_ok; else note_bad "inspect exits successfully"; fi
grep -v '^0\sMER_COUNT' inspect.txt | grep -v '^0\sOUTDEGREE' | grep -v '^0\sNUM_COLORS' > inspect.filtered
check "inspect output" norm_text inspect.filtered "inspect.txt"

# ------------------------------------------------------------------- dist

run_dist() { # run_dist <golden-name> <label> [extra args...]
  local name="$1" label="$2"
  shift 2
  "${KREPP[@]}" dist -i idx "$@" -o out.tsv 2>/dev/null
  check "$label" norm_dist_tsv out.tsv "$name"
}

run_dist dist_reads.tsv "dist reads" -q query_reads.fq

# --no-multi reports the single best reference. Which reference wins a tie
# depends on hash-map iteration order today, so the check is semantic: one row
# per query, and a distance that matches the best distance of the multi output.
"${KREPP[@]}" dist -i idx -q query_reads.fq --no-multi -o out_single.tsv 2>/dev/null
if awk -F'\t' '
  FNR==NR && $0 !~ /^#/ && $1 != "SEQ_ID" && NF==4 { if (!($1 in best) || $3+0 < best[$1]) best[$1] = $3+0; next }
  $0 !~ /^#/ && $1 != "SEQ_ID" && NF==4 {
    rows++
    if (!($1 in best)) { print "unknown query " $1 > "/dev/stderr"; bad = 1; next }
    if (($3+0) - best[$1] > 1e-9) { print "distance mismatch for " $1 > "/dev/stderr"; bad = 1 }
    seen[$1]++
  }
  END {
    for (q in seen) if (seen[q] != 1) { print "duplicate rows for " q > "/dev/stderr"; bad = 1 }
    if (rows == 0) { print "no rows" > "/dev/stderr"; bad = 1 }
    exit bad
  }' out.tsv out_single.tsv 2> "$tmp/awk.err"; then
  # out_single.tsv must cover the same queries as the multi output.
  if [ "$(cut -f1 out_single.tsv | grep -vc '^#')" -eq "$(cut -f1 out.tsv | grep -v '^#' | sort -u | wc -l | tr -d ' ')" ]; then
    note_ok
  else
    note_bad "dist reads --no-multi (query coverage)"
    cat "$tmp/awk.err" >&2
  fi
else
  note_bad "dist reads --no-multi"
  cat "$tmp/awk.err" >&2
fi
run_dist dist_reads_filter.tsv "dist reads --filter" -q query_reads.fq --filter
run_dist dist_reads_max.tsv "dist reads --dist-max" -q query_reads.fq --dist-max 0.1
run_dist dist_reads_summary.tsv "dist reads --summarize" -q query_reads.fq --summarize
run_dist dist_contigs.tsv "dist contigs --hdist-th 2" -q query_contigs.fa --hdist-th 2

# The p-value column: four fields per row, a probability in [0, 1], zero for the
# closest hit of each query (it is compared with itself), and above zero for at
# least one row of a query whose hits span more than one distance.
"${KREPP[@]}" dist -i idx -q query_reads.fq -o out.tsv 2>/dev/null
if awk -F'\t' '
  /^#/ { next }
  $1 == "SEQ_ID" { if (NF != 4) { print "header has " NF " fields" > "/dev/stderr"; bad = 1 } next }
  NF != 4 { print "row with " NF " fields: " $0 > "/dev/stderr"; bad = 1; next }
  {
    p = $4 + 0
    if (!(p >= 0 && p <= 1)) { print "p-value out of range: " $0 > "/dev/stderr"; bad = 1 }
    if ($3 == "NaN") next
    if (!($1 in best) || ($3 + 0) < best[$1]) { best[$1] = $3 + 0; best_p[$1] = p }
    distinct[$1 "\t" ($3 + 0)] = 1
    if (p > 0) { above_zero[$1] = 1 }
  }
  END {
    # Tied distances are not distinguishable from the best hit, so only a query
    # whose rows span more than one distance must have a row above zero.
    for (k in distinct) { split(k, a, "\t"); ndist[a[1]]++ }
    for (q in best_p) {
      if (best_p[q] != 0) { print "closest hit of " q " has P_VALUE " best_p[q] > "/dev/stderr"; bad = 1 }
      if (ndist[q] > 1 && !(q in above_zero)) { print "no row above P_VALUE 0 for " q > "/dev/stderr"; bad = 1 }
    }
    if (!length(best_p)) { print "no rows" > "/dev/stderr"; bad = 1 }
    exit bad
  }' out.tsv 2> "$tmp/awk.err"; then
  note_ok
else
  note_bad "dist p-value column"
  cat "$tmp/awk.err" >&2
fi

# An index built by an older release must stay readable and match the same
# k-mers. test/index_bench was written by krepp v0.8.5 but is gitignored (it is
# a benchmark artefact, not a source file), so the check is skipped when it is
# not there - the format itself is still pinned by the size/checksum goldens
# above, which were produced by the previous release.
if [ -d "$root/test/index_bench" ]; then
  "${KREPP[@]}" dist -i "$root/test/index_bench" -q "$root/test/query_toy.fq" -o out.tsv 2>/dev/null
  check "dist over a legacy index" norm_dist_tsv out.tsv "dist_legacy.tsv"
else
  note_skip "test/index_bench is absent, so the legacy-index query is not compared"
fi

# ------------------------------------------------------------------ place

run_place() { # run_place <golden-name> <label> [extra args...]
  local name="$1" label="$2"
  shift 2
  "${KREPP[@]}" place -i idx "$@" -o out.txt 2>/dev/null
  check "$label" norm_place_tsv out.txt "$name"
}

run_place place_reads.tsv "place reads --tabular" -q query_reads.fq --tabular

# Same for the single-placement mode: one row per query and an edge that the
# multi run also reported.
"${KREPP[@]}" place -i idx -q query_reads.fq --tabular --no-multi -o out_single.tsv 2>/dev/null
if awk -F'\t' '
  FNR==NR && $0 !~ /^#/ && NF==5 { key = $1 "\t" $3; edges[key] = 1; queries[$1] = 1; next }
  $0 !~ /^#/ && NF==5 {
    rows++
    key = $1 "\t" $3
    if (!(key in edges)) { print "edge " $3 " was not a candidate for " $1 > "/dev/stderr"; bad = 1 }
    seen[$1]++
  }
  END {
    for (q in queries) if (seen[q] != 1) { print "expected exactly one row for " q > "/dev/stderr"; bad = 1 }
    if (rows == 0) { print "no rows" > "/dev/stderr"; bad = 1 }
    exit bad
  }' out.txt out_single.tsv 2> "$tmp/awk2.err"; then
  note_ok
else
  note_bad "place reads --tabular --no-multi"
  cat "$tmp/awk2.err" >&2
fi
run_place place_contigs.tsv "place contigs --tabular --no-filter" -q query_contigs.fa --tabular --no-filter
run_place place_tau0.tsv "place reads --tau 0" -q query_reads.fq --tau 0 --tabular
run_place place_lineage.tsv "place reads on a lineage" -q query_reads.fq -l lineages_toy.txt --tabular

"${KREPP[@]}" place -i idx -q query_reads.fq -o out.jplace 2>/dev/null
if [ "$HAVE_PYTHON" = "1" ]; then
  check "place reads jplace" norm_jplace out.jplace "place_reads.jplace"
else
  echo "skip: python3 not found, jplace comparison not run" >&2
fi

# ------------------------------------------------------------------- sketch

if "${KREPP[@]}" sketch -i references_toy/G000341695.fna -o ref.sketch -k 27 -w 35 -h 11 > sketch.log 2>&1; then
  note_ok
else
  note_bad "sketch exits successfully"
  cat sketch.log >&2
fi
check "sketch bytes" norm_sizehash ref.sketch "ref.sketch.sizehash"

if "${KREPP[@]}" seek -i ref.sketch -q query_reads.fq -o out.tsv 2>/dev/null; then note_ok; else note_bad "seek exits successfully"; fi
check "seek reads" norm_tsv out.tsv "seek_reads.tsv"

# ------------------------------------------------------------ failure paths

if "${KREPP[@]}" dist -i idx -q "$tmp/does-not-exist.fq" > bad.log 2>&1; then
  note_bad "missing query file is rejected"
else
  note_ok
fi

if "${KREPP[@]}" place -i idx -q query_reads.fq --tau 9 > bad2.log 2>&1; then
  note_bad "tau above the hamming threshold is rejected"
else
  note_ok
fi

# An index directory holding no index files used to reach the query loop with a
# null LSH and crash.
mkdir -p empty_index
if "${KREPP[@]}" dist -i empty_index -q query_reads.fq > bad3.log 2>&1; then
  note_bad "an empty index directory is rejected"
else
  if grep -q "No index found" bad3.log; then note_ok; else note_bad "empty index directory message"; fi
fi

# ------------------------------------------------------------------- report

echo
if [ "$update" = "1" ]; then
  echo "regression: goldens updated ($pass checks)"
  exit 0
fi
if [ "$skip" -ne 0 ]; then
  echo "regression: $pass passed, $fail failed, $skip skipped"
else
  echo "regression: $pass passed, $fail failed"
fi
if [ "$fail" -ne 0 ]; then
  printf '  - %s\n' "${failures[@]}" >&2
  exit 1
fi
echo "regression: OK"
