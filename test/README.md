# Tests

Three suites, all wired into `make test`:

| target | what it runs |
|---|---|
| `make test-unit` | `test/unit/*.cpp` — the doctest suite, linked against the whole library |
| `make test-regression` | `test/regression/run_regression.sh` — CLI golden-file tests |
| `make test-threads` | `test/omp_index_test.cpp` — the index builder under threads |

```bash
make test              # everything
make test-unit         # just the doctest suite
./build/krepp_tests --list-test-cases
./build/krepp_tests -ts=query            # one suite
./build/krepp_tests -tc="dist reports*"  # one case
make coverage          # clang source-based coverage of src/, report on stdout
```

`make coverage` needs `xcrun llvm-profdata`/`llvm-cov` (Xcode command line tools)
and rebuilds the tree with instrumentation, so run `make clean && make` when you
are done with it.

## Unit and integration suite (`test/unit`)

The doctest framework is a submodule (`external/doctest`, pinned to v2.4.11), so a
checkout needs `git submodule update --init` - the clone instructions in the top level
README already pass `--recurse-submodules`.

`test_helpers.hpp` holds the shared fixtures: temporary directories, FASTA/FASTQ
writers, deterministic random DNA, a corpus-to-index builder, an index loader
that mirrors `TargetIndex::load_index()`, query runners that mirror
`QueryIndex::estimate_distances()`/`place_sequences()`, and a throwing
`error_exit` handler for testing failure paths.

| file | covers |
|---|---|
| `test_common.cpp` | nucleotide tables, encodings, bit tricks, hashes, `error_exit`, `ErrorRelay` |
| `test_lshf.cpp` | LSH position sampling, masks, `compute_hash`/`drop_ppos_*`, golden hash values |
| `test_table.cpp` | `DynHT`/`SDynHT` fill, merge and prune, `FlatHT`/`SFlatHT` layout and round-trip |
| `test_record.cpp` | clade records, colour decomposition, `CRecord` round-trip |
| `test_phytree.cpp` | Newick parsing/printing, jplace decorations, traversal, lineages |
| `test_rqseq.cpp` | `RSeq` minimizer extraction against an independent oracle, `QSeq` batching |
| `test_hll_llh.cpp` | HyperLogLog accuracy, the likelihood model, `Minfo` accumulators |
| `test_index.cpp` | input parsing, both build strategies, validation, metadata, on-disk layout |
| `test_query.cpp` | distance estimation and placement over a real index, filters, NA paths |
| `test_sketch_seek.cpp` | sketch creation/loading and the seek path |

Conventions worth keeping:

* tests are single threaded and seed `gen` explicitly, because the LSH positions
  are drawn from it; index builds are byte-reproducible for a fixed seed (there
  is a test for that);
* assertions on the reported distances are qualitative (ordering, thresholds)
  rather than exact floats, so the numerics can evolve;
* assertions on encodings and hashes *are* exact: both are part of the on-disk
  format and must not drift;
* the depth, colour-decoding and table-move tests assert the *correct*
  behaviour, not the behaviour the code happened to have: those were bugs, they
  are fixed, and the tests now guard the fix.

## Regression suite (`test/regression`)

`run_regression.sh` builds an index from the toy references, runs every
subcommand and compares against the goldens in `golden/`. It normalizes the
parts of the output that are not stable between runs today (row order, the
like-weight-ratio, the invocation string) and compares everything else exactly,
including the size and checksum of the binary index parts. It also runs `dist`
against `test/index_bench`, an index written by krepp v0.8.5, which pins
backwards compatibility of the format and of the hash functions.

```bash
bash test/regression/run_regression.sh ./krepp
UPDATE_GOLDEN=1 bash test/regression/run_regression.sh ./krepp   # refresh goldens
```

The suite unpacks `test/references_toy.tar.gz` on demand (the extracted `.fna`
files are gitignored). One check needs a fixture that is not in the repository:
`test/index_bench`, an index written by krepp v0.8.5, which is only used to
compare `dist` against an index from an older release. When it is missing the
check is skipped and reported as such; the on-disk format is still pinned by the
size/checksum goldens, which were produced by the previous release.

Refresh the goldens only when a change to the *expected* output is intended, and
say why in the commit message.
