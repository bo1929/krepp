/* krepp must keep working after the OpenMP runtime has shut down and started
 * again within one process. Build and run with `make test-omp-runtime`.
 *
 * LLVM's libomp, when linked statically (a host that cannot depend on a shared
 * libomp at run time does this, e.g. a DuckDB extension on macOS), shuts itself
 * down when the last thread that used it exits, and initializes again on the
 * next call. A `#pragma omp critical` leaves its lock in a variable the compiler
 * places in krepp's own binary: the shutdown frees the lock, the variable keeps
 * pointing at it, and the next critical section dereferences a freed lock and
 * crashes. A host that runs krepp on short-lived threads (one per query, say)
 * hits exactly this sequence.
 *
 * So the main thread here never touches OpenMP. Each round runs on a thread of
 * its own, with several OpenMP threads: it builds an index from per-reference
 * files and another from one multi-record FASTA (index_files() and
 * index_sequences()), then loads the first and places and measures queries
 * through api::Session. Every round's thread exit is the runtime's last-thread
 * exit, so every round after the first starts a fresh runtime and enters each of
 * those paths' locks again. The critical sections only the CLI drivers reach
 * (QueryIndex in krepp.cpp, SBatch::seek_sequences) are not covered here.
 *
 * Only a static libomp shuts down this way. Linked against a shared libomp or
 * libgomp, or built without OpenMP, the test passes either way. To see it fail,
 * build from clean (the makefile does not relink when only LDOMP changes) and
 * link statically:
 *   make clean && make test-omp-runtime LDOMP=/path/to/libomp.a OMP_INC=-I/path/to/include
 *
 * The test is self-contained: it writes its corpus under a directory named for
 * this process and removes it at the end. A crash leaves that directory behind. */

#include "../src/api.hpp"
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <unistd.h>

static const uint32_t NUM_REFERENCES = 4;
static const uint32_t REFERENCE_LENGTH = 20000;
static const uint32_t NUM_QUERIES = 8;
static const uint32_t QUERY_LENGTH = 150;
static const uint32_t NUM_THREADS = 4;
static const uint32_t NUM_ROUNDS = 3;
static const uint32_t SEED = 20261118;

static uint32_t num_checked = 0;
static uint32_t num_failed = 0;

static void expect(bool condition, const std::string& what)
{
  num_checked++;
  std::cerr << (condition ? "[PASS] " : "[FAIL] ") << what << std::endl;
  if (!condition) num_failed++;
}

static void on_thread_that_exits(const std::function<void()>& fn)
{
  std::thread worker(fn);
  worker.join();
}

/* Four random references under ((REF0,REF1),(REF2,REF3)), and queries cut from
 * them, so that every query has somewhere to go. */
static void write_corpus(const std::filesystem::path& root)
{
  std::mt19937 rgen(SEED);
  std::uniform_int_distribution<uint32_t> nt(0, 3);
  std::filesystem::create_directories(root / "references");
  std::ofstream map_stream(root / "input_map.tsv");
  std::ofstream query_stream(root / "queries.fa");
  std::ofstream fastx_stream(root / "references.fna");
  for (uint32_t rix = 0; rix < NUM_REFERENCES; ++rix) {
    const std::string name = "REF" + std::to_string(rix);
    std::string sequence;
    for (uint32_t bix = 0; bix < REFERENCE_LENGTH; ++bix) {
      sequence.push_back("ACGT"[nt(rgen)]);
    }
    const std::filesystem::path fasta = root / "references" / (name + ".fna");
    std::ofstream(fasta) << ">" << name << "\n" << sequence << "\n";
    map_stream << name << "\t" << fasta.string() << "\n";
    fastx_stream << ">" << name << "\n" << sequence << "\n";
    for (uint32_t qix = 0; qix < NUM_QUERIES / NUM_REFERENCES; ++qix) {
      query_stream << ">q" << rix << "_" << qix << "\n" << sequence.substr(qix * 1000, QUERY_LENGTH) << "\n";
    }
  }
  std::ofstream(root / "backbone.nwk") << "((REF0:0.1,REF1:0.1):0.1,(REF2:0.1,REF3:0.1):0.1);\n";
}

/* With a guide tree and a reference map the build takes index_files(); from one
 * FASTA and no tree, index_sequences(). */
static void build_index(const std::filesystem::path& input,
                        const std::filesystem::path& nwk_path,
                        const std::filesystem::path& index_dir)
{
  IndexConfig config;
  config.input = input;
  config.index_dir = index_dir;
  config.nwk_path = nwk_path;
  config.k = 20;
  config.w = static_cast<uint8_t>(26);
  config.h = static_cast<uint8_t>(9);
  config.m = 2;
  set_num_threads(NUM_THREADS);
  gen.seed(SEED);
  IndexMultiple index(config);
  index.set_nrows();
  index.set_lshf();
  index.read_input_file();
  index.obtain_build_tree();
  index.build_index();
  index.save_index();
}

int main()
{
  const std::filesystem::path root =
    std::filesystem::temp_directory_path() / ("krepp-omp-runtime-test-" + std::to_string(static_cast<long>(getpid())));
  std::error_code ignored;
  std::filesystem::remove_all(root, ignored);
  write_corpus(root);

  for (uint32_t round = 0; round < NUM_ROUNDS; ++round) {
    const std::filesystem::path files_index = root / ("index-files-" + std::to_string(round));
    const std::filesystem::path sequences_index = root / ("index-sequences-" + std::to_string(round));
    uint32_t placed = 0;
    uint32_t measured = 0;
    on_thread_that_exits([&] {
      build_index(root / "input_map.tsv", root / "backbone.nwk", files_index);
      build_index(root / "references.fna", "", sequences_index);
      api::Options options;
      options.num_threads = NUM_THREADS;
      api::Session session(files_index, options);
      for (const auto& result : session.place((root / "queries.fa").string())) {
        placed += result.placed;
      }
      for (const auto& result : session.dist((root / "queries.fa").string())) {
        measured += !result.na;
      }
    });
    const std::string when = round == 0 ? "in the process's first runtime" : "after the runtime restarted";
    const std::string prefix = "round " + std::to_string(round) + ", " + when + ": ";
    expect(!std::filesystem::is_empty(files_index), prefix + "an index builds from per-reference files");
    expect(!std::filesystem::is_empty(sequences_index), prefix + "an index builds from one FASTA");
    expect(placed == NUM_QUERIES, prefix + "every query is placed");
    expect(measured == NUM_QUERIES, prefix + "every query gets a distance");
  }

  std::filesystem::remove_all(root, ignored);
  std::cerr << num_checked - num_failed << "/" << num_checked << " checks passed" << std::endl;
  return num_failed == 0 ? 0 : 1;
}
