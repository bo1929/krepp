/* Tests for the parts of the index builder that only misbehave under threads.
 * Build and run with `make test`.
 *
 * Two things are pinned down here:
 *  - a handler installed with set_error_handler may throw, and the exception has
 *    to reach the caller of build_index() rather than escaping an OpenMP region,
 *    where leaving a structured block is undefined;
 *  - the number of threads must not change which k-mers land in the index.
 *
 * Both build strategies are covered, because they parallelise differently and
 * are guarded separately: index_files() spreads a task graph over the backbone,
 * index_sequences() slices one file across a worksharing loop. So is each of the
 * five critical regions the build path reports progress from, which is why the
 * log sink below throws on a chosen line rather than on the first write.
 *
 * The test is self-contained: it generates its own corpus under a directory
 * named for this process and removes it again at the end. */

#include "../src/index.hpp"
#include <atomic>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <unistd.h>
#if defined(_OPENMP) && _WOPENMP == 1
  #include "omp.h"
#endif

/* Small enough to run in seconds, wide enough that the task graph of
 * index_files() and the slices of index_sequences() have something to spread
 * over threads. */
static const uint32_t NUM_REFERENCES = 8;
static const uint32_t REFERENCE_LENGTH = 50000;
static const uint32_t BUILD_THREADS = 4;
/* set_lshf() draws the LSH positions from `gen`, so every build in this process
 * has to start from the same state or two indexes are not comparable. */
static const uint32_t LSH_SEED = 20260908;

static uint32_t num_checked = 0;
static uint32_t num_failed = 0;

static void expect(bool condition, const std::string& what)
{
  num_checked++;
  std::cerr << (condition ? "[PASS] " : "[FAIL] ") << what << std::endl;
  if (!condition) num_failed++;
}

static void skip(const std::string& what) { std::cerr << "[SKIP] " << what << std::endl; }

class BuildError : public std::runtime_error
{
public:
  BuildError(const std::string& msg, int code)
    : std::runtime_error(msg)
    , code(code)
  {}
  int code;
};

/* What the handler saw when it fired. omp_in_parallel() on its own is not worth
 * asserting on: a runtime that gives the build a one-thread team - setting
 * OMP_THREAD_LIMIT is enough - runs the region inactive, where it reports false
 * even though the unwind worked. The team size tells those two apart, and is
 * also what says whether a "four-thread" build was threaded at all. */
static std::atomic<bool> raised_in_parallel = {false};
static std::atomic<int> raised_team_size = {0};
static int build_team_size = 0;

static void reset_raise_state()
{
  raised_in_parallel = false;
  raised_team_size = 0;
}

static void record_raise()
{
#if defined(_OPENMP) && _WOPENMP == 1
  if (omp_in_parallel()) raised_in_parallel = true;
  raised_team_size = omp_get_num_threads();
#endif
}

/* The handler an embedder would install: turn a fatal error into an exception so
 * that a bad input file cannot take the hosting process down with it. */
static void install_throwing_handler()
{
  set_error_handler([](const std::string& msg, int code) {
    record_raise();
    throw BuildError(msg, code);
  });
}

static std::string balanced_nwk(vec_str_iter first, vec_str_iter last)
{
  if (std::distance(first, last) == 1) {
    return *first + ":0.1";
  }
  vec_str_iter mid = std::next(first, std::distance(first, last) / 2);
  return "(" + balanced_nwk(first, mid) + "," + balanced_nwk(mid, last) + "):0.1";
}

struct Corpus
{
  std::filesystem::path root;
  std::filesystem::path input_map; /* name-to-path TSV; drives index_files() */
  std::filesystem::path fastx;     /* one multi-record FASTA; drives index_sequences() */
  std::filesystem::path nwk_path;
  vec<std::string> names;
  std::filesystem::path path_of(const std::string& name) const { return root / "references" / (name + ".fna"); }
};

/* Random ACGT references under a balanced binary backbone; krepp rejects
 * unifurcations, so the tree has to branch two ways at every internal node. The
 * same sequences are written twice, once per reference and once concatenated,
 * so that both build strategies see identical input. */
static Corpus write_corpus(const std::filesystem::path& root)
{
  Corpus corpus;
  corpus.root = root;
  std::filesystem::create_directories(root / "references");
  std::mt19937 rgen(LSH_SEED);
  std::uniform_int_distribution<uint32_t> nt(0, 3);
  corpus.input_map = root / "input_map.tsv";
  corpus.fastx = root / "references.fna";
  std::ofstream map_stream(corpus.input_map);
  std::ofstream fastx_stream(corpus.fastx);
  for (uint32_t rix = 0; rix < NUM_REFERENCES; ++rix) {
    std::string name = "REF" + std::to_string(rix);
    std::string sequence;
    sequence.reserve(REFERENCE_LENGTH);
    for (uint32_t bix = 0; bix < REFERENCE_LENGTH; ++bix) {
      sequence.push_back("ACGT"[nt(rgen)]);
    }
    std::ofstream fasta_stream(corpus.path_of(name));
    fasta_stream << ">" << name << "\n" << sequence << "\n";
    fasta_stream.close();
    /* Checked: a silently empty corpus would make every negative test below
     * pass for the wrong reason. */
    if (!fasta_stream) {
      throw std::runtime_error("failed to write " + corpus.path_of(name).string());
    }
    fastx_stream << ">" << name << "\n" << sequence << "\n";
    map_stream << name << "\t" << corpus.path_of(name).string() << "\n";
    corpus.names.push_back(name);
  }
  map_stream.close();
  fastx_stream.close();
  corpus.nwk_path = root / "backbone.nwk";
  std::ofstream nwk_stream(corpus.nwk_path);
  nwk_stream << balanced_nwk(corpus.names.begin(), corpus.names.end()) << ";\n";
  nwk_stream.close();
  if (!map_stream || !fastx_stream || !nwk_stream) {
    throw std::runtime_error("failed to write the generated corpus");
  }
  return corpus;
}

/* Writes a name-to-path map. `missing` is listed but pointed at a file that does
 * not exist; `omitted` is left out of the map altogether, which is what makes
 * build_for_subtree take its "Genome skipped" branch for a backbone leaf. */
static std::filesystem::path
write_map(const Corpus& corpus, const std::filesystem::path& path, const std::string& missing, const std::string& omitted)
{
  std::ofstream map_stream(path);
  for (const std::string& name : corpus.names) {
    if (name == omitted) continue;
    const std::filesystem::path target =
      name == missing ? corpus.root / "references" / "no-such-reference.fna" : corpus.path_of(name);
    map_stream << name << "\t" << target.string() << "\n";
  }
  map_stream.close();
  if (!map_stream) throw std::runtime_error("failed to write " + path.string());
  return path;
}

/* The smallest configuration validate_configuration() accepts, so that the
 * k-mer table stays a few megabytes rather than a few hundred. */
static IndexConfig make_config(const std::filesystem::path& input, const std::filesystem::path& index_dir)
{
  IndexConfig config;
  config.input = input;
  config.index_dir = index_dir;
  config.k = 20;
  config.w = static_cast<uint8_t>(26);
  config.h = static_cast<uint8_t>(9);
  config.m = 2;
  return config;
}

static IndexConfig per_file_config(const Corpus& corpus, const std::filesystem::path& index_dir)
{
  IndexConfig config = make_config(corpus.input_map, index_dir);
  config.nwk_path = corpus.nwk_path;
  return config;
}

/* A guide tree is rejected for per-sequence indexing, so this one has none. */
static IndexConfig per_sequence_config(const Corpus& corpus, const std::filesystem::path& index_dir)
{
  return make_config(corpus.fastx, index_dir);
}

/* `disturb` runs after the input has been read but before the build, which is
 * how the per-sequence test reaches error_exit from inside the worksharing loop
 * rather than from the serial scan that precedes it. */
static void build_index_at(const IndexConfig& config, uint32_t nthreads, const std::function<void()>& disturb = {})
{
  set_num_threads(nthreads);
  gen.seed(LSH_SEED);
  IndexMultiple index(config);
  index.set_nrows();
  index.set_lshf();
  index.read_input_file();
  index.obtain_build_tree();
  if (disturb) disturb();
  index.build_index();
  index.save_index();
}

/* The index file names carry an -m<m>r<r> suffix; look them up by prefix rather
 * than rebuilding that string here. Missing is a hard error, never an empty
 * string: two absent files must not compare equal and report a pass. */
static std::filesystem::path find_by_prefix(const std::filesystem::path& index_dir, const std::string& prefix)
{
  for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    std::string name = entry.path().filename().string();
    if (name.rfind(prefix, 0) == 0 && entry.path().extension().empty()) {
      return entry.path();
    }
  }
  throw std::runtime_error("no index file named \"" + prefix + "*\" under " + index_dir.string());
}

static std::string read_index_file(const std::filesystem::path& index_dir, const std::string& prefix)
{
  std::filesystem::path path = find_by_prefix(index_dir, prefix);
  std::ifstream stream(path, std::ifstream::binary);
  std::string content((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
  if (content.empty()) throw std::runtime_error("index file is empty: " + path.string());
  return content;
}

/* cmer is a k-mer count followed by that many (encoding, subset) pairs. Only the
 * encodings are read back: union_table() merges children in completion order, so
 * a threaded build can number the subsets differently without indexing different
 * k-mers, and it is the k-mers this is about. */
static vec<uint32_t> kmer_encodings_of(const std::filesystem::path& index_dir)
{
  std::filesystem::path path = find_by_prefix(index_dir, "cmer");
  std::ifstream stream(path, std::ifstream::binary);
  uint64_t nkmers = 0;
  stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  vec<uint32_t> encodings;
  for (uint64_t ix = 0; ix < nkmers && stream; ++ix) {
    cmer_t cmer;
    stream.read(reinterpret_cast<char*>(&cmer), sizeof(cmer_t));
    encodings.push_back(cmer.first);
  }
  if (!stream || encodings.empty()) {
    throw std::runtime_error("could not read k-mer encodings from " + path.string());
  }
  return encodings;
}

static std::string kmer_count_of(const std::filesystem::path& index_dir)
{
  std::filesystem::path info_path = find_by_prefix(index_dir, "metadata");
  info_path += ".txt";
  std::ifstream stream(info_path);
  std::string info((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
  const std::string key = "total_num_kmers: ";
  size_t bix = info.find(key);
  if (bix == std::string::npos) throw std::runtime_error("no total_num_kmers in " + info_path.string());
  bix += key.size();
  return info.substr(bix, info.find('\n', bix) - bix);
}

/* Asserted only when the runtime actually gave the build a team to spread over;
 * an inactive region is not evidence either way. */
static void expect_raised_in_parallel(const std::string& what)
{
#if defined(_OPENMP) && _WOPENMP == 1
  if (raised_team_size.load() > 1) {
    expect(raised_in_parallel.load(), what);
  } else {
    skip(what + " (the runtime ran the region on a one-thread team)");
  }
#else
  skip(what + " (compiled without OpenMP)");
#endif
}

/* set_num_threads() is what an embedder calls, and if it did nothing the
 * thread-count comparisons below would silently compare two serial builds. */
static void test_set_num_threads()
{
  set_num_threads(BUILD_THREADS);
  expect(num_threads == BUILD_THREADS, "set_num_threads sets the thread count an embedder asked for");
  set_num_threads(0);
  expect(num_threads == 1, "set_num_threads clamps zero, which omp_set_num_threads does not accept");
}

/* ErrorRelay on its own, without an index in the way: many threads raising at
 * once must yield exactly one exception, with its type intact, and work queued
 * behind the first failure must not run. */
static void test_error_relay_hands_back_one_exception()
{
  ErrorRelay relay;
  const int32_t niter = 256;
#if defined(_OPENMP) && _WOPENMP == 1
  #pragma omp parallel for num_threads(8) schedule(static)
#endif
  for (int32_t ix = 0; ix < niter; ++ix) {
    relay.guard([&] { throw BuildError("relay-" + std::to_string(ix), ix); });
  }
  expect(relay.check_error(), "the relay reports an error once its region has closed");
  std::string message;
  try {
    relay.rethrow_error();
  } catch (const BuildError& error) {
    message = error.what();
  }
  expect(message.rfind("relay-", 0) == 0, "the relay rethrows one of the raised exceptions with its type intact");

  /* Serial, so the drain is a fact rather than a race. */
  ErrorRelay serial_relay;
  uint32_t entered = 0;
  for (int32_t ix = 0; ix < 8; ++ix) {
    serial_relay.guard([&] {
      entered++;
      throw BuildError("drain", ix);
    });
  }
  expect(entered == 1, "guard() runs nothing once an error has been captured");
}

/* index_files(): a reference the builder cannot open. RSeq's constructor calls
 * error_exit from inside the task graph, which is the case that used to force an
 * embedder to compile without OpenMP. */
static void test_unreadable_reference_raises(const Corpus& corpus)
{
  IndexConfig config = per_file_config(corpus, corpus.root / "index-broken-map");
  config.input = write_map(corpus, corpus.root / "broken_map.tsv", corpus.names.back(), "");
  reset_raise_state();
  std::string message;
  try {
    build_index_at(config, BUILD_THREADS);
  } catch (const BuildError& error) {
    message = error.what();
  }
  expect(message.find("no-such-reference.fna") != std::string::npos,
         "index_files: an unreadable reference reaches the caller naming the file, instead of exiting");
  expect_raised_in_parallel("index_files: the error was raised from inside a parallel region");
  /* Whether the runtime gave this build a real team decides if the thread-count
   * comparisons later on mean anything. */
  build_team_size = raised_team_size.load();
}

/* index_sequences(): the same failure on the other build strategy. The input is
 * moved aside after read_input_file() has scanned it, so that the open fails in
 * fill_slice() inside the worksharing loop rather than during the serial scan. */
static void test_unreadable_fastx_raises(const Corpus& corpus)
{
  std::filesystem::path staged = corpus.root / "staged.fna";
  std::filesystem::path hidden = corpus.root / "staged-moved-aside.fna";
  std::filesystem::copy_file(corpus.fastx, staged, std::filesystem::copy_options::overwrite_existing);
  IndexConfig config = per_sequence_config(corpus, corpus.root / "index-broken-fastx");
  config.input = staged;
  reset_raise_state();
  std::string message;
  try {
    build_index_at(config, BUILD_THREADS, [&] { std::filesystem::rename(staged, hidden); });
  } catch (const BuildError& error) {
    message = error.what();
  }
  expect(message.find("staged.fna") != std::string::npos,
         "index_sequences: an unreadable input reaches the caller naming the file, instead of exiting");
  expect_raised_in_parallel("index_sequences: the error was raised from inside a parallel region");
}

/* set_error_handler takes any callable, so the relay's catch(...) has to hand
 * back whatever was thrown. The failure has to be one that happens inside the
 * task graph: an input the serial scan rejects would never reach the relay. */
static void test_handler_may_throw_any_type(const Corpus& corpus)
{
  set_error_handler([](const std::string&, int code) {
    record_raise();
    throw code;
  });
  IndexConfig config = per_file_config(corpus, corpus.root / "index-any-type");
  config.input = corpus.root / "broken_map.tsv";
  reset_raise_state();
  bool caught_int = false;
  try {
    build_index_at(config, BUILD_THREADS);
  } catch (int) {
    caught_int = true;
  } catch (...) {
  }
  install_throwing_handler();
  expect(caught_int, "an exception that does not derive from std::exception is handed back unchanged");
  expect_raised_in_parallel("the any-type exception was raised from inside a parallel region");
}

/* Stands in for an embedder that has redirected krepp's chatter into its own log
 * and had that log fail. It throws on the first write containing `trigger` and
 * passes everything else through, so each progress critical can be singled out -
 * throwing on the first write of any kind would only ever reach the first one,
 * because std::cerr goes bad and stops calling the sink. */
class ThrowingSink : public std::streambuf
{
public:
  explicit ThrowingSink(std::string trigger)
    : trigger(std::move(trigger))
  {}
  bool fired = false;

protected:
  std::streamsize xsputn(const char* s, std::streamsize n) override
  {
    if (!fired && std::string(s, static_cast<size_t>(n)).find(trigger) != std::string::npos) {
      fired = true;
      throw BuildError("embedder log sink failed", 7);
    }
    return n;
  }
  int_type overflow(int_type c) override { return c; }

private:
  std::string trigger;
};

static void expect_sink_raise(const IndexConfig& config, const std::string& trigger, const std::string& what)
{
  ThrowingSink sink(trigger);
  std::streambuf* saved_buffer = std::cerr.rdbuf();
  std::ios_base::iostate saved_mask = std::cerr.exceptions();
  std::string message;
  std::cerr.rdbuf(&sink);
  std::cerr.exceptions(std::ios_base::badbit);
  try {
    build_index_at(config, BUILD_THREADS);
  } catch (const BuildError& error) {
    message = error.what();
  } catch (...) {
    message = "<an exception that was not a BuildError>";
  }
  /* rdbuf() first: it clears the pending badbit, which exceptions() would
   * otherwise rethrow while the sink is still installed and about to die. */
  std::cerr.rdbuf(saved_buffer);
  std::cerr.exceptions(saved_mask);
  expect(sink.fired && message.find("log sink") != std::string::npos, what);
}

/* A relay that kept the first failure would make every later build skip its work
 * and rethrow the stale exception, which is a trap for an embedder that retries. */
static void test_a_failed_build_can_be_retried(const Corpus& corpus)
{
  const std::string missing = corpus.names.back();
  std::filesystem::path late = corpus.root / "references" / "late-arrival.fna";
  std::ofstream map_stream(corpus.root / "retry_map.tsv");
  for (const std::string& name : corpus.names) {
    map_stream << name << "\t" << (name == missing ? late : corpus.path_of(name)).string() << "\n";
  }
  map_stream.close();
  std::error_code ignored;
  std::filesystem::remove(late, ignored);

  IndexConfig config = per_file_config(corpus, corpus.root / "index-retry");
  config.input = corpus.root / "retry_map.tsv";
  set_num_threads(BUILD_THREADS);
  gen.seed(LSH_SEED);
  IndexMultiple index(config);
  index.set_nrows();
  index.set_lshf();
  index.read_input_file();
  index.obtain_build_tree();

  bool first_raised = false;
  try {
    index.build_index();
  } catch (const BuildError&) {
    first_raised = true;
  }
  expect(first_raised, "the first build of a corpus with a missing reference raises");

  std::filesystem::copy_file(corpus.path_of(missing), late);
  std::string message;
  try {
    index.build_index();
    index.save_index();
  } catch (const std::exception& error) {
    message = error.what();
  }
  if (!message.empty()) std::cerr << "[INFO] retry failed with: " << message << std::endl;
  expect(message.empty(), "the same IndexMultiple builds successfully once the reference exists");
  if (!message.empty()) return;

  /* A retry that silently produced a different index would otherwise pass. */
  IndexConfig fresh = per_file_config(corpus, corpus.root / "index-retry-fresh");
  fresh.input = corpus.root / "retry_map.tsv";
  build_index_at(fresh, BUILD_THREADS);
  expect(kmer_encodings_of(config.index_dir) == kmer_encodings_of(fresh.index_dir) &&
           read_index_file(config.index_dir, "inc") == read_index_file(fresh.index_dir, "inc"),
         "the retried index holds the same k-mers as one built from scratch");
}

/* Whatever the thread count, the same k-mers have to end up in the same buckets.
 * The colour arrays (cmer's subset column and crecord) are deliberately not
 * compared: union_table() merges children in completion order, so a threaded
 * build numbers the subsets differently. That changes those bytes without
 * changing which references a k-mer resolves to. */
static void
test_thread_count_does_not_change_kmers(const IndexConfig& serial, const IndexConfig& threaded, const std::string& label)
{
  std::string message;
  try {
    build_index_at(serial, 1);
    build_index_at(threaded, BUILD_THREADS);
  } catch (const std::exception& error) {
    message = error.what();
  }
  if (!message.empty()) std::cerr << "[INFO] unexpected error: " << message << std::endl;
  expect(message.empty(), label + ": a well-formed corpus builds on one thread and on four");
  if (!message.empty()) return;
  expect(read_index_file(serial.index_dir, "inc") == read_index_file(threaded.index_dir, "inc"),
         label + ": the bucket offset array is identical across thread counts");
  expect(kmer_encodings_of(serial.index_dir) == kmer_encodings_of(threaded.index_dir),
         label + ": the indexed k-mers themselves are identical across thread counts");
  expect(kmer_count_of(serial.index_dir) == kmer_count_of(threaded.index_dir),
         label + ": the same number of k-mers is indexed on one thread and on four");
}

int main()
{
  install_throwing_handler();
  /* Named for this process: two concurrent runs must not delete each other's
   * tree and leave behind something that looks like the bug under test. */
  std::filesystem::path root =
    std::filesystem::temp_directory_path() / ("krepp-omp-index-test-" + std::to_string(static_cast<long>(getpid())));
  std::error_code ignored;
  std::filesystem::remove_all(root, ignored);
  try {
    std::filesystem::create_directories(root);
    Corpus corpus = write_corpus(root);

    test_set_num_threads();
    test_error_relay_hands_back_one_exception();
    test_unreadable_reference_raises(corpus);
    test_unreadable_fastx_raises(corpus);
    test_handler_may_throw_any_type(corpus);

    /* One per progress critical in the build path: two in build_for_subtree's
     * leaf branches, one in its internal-node branch, and two more in the
     * fill_slice() loop that index_sequences() drives. */
    IndexConfig skipped = per_file_config(corpus, root / "index-skipped-leaf");
    skipped.input = write_map(corpus, root / "skipped_map.tsv", "", corpus.names.back());
    expect_sink_raise(per_file_config(corpus, root / "index-sink-leaf"),
                      "Leaf node:",
                      "index_files: a throwing log sink at a leaf report raises instead of terminating");
    expect_sink_raise(per_file_config(corpus, root / "index-sink-internal"),
                      "Internal node:",
                      "index_files: a throwing log sink at an internal report raises instead of terminating");
    expect_sink_raise(
      skipped, "Genome skipped:", "index_files: a throwing log sink at a skipped genome raises instead of terminating");
    expect_sink_raise(per_sequence_config(corpus, root / "index-sink-seq-leaf"),
                      "Leaf node:",
                      "index_sequences: a throwing log sink at a leaf report raises instead of terminating");
    expect_sink_raise(per_sequence_config(corpus, root / "index-sink-seq-internal"),
                      "Internal node:",
                      "index_sequences: a throwing log sink at an internal report raises instead of terminating");

    /* Runs after the failed builds on purpose: a raise must leave the process
     * able to build again rather than in a state only exit() could clean up. */
    test_a_failed_build_can_be_retried(corpus);

    if (build_team_size > 1) {
      test_thread_count_does_not_change_kmers(
        per_file_config(corpus, root / "index-files-t1"), per_file_config(corpus, root / "index-files-t4"), "index_files");
      test_thread_count_does_not_change_kmers(per_sequence_config(corpus, root / "index-seqs-t1"),
                                              per_sequence_config(corpus, root / "index-seqs-t4"),
                                              "index_sequences");
    } else {
      /* Comparing two builds that both ran on one thread would say nothing about
       * thread counts, so say that instead of printing four passes. */
      skip("thread-count invariance (the build ran on a one-thread team; nothing to compare)");
      build_index_at(per_file_config(corpus, root / "index-files-serial"), 1);
      build_index_at(per_sequence_config(corpus, root / "index-seqs-serial"), 1);
      expect(!kmer_encodings_of(root / "index-files-serial").empty() &&
               !kmer_encodings_of(root / "index-seqs-serial").empty(),
             "both build strategies still produce an index holding k-mers");
    }
  } catch (const std::exception& error) {
    std::cerr << "[ERROR] the test itself could not run: " << error.what() << std::endl;
    num_failed++;
  }
  std::filesystem::remove_all(root, ignored);
  std::cerr << std::endl
            << (num_failed ? "FAILED: " : "PASSED: ") << num_checked << " checks, " << num_failed << " failing" << std::endl;
  return num_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
