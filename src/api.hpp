#ifndef _API_H
#define _API_H

// Embedding facade: the same work the `dist` and `place` subcommands do, but
// returning plain data structures instead of writing streams. The CLI
// (krepp.cpp) is unaffected and keeps driving IBatch directly; this adds a
// second driver over the same components, plus the TargetIndex method bodies
// (moved out of krepp.cpp so the module does not link the CLI layer).
//
// Error handling: Session installs a throwing error handler (scoped) and
// routes the parallel load/query regions through ErrorRelay, so corrupt
// indexes, bad FASTX files and query-time record errors surface as C++
// exceptions (RuntimeError in Python) instead of exiting the host process.
// Invalid paths/option values throw std::invalid_argument (ValueError).
// Session entry points are serialised by a process-wide mutex because the
// error handler and thread/mmap settings are globals.

#include "krepp.hpp"

namespace api {

  struct Options
  {
    uint32_t hdist_th = 4;                                      // --hdist-th [0,31]
    double chisq_value = 2.706;                                 // --chisq
    double dist_max = std::numeric_limits<double>::quiet_NaN(); // --dist-max [1e-8, 0.33]
    uint32_t tau = 2;                                           // --tau
    bool multi = true;                                          // --multi
    bool p_value = false;                                       // dist --p-value
    bool dist_filter = false;                                   // dist --filter (off by default)
    bool place_filter = true;                                   // place --filter (on by default)
    uint32_t num_threads = 1;                                   // --num-threads
    bool verbose = false;                                       // --verbose
    bool use_mmap = true;                                       // --mmap / --no-mmap
    std::filesystem::path lineage_path;                         // place -l
    std::filesystem::path nwk_path;                             // place -t
  };

  class Session : public TargetIndex
  {
  public:
    // Loads the index directory (the -i of dist/place). Throws
    // std::invalid_argument for missing paths or out-of-range options.
    Session(const std::filesystem::path& index_dir_path, const Options& options = Options{});

    // Structured equivalents of `krepp dist` / `krepp place`, plus a per-query
    // k-mer match summary. All take a FASTA/FASTQ path (gzip and URLs allowed,
    // same as the CLI's -q). With num_threads > 1 the result order across
    // batches is nondeterministic, exactly like the CLI's output order.
    vec<query_dist_t> dist(const std::string& query_path);
    vec<query_place_t> place(const std::string& query_path);
    vec<query_match_t> matches(const std::string& query_path);

    uint64_t get_total_qseq() const { return total_qseq_count; }
    const Options& get_options() const { return opts; }

  private:
    // Mirrors the CLI driver loop (omp single reader + per-batch tasks); the
    // collector runs inside the task, results are appended under a critical
    // section, so the heavy work stays parallel.
    template<typename T, typename BatchFn>
    vec<T> drive(const std::string& query_path, bool no_filter, bool p_value_flag, BatchFn&& batch_fn);

    void prepare_placement_tree(); // read_lineages() or ensure_backbone(), once
    static void check_query_path(const std::string& query_path);

    Options opts;
    uint64_t total_qseq_count = 0;
    bool tree_ready = false;
  };

} // namespace api

#endif
