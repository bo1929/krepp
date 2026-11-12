#include "api.hpp"

#include <mutex>

// ---------------------------------------------------------------------------
// TargetIndex methods, moved verbatim from krepp.cpp so that the embedding
// facade links them without pulling in the CLI layer (krepp.o has main()).
// ---------------------------------------------------------------------------

void TargetIndex::read_lineages()
{
  qtree = std::make_shared<Tree>();
  std::ifstream lineage_stream(lineage_path);
  CHECK_STREAM_OR_EXIT(lineage_stream, (std::string("Error opening ") + lineage_path.string()));
  qtree->parse_lineages(lineage_stream);
  qtree->reset_traversal();
  index->get_tree()->map_to_qtree(qtree);
  lineage_stream.close();
}

void TargetIndex::ensure_backbone()
{
  if (!nwk_path.empty()) {
    qtree = std::make_shared<Tree>();
    std::ifstream tree_stream(nwk_path);
    CHECK_STREAM_OR_EXIT(tree_stream, (std::string("Error opening ") + nwk_path.string()));
    qtree->load(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to read the given backbone tree!");
    qtree->reset_traversal();
    index->get_tree()->map_to_qtree(qtree);
    tree_stream.close();
  } else if (index->check_wbackbone()) {
    qtree = index->get_tree();
  } else {
    error_exit("Given index lacks a tree and no backbone tree is provided...");
  }
}

void TargetIndex::load_index(ErrorRelay* relay)
{
  const std::set<std::string> lall{"cmer", "crecord", "inc", "metadata", "tree", "reflist"};
  const std::set<std::string> lall_wbackbone{"cmer", "crecord", "inc", "metadata", "tree"};
  const std::set<std::string> lall_wobackbone{"cmer", "crecord", "inc", "metadata", "reflist"};
  node_phmap<std::string, std::set<std::string>> suffix_to_ltype;
  for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    std::string filename, ltype, mrcfg, fracv;
    size_t pos1, pos2;
    filename = entry.path().filename();
    pos1 = filename.find("-", 0);
    pos2 = filename.find("-", pos1 + 1);
    ltype = filename.substr(0, pos1);
    if (lall.find(ltype) != lall.end()) {
      if (entry.path().extension().empty()) {
        mrcfg = filename.substr(pos1, pos2 - pos1);
        fracv = filename.substr(pos2, filename.size() - pos2);
        suffix_to_ltype[mrcfg + fracv].insert(ltype);
      }
    } else {
    }
  }
  std::vector<std::string> suffixes;
  for (auto const& [suffix, ltypes] : suffix_to_ltype) {
    suffixes.push_back(suffix);
  }
  if (suffixes.empty()) {
    error_exit("No index found in " + index_dir.string() +
               " (expected cmer/inc/crecord/metadata files, optionally with tree/reflist).");
  }
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (uint32_t lix = 0; lix < suffixes.size(); ++lix) {
    // Same body as before; under a relay the work is guarded so an
    // error_exit() with a throwing handler (embedding API) is latched instead
    // of escaping the OpenMP region (which would be std::terminate).
    auto load_one_partial = [&] {
      const std::set<std::string>& ltype = suffix_to_ltype[suffixes[lix]];
      bool wobackbone = std::includes(ltype.begin(), ltype.end(), lall_wobackbone.begin(), lall_wobackbone.end());
      bool wbackbone = std::includes(ltype.begin(), ltype.end(), lall_wbackbone.begin(), lall_wbackbone.end());
      if (wbackbone) {
        index->load_partial_tree(suffixes[lix]);
        index->load_partial_index(suffixes[lix]);
      } else if (wobackbone) {
        index->generate_partial_tree(suffixes[lix]);
        index->load_partial_index(suffixes[lix]);
      } else {
        error_exit("There is a partial index with a missing file!");
      }
    };
    if (relay) {
      relay->guard(load_one_partial);
    } else {
      load_one_partial();
    }
  }
  if (relay) {
    relay->rethrow_error(); // on the main thread, before any further work
  }
  index->make_rho_partial();
  if (verbose) {
    index->display_load_stats(std::cerr);
  }
}

// ---------------------------------------------------------------------------
// Session
// ---------------------------------------------------------------------------

namespace api {

  namespace {
    bool looks_like_url(const std::string& s)
    {
      return s.rfind("http://", 0) == 0 || s.rfind("https://", 0) == 0 || s.rfind("ftp://", 0) == 0;
    }

    // While alive, error_exit() throws std::runtime_error instead of calling
    // std::exit, so embedders get exceptions. Combined with ErrorRelay (which
    // latches throws inside OpenMP regions and rethrows them on the main
    // thread), this is what keeps a corrupt index from killing the host
    // process. The handler is process-global, hence the mutex below.
    class ScopedThrowingErrors
    {
    public:
      ScopedThrowingErrors()
      {
        set_error_handler([](const std::string& msg, int) { throw std::runtime_error(msg); });
      }
      ~ScopedThrowingErrors() { set_error_handler(error_handler_t()); }
      ScopedThrowingErrors(const ScopedThrowingErrors&) = delete;
      ScopedThrowingErrors& operator=(const ScopedThrowingErrors&) = delete;
    };

    // Serialises Session entry points: the error handler, num_threads and
    // use_mmap are process globals, so two Sessions must not run at once.
    // (OpenMP parallelism *within* one call is unaffected.)
    std::mutex& session_mutex()
    {
      static std::mutex mtx;
      return mtx;
    }
  } // namespace

  Session::Session(const std::filesystem::path& index_dir_path, const Options& options)
    : opts(options)
  {
    // Same validations the CLI applies through CLI11 checks, as exceptions.
    if (opts.hdist_th > 31) {
      throw std::invalid_argument("hdist_th must be in [0, 31]");
    }
    if (!std::isnan(opts.dist_max) && (opts.dist_max < 1e-8 || opts.dist_max > 0.33)) {
      throw std::invalid_argument("dist_max must be in [1e-8, 0.33]");
    }
    if (!std::filesystem::is_directory(index_dir_path)) {
      throw std::invalid_argument("index directory not found: " + index_dir_path.string());
    }
    if (!opts.lineage_path.empty() && !std::filesystem::is_regular_file(opts.lineage_path)) {
      throw std::invalid_argument("lineage file not found: " + opts.lineage_path.string());
    }
    if (!opts.nwk_path.empty() && !std::filesystem::is_regular_file(opts.nwk_path)) {
      throw std::invalid_argument("nwk file not found: " + opts.nwk_path.string());
    }
    std::lock_guard<std::mutex> lock(session_mutex());
    index_dir = index_dir_path;
    lineage_path = opts.lineage_path;
    nwk_path = opts.nwk_path;
    // The globals the CLI's app callback sets before any subcommand runs.
    set_num_threads(opts.num_threads);
    verbose = opts.verbose;
    use_mmap = opts.use_mmap;
    ScopedThrowingErrors handler_scope;
    ErrorRelay relay;
    index = std::make_shared<Index>(index_dir);
    load_index(&relay); // throws std::runtime_error on corrupt/missing files
  }

  void Session::check_query_path(const std::string& query_path)
  {
    if (looks_like_url(query_path)) {
      return;
    }
    if (!std::filesystem::is_regular_file(query_path)) {
      throw std::invalid_argument("query file not found: " + query_path);
    }
  }

  // Called with session_mutex() held (see place()).
  void Session::prepare_placement_tree()
  {
    if (tree_ready) {
      return;
    }
    // Mirrors QueryIndex::validate_configuration_place().
    if (opts.hdist_th < opts.tau) {
      throw std::invalid_argument("tau must be <= hdist_th");
    }
    if (lineage_path.empty() && nwk_path.empty() && !index->check_wbackbone()) {
      // Preempt ensure_backbone()'s error_exit so embedders get an exception.
      throw std::invalid_argument("index lacks a tree and no backbone tree or lineage file is provided");
    }
    // Main-thread only: a throwing handler turns CHECK_STREAM_OR_EXIT /
    // error_exit failures into exceptions directly; no relay needed.
    ScopedThrowingErrors handler_scope;
    if (!lineage_path.empty()) {
      read_lineages();
    } else {
      ensure_backbone();
    }
    tree_ready = true;
  }

  // Called with session_mutex() held by dist()/place()/matches().
  template<typename T, typename BatchFn>
  vec<T> Session::drive(const std::string& query_path, bool no_filter, bool p_value_flag, BatchFn&& batch_fn)
  {
    check_query_path(query_path);
    ScopedThrowingErrors handler_scope;
    ErrorRelay relay;
    vec<T> results;
#if defined(_OPENMP) && _WOPENMP == 1
    // This session's own setting, not the mutable global: another Session
    // constructed later must not change how this one queries.
    omp_set_num_threads(static_cast<int>(opts.num_threads));
#endif
    qseq_sptr_t qs = std::make_shared<QSeq>(query_path);
#pragma omp parallel shared(qs)
    {
#pragma omp single
      {
        // The reader loop runs inside omp single; guard it so an
        // error_exit() throw (corrupt FASTX, gz errors) is latched instead
        // of escaping the region. Once latched, stop reading and let the
        // already-spawned tasks drain as no-ops.
        relay.guard([&] {
          bool cont_reading = false;
          while (!relay.check_error() && ((cont_reading = qs->read_next_batch()) || !qs->is_batch_finished())) {
            total_qseq_count += qs->get_cbatch_size();
            auto ib = std::make_shared<IBatch>(
              index, qs, opts.hdist_th, opts.chisq_value, opts.dist_max, opts.tau, no_filter, opts.multi, false, p_value_flag);
#pragma omp task firstprivate(ib)
            {
              relay.guard([&] {
                vec<T> local = batch_fn(*ib);
#pragma omp critical(krepp_api_results)
                {
                  results.insert(
                    results.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
                }
              });
            }
          }
        });
#pragma omp taskwait
      }
    }
    relay.rethrow_error();
    return results;
  }

  vec<query_dist_t> Session::dist(const std::string& query_path)
  {
    std::lock_guard<std::mutex> lock(session_mutex());
    return drive<query_dist_t>(
      query_path, !opts.dist_filter, opts.p_value, [](IBatch& ib) { return ib.dist_batch(); });
  }

  vec<query_place_t> Session::place(const std::string& query_path)
  {
    std::lock_guard<std::mutex> lock(session_mutex());
    prepare_placement_tree();
    return drive<query_place_t>(
      query_path, !opts.place_filter, false, [](IBatch& ib) { return ib.place_batch(); });
  }

  vec<query_match_t> Session::matches(const std::string& query_path)
  {
    std::lock_guard<std::mutex> lock(session_mutex());
    // Unfiltered view of the k-mer matches: every reference that had at least
    // one k-mer within hdist_th, with its histogram and distance.
    return drive<query_match_t>(
      query_path, true, false, [](IBatch& ib) { return ib.match_batch(); });
  }

  // Explicit instantiation is not needed: drive is only used from this TU and
  // the definitions above precede all uses.

} // namespace api
