#ifndef KREPP_TEST_HELPERS_HPP
#define KREPP_TEST_HELPERS_HPP

/* Shared fixtures: temp dirs, corpora, an index builder/loader and query
 * runners. Single threaded, with the global LSH generator seeded explicitly. */

#include "doctest/doctest.h"

#include "common.hpp"
#include "hdhistllh.hpp"
#include "hyperloglog.hpp"
#include "index.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "seek.hpp"
#include "sketch.hpp"
#include "table.hpp"

#include <atomic>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace ktest {

/* ------------------------------------------------------------------ paths */

inline const std::filesystem::path& temp_root()
{
  static const std::filesystem::path root = [] {
    std::filesystem::path p = std::filesystem::temp_directory_path() / "krepp-tests";
    std::filesystem::create_directories(p);
    return p;
  }();
  return root;
}

/* Removed recursively when it dies. */
class TempDir
{
public:
  explicit TempDir(const std::string& tag = "t")
    : path_(temp_root() / (tag + "-" + std::to_string(++counter_)))
  {
    std::filesystem::remove_all(path_);
    std::filesystem::create_directories(path_);
  }
  ~TempDir()
  {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
  }
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;

  const std::filesystem::path& path() const { return path_; }
  std::filesystem::path operator/(const std::string& leaf) const { return path_ / leaf; }
  std::string str() const { return path_.string(); }

private:
  std::filesystem::path path_;
  inline static std::atomic<uint32_t> counter_{0};
};

/* ------------------------------------------------------------- file utils */

inline void spit(const std::filesystem::path& p, const std::string& contents)
{
  std::filesystem::create_directories(p.parent_path());
  std::ofstream out(p, std::ios::binary);
  REQUIRE_MESSAGE(out.good(), "cannot write " << p.string());
  out << contents;
}

inline std::string slurp(const std::filesystem::path& p)
{
  std::ifstream in(p, std::ios::binary);
  REQUIRE_MESSAGE(in.good(), "cannot read " << p.string());
  return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

inline bool file_exists(const std::filesystem::path& p) { return std::filesystem::exists(p); }

/* gzip-compressed input, as real query files usually are. */
inline void write_gzip(const std::filesystem::path& p, const std::string& data)
{
  gzFile f = gzopen(p.string().c_str(), "wb");
  REQUIRE_MESSAGE(f != nullptr, "cannot open " << p.string() << " for gzip writing");
  const int written = gzwrite(f, data.data(), static_cast<unsigned>(data.size()));
  CHECK(written == static_cast<int>(data.size()));
  gzclose(f);
}

inline std::vector<std::string> read_lines(const std::filesystem::path& p)
{
  std::vector<std::string> lines;
  std::istringstream in(slurp(p));
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    lines.push_back(line);
  }
  return lines;
}

/* --------------------------------------------------------- stream capture */

/* Redirects an ostream into a buffer for the object's lifetime. */
class CaptureStream
{
public:
  explicit CaptureStream(std::ostream& os)
    : os_(os)
    , old_(os.rdbuf(buffer_.rdbuf()))
  {
  }
  ~CaptureStream() { os_.rdbuf(old_); }
  CaptureStream(const CaptureStream&) = delete;
  CaptureStream& operator=(const CaptureStream&) = delete;

  std::string str() const { return buffer_.str(); }
  std::stringstream& buffer() { return buffer_; }

private:
  std::ostream& os_;
  std::stringstream buffer_;
  std::streambuf* old_;
};

/* Index builds chatter on stderr. */
struct SilenceStderr
{
  CaptureStream cap{std::cerr};
};

/* stderr of the last build_index() call. */
inline std::string& last_build_log()
{
  static std::string log;
  return log;
}

/* ------------------------------------------------------------- sequences */

/* Deterministic DNA; `bad_every` injects an 'N' every n bases. */
inline std::string rand_dna(size_t n, uint64_t seed, uint32_t bad_every = 0)
{
  static const char alphabet[4] = {'A', 'C', 'G', 'T'};
  std::mt19937_64 rng(seed);
  std::string s;
  s.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (bad_every && ((i + 1) % bad_every == 0)) {
      s.push_back('N');
    } else {
      s.push_back(alphabet[rng() & 3]);
    }
  }
  return s;
}

/* Substitutions at a fixed rate. */
inline std::string mutate(const std::string& s, double rate, uint64_t seed)
{
  static const char alphabet[4] = {'A', 'C', 'G', 'T'};
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> pick(0.0, 1.0);
  std::string out = s;
  for (char& c : out) {
    if (pick(rng) < rate) {
      const char alt = alphabet[rng() & 3];
      if (alt != c) c = alt;
    }
  }
  return out;
}

inline std::string dna_revcomp(const std::string& s)
{
  std::string out;
  out.reserve(s.size());
  for (auto it = s.rbegin(); it != s.rend(); ++it) {
    switch (*it) {
      case 'A': out.push_back('T'); break;
      case 'C': out.push_back('G'); break;
      case 'G': out.push_back('C'); break;
      case 'T': out.push_back('A'); break;
      default: out.push_back(*it); break;
    }
  }
  return out;
}

inline std::string wrap(const std::string& seq, size_t width = 60)
{
  std::string out;
  for (size_t i = 0; i < seq.size(); i += width) {
    out += seq.substr(i, width);
    out.push_back('\n');
  }
  return out;
}

inline void write_fasta(const std::filesystem::path& p, const std::string& name, const std::string& seq)
{
  spit(p, ">" + name + "\n" + wrap(seq));
}

inline void write_fastq(const std::filesystem::path& p, const std::vector<std::pair<std::string, std::string>>& recs)
{
  std::string out;
  for (const auto& [name, seq] : recs) {
    out += "@" + name + "\n" + seq + "\n+\n" + std::string(seq.size(), 'I') + "\n";
  }
  spit(p, out);
}

/* ------------------------------------------------------------ index build */

struct Reference
{
  std::string name;
  std::string seq;
};

/* One FASTA per reference plus the TSV map. */
inline std::filesystem::path write_reference_map(const std::filesystem::path& dir, const std::vector<Reference>& refs)
{
  std::filesystem::create_directories(dir);
  std::string map;
  for (const Reference& r : refs) {
    const std::filesystem::path p = dir / (r.name + ".fna");
    write_fasta(p, r.name, r.seq);
    map += r.name + "\t" + p.string() + "\n";
  }
  const std::filesystem::path map_path = dir / "input_map.tsv";
  spit(map_path, map);
  return map_path;
}

/* Mirrors the order main() drives IndexMultiple in. */
struct BuildOptions
{
  uint8_t k = 21;
  uint8_t w = 25;
  uint8_t h = 9;
  uint32_t m = 4;
  uint32_t r = 1;
  bool frac = true;
  uint32_t sdust_t = 0;
  uint32_t sdust_w = 0;
  uint32_t seed = 42;
  std::filesystem::path nwk_path; // empty: generate a balanced tree
};

inline std::filesystem::path
build_index(const std::filesystem::path& index_dir, const std::filesystem::path& input, const BuildOptions& opt = {})
{
  IndexConfig config;
  config.input = input;
  config.index_dir = index_dir;
  config.nwk_path = opt.nwk_path;
  config.k = opt.k;
  config.w = opt.w;
  config.h = opt.h;
  config.m = opt.m;
  config.r = opt.r;
  config.frac = opt.frac;
  config.sdust_t = opt.sdust_t;
  config.sdust_w = opt.sdust_w;

  const uint32_t old_threads = num_threads;
  set_num_threads(1); // index builds must be reproducible in the tests
  gen.seed(opt.seed);
  {
    SilenceStderr quiet;
    IndexMultiple builder(config);
    builder.set_nrows();
    builder.set_lshf();
    builder.read_input_file();
    builder.obtain_build_tree();
    builder.build_index();
    builder.save_index();
    last_build_log() = quiet.cap.str();
  }
  set_num_threads(old_threads);
  return index_dir;
}

/* Convenience: build from in-memory references. */
inline std::filesystem::path build_index_from_refs(const std::filesystem::path& index_dir,
                                                   const std::vector<Reference>& refs,
                                                   const BuildOptions& opt = {})
{
  const std::filesystem::path input = write_reference_map(index_dir / "input", refs);
  return build_index(index_dir, input, opt);
}

/* The suffix IndexMultiple derives from the LSH configuration. */
inline std::string index_suffix(uint32_t m, uint32_t r, bool frac)
{
  return "-m" + std::to_string(m) + "r" + std::to_string(r) + (frac ? "-frac" : "-no_frac");
}

/* Overwrite a uint32 field of the metadata layout
 * k(1) w(1) h(1) m(4) r(4) frac(1) nrows(4) ppos(h) npos(k-h). */
inline void patch_metadata_field(const std::filesystem::path& metadata_path, size_t offset, uint32_t value)
{
  std::fstream file(metadata_path, std::ios::in | std::ios::out | std::ios::binary);
  REQUIRE_MESSAGE(file.good(), "cannot open " << metadata_path.string());
  file.seekp(static_cast<std::streamoff>(offset));
  file.write(reinterpret_cast<const char*>(&value), sizeof(value));
  file.close();
  CHECK(file.good());
}

inline constexpr size_t kMetadataOffsetM = 3;
inline constexpr size_t kMetadataOffsetR = 7;

/* Offset of the configuration block in a sketch file. */
inline size_t sketch_config_offset(const std::filesystem::path& sketch_path)
{
  const std::string data = slurp(sketch_path);
  REQUIRE(data.size() > 16);
  uint64_t nkmers = 0;
  std::memcpy(&nkmers, data.data(), sizeof(nkmers));
  uint32_t nrows = 0;
  std::memcpy(&nrows, data.data() + 8 + nkmers * sizeof(enc_t), sizeof(nrows));
  return 8 + nkmers * sizeof(enc_t) + sizeof(nrows) + static_cast<size_t>(nrows) * sizeof(inc_t);
}

/* ------------------------------------------------------------- index load */

/* Mirrors TargetIndex::load_index(). */
inline index_sptr_t load_index_dir(const std::filesystem::path& index_dir)
{
  const std::set<std::string> lall{"cmer", "crecord", "inc", "metadata", "tree", "reflist"};
  const std::set<std::string> lall_wbackbone{"cmer", "crecord", "inc", "metadata", "tree"};
  const std::set<std::string> lall_wobackbone{"cmer", "crecord", "inc", "metadata", "reflist"};
  node_phmap<std::string, std::set<std::string>> suffix_to_ltype;
  for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    const std::string filename = entry.path().filename();
    const size_t pos1 = filename.find("-", 0);
    if (pos1 == std::string::npos) continue;
    const size_t pos2 = filename.find("-", pos1 + 1);
    const std::string ltype = filename.substr(0, pos1);
    if (lall.find(ltype) == lall.end()) continue;
    if (!entry.path().extension().empty()) continue;
    const std::string suffix = filename.substr(pos1, pos2 - pos1) + filename.substr(pos2);
    suffix_to_ltype[suffix].insert(ltype);
  }
  auto index = std::make_shared<Index>(index_dir);
  for (const auto& [suffix, ltypes] : suffix_to_ltype) {
    const bool wobackbone = std::includes(ltypes.begin(), ltypes.end(), lall_wobackbone.begin(), lall_wobackbone.end());
    const bool wbackbone = std::includes(ltypes.begin(), ltypes.end(), lall_wbackbone.begin(), lall_wbackbone.end());
    if (wbackbone) {
      index->load_partial_tree(suffix);
      index->load_partial_index(suffix);
    } else if (wobackbone) {
      index->generate_partial_tree(suffix);
      index->load_partial_index(suffix);
    } else {
      error_exit("There is a partial index with a missing file!");
    }
  }
  index->make_rho_partial();
  return index;
}

/* ------------------------------------------------------------- load modes */

/* Selects the array load path (file mapping or reading) for a scope. */
class UseMmap
{
public:
  explicit UseMmap(bool enabled)
    : previous(use_mmap)
  {
    use_mmap = enabled;
  }
  ~UseMmap() { use_mmap = previous; }
  UseMmap(const UseMmap&) = delete;
  UseMmap& operator=(const UseMmap&) = delete;

private:
  bool previous;
};

/* ---------------------------------------------------------- error handling */

/* Turns error_exit() into an exception so failure paths are testable. */
class ThrowingErrorHandler
{
public:
  struct Error : public std::runtime_error
  {
    Error(const std::string& msg, int code)
      : std::runtime_error(msg)
      , code(code)
    {
    }
    int code;
  };

  ThrowingErrorHandler() { set_error_handler([](const std::string& msg, int code) { throw Error(msg, code); }); }
  ~ThrowingErrorHandler() { set_error_handler(error_handler_t()); }
  ThrowingErrorHandler(const ThrowingErrorHandler&) = delete;
  ThrowingErrorHandler& operator=(const ThrowingErrorHandler&) = delete;

  /* The message fn failed with, or empty. */
  template<typename Fn>
  static std::string catches(Fn&& fn)
  {
    try {
      fn();
    } catch (const Error& e) {
      return std::string(e.what());
    }
    return std::string();
  }
};

/* Reaches the protected BaseLSH configuration. */
class TestLSH : public BaseLSH
{
public:
  void configure(uint8_t k_, uint8_t w_, uint8_t h_, uint32_t m_, uint32_t r_, bool frac_)
  {
    k = k_;
    w = w_;
    h = h_;
    m = m_;
    r = r_;
    frac = frac_;
    sdust_t = 0;
    sdust_w = 0;
  }
  void build()
  {
    set_nrows();
    set_lshf();
  }
  lshf_sptr_t lsh() const { return lshf; }
  uint8_t win() const { return w; }
  uint8_t kmer() const { return k; }
  uint8_t npos() const { return h; }
  uint32_t rows() const { return nrows; }
};

/* ------------------------------------------------------------------ query */

/* Runs all queries through the distance path, exactly like
 * QueryIndex::estimate_distances() does for a single-threaded run. */
inline std::string dist_queries(index_sptr_t index,
                                const std::string& query_path,
                                uint32_t hdist_th = 4,
                                double chisq_value = 2.706,
                                double dist_max = std::numeric_limits<double>::quiet_NaN(),
                                uint32_t tau = 2,
                                bool no_filter = true,
                                bool multi = true,
                                bool summarize = false)
{
  auto qs = std::make_shared<QSeq>(query_path);
  strstream out;
  out.precision(STRSTREAM_PRECISION);
  out << std::fixed;
  while (qs->read_next_batch() || !qs->is_batch_finished()) {
    IBatch ib(index, qs, hdist_th, chisq_value, dist_max, tau, no_filter, multi, summarize);
    strstream batch;
    ib.estimate_distances(batch);
    out << batch.rdbuf();
  }
  return out.str();
}

/* Runs all queries through the placement path, like QueryIndex::place_sequences()
 * for a single-threaded run (without the jplace wrapper, which is CLI-side). */
inline std::string place_queries(index_sptr_t index,
                                 const std::string& query_path,
                                 bool tabular = true,
                                 uint32_t hdist_th = 4,
                                 double chisq_value = 2.706,
                                 uint32_t tau = 2,
                                 bool no_filter = false,
                                 bool multi = true,
                                 bool summarize = false)
{
  auto qs = std::make_shared<QSeq>(query_path);
  strstream out;
  out.precision(STRSTREAM_PRECISION);
  out << std::fixed;
  while (qs->read_next_batch() || !qs->is_batch_finished()) {
    IBatch ib(
      index, qs, hdist_th, chisq_value, std::numeric_limits<double>::quiet_NaN(), tau, no_filter, multi, summarize);
    strstream batch;
    ib.place_sequences(batch, tabular);
    out << batch.rdbuf();
  }
  return out.str();
}

/* Splits TSV text into data lines (comments dropped) so comparisons do not
 * depend on the order the hash maps happen to iterate in. */
inline std::vector<std::string> sorted_lines(const std::string& text, bool skip_comments = true)
{
  std::vector<std::string> lines;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (skip_comments && line[0] == '#') continue;
    lines.push_back(line);
  }
  std::sort(lines.begin(), lines.end());
  return lines;
}

/* ---------------------------------------------------------------- sketch */

/* Mirrors SketchSingle::create_sketch()/save_sketch(). */
inline std::filesystem::path
write_sketch(const std::filesystem::path& sketch_path, const std::string& input, const BuildOptions& opt)
{
  TestLSH lsh;
  lsh.configure(opt.k, opt.w, opt.h, opt.m, opt.r, opt.frac);
  gen.seed(opt.seed);
  lsh.build();
  auto rs = std::make_shared<RSeq>(input, lsh.lsh(), lsh.win(), opt.r, opt.frac, opt.sdust_t, opt.sdust_w);
  auto sdynht = std::make_shared<SDynHT>();
  sdynht->fill_table(lsh.rows(), rs);
  auto sflatht = std::make_shared<SFlatHT>(sdynht);
  std::ofstream out(sketch_path, std::ofstream::binary);
  sflatht->save(out);
  lsh.save_configuration(out);
  const double rho = rs->get_rho();
  out.write(reinterpret_cast<const char*>(&rho), sizeof(double));
  out.close();
  CHECK(out.good());
  return sketch_path;
}

/* ------------------------------------------------------------------ trees */

inline tree_sptr_t make_tree(const std::vector<std::string>& names)
{
  auto tree = std::make_shared<Tree>();
  std::vector<std::string> mutable_names = names;
  tree->generate_tree(mutable_names);
  tree->reset_traversal();
  return tree;
}

inline node_sptr_t find_leaf(tree_sptr_t tree, const std::string& name)
{
  tree->reset_traversal();
  node_sptr_t nd;
  while ((nd = tree->next_post_order())) {
    if (nd->check_leaf() && nd->get_name() == name) return nd;
  }
  return nullptr;
}

inline void collect_leaves(node_sptr_t nd, std::vector<std::string>& out)
{
  if (nd->check_leaf()) {
    out.push_back(nd->get_name());
    return;
  }
  for (tuint_t i = 0; i < nd->get_nchildren(); ++i) collect_leaves(*std::next(nd->get_children(), i), out);
}

/* Every leaf under a node, sorted and deduplicated. */
inline std::vector<std::string> clade_leaves(node_sptr_t nd)
{
  std::vector<std::string> names;
  collect_leaves(nd, names);
  std::sort(names.begin(), names.end());
  names.erase(std::unique(names.begin(), names.end()), names.end());
  return names;
}

/* A colour is stored as a set of clade roots, so the leaves it covers are the
 * leaves below those roots. */
inline std::vector<std::string> color_leaves(crecord_sptr_t crecord, se_t se)
{
  vec<node_sptr_t> subset_v;
  crecord->decode_se(se, subset_v);
  std::vector<std::string> names;
  for (const node_sptr_t& nd : subset_v) collect_leaves(nd, names);
  std::sort(names.begin(), names.end());
  names.erase(std::unique(names.begin(), names.end()), names.end());
  return names;
}

inline std::vector<std::string> color_leaves_sh(record_sptr_t record, sh_t sh)
{
  vec<node_sptr_t> subset_v;
  record->decode_sh(sh, subset_v);
  std::vector<std::string> names;
  for (const node_sptr_t& nd : subset_v) collect_leaves(nd, names);
  std::sort(names.begin(), names.end());
  names.erase(std::unique(names.begin(), names.end()), names.end());
  return names;
}

/* Post-order walk of the whole tree, as a vector. */
inline std::vector<node_sptr_t> post_order(tree_sptr_t tree)
{
  std::vector<node_sptr_t> nodes;
  tree->reset_traversal();
  node_sptr_t nd;
  while ((nd = tree->next_post_order())) nodes.push_back(nd);
  return nodes;
}

} // namespace ktest

#endif
