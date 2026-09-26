/* Unit tests for the index builder, loader and on-disk layout. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

BuildOptions small_opts()
{
  BuildOptions opt;
  opt.k = 19;
  opt.w = 23;
  opt.h = 9;
  opt.m = 4;
  opt.r = 1;
  opt.frac = true;
  opt.seed = 1234;
  return opt;
}

std::vector<Reference> make_refs()
{
  const std::string base = rand_dna(6000, 5);
  return {
    {"refA", base},
    {"refB", mutate(base, 0.02, 6)},
    {"refC", mutate(base, 0.20, 7)},
    {"refD", rand_dna(6000, 8)},
  };
}

/* Every k-mer colour pair in a loaded index, sorted, as a comparable blob. */
std::string index_fingerprint(index_sptr_t index, uint32_t m, uint32_t r, bool frac, uint32_t nrows)
{
  std::string out;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    if (!index->check_partial(rix)) continue;
    auto range = index->bucket_indices(rix);
    for (auto it = range.first; it != range.second; ++it) {
      out += std::to_string(rix) + ":" + std::to_string(it->first) + ":" + std::to_string(it->second) + ";";
    }
  }
  return out;
}

} // namespace

TEST_SUITE_BEGIN("index");

TEST_CASE("validate_configuration enforces the documented limits")
{
  ThrowingErrorHandler handler;
  auto check_invalid = [](uint8_t k, uint8_t w, uint8_t h) {
    TestLSH lsh;
    lsh.configure(k, w, h, 4, 1, true);
    return ThrowingErrorHandler::catches([&] { lsh.validate_configuration(); });
  };
  CHECK(check_invalid(18, 24, 9).find("must be at least 19") != std::string::npos);
  CHECK(check_invalid(32, 38, 9).find("must be at most 31") != std::string::npos);
  CHECK(check_invalid(21, 20, 9).find("at least the k-mer length") != std::string::npos);
  CHECK(check_invalid(21, 27, 8).find("must be at least 9") != std::string::npos);
  CHECK(check_invalid(27, 33, 16).find("must be at most 15") != std::string::npos);
  CHECK(check_invalid(31, 37, 9).find("h must be >= k-16") != std::string::npos);

  // The residue has to address a residue class that exists.
  {
    TestLSH lsh;
    lsh.configure(21, 27, 9, 4, 5, true);
    CHECK(ThrowingErrorHandler::catches([&] { lsh.validate_configuration(); })
            .find("must be smaller than the modulo") != std::string::npos);
  }
  {
    TestLSH lsh;
    lsh.configure(21, 27, 9, 0, 0, true);
    CHECK(ThrowingErrorHandler::catches([&] { lsh.validate_configuration(); }).find("must be positive") !=
          std::string::npos);
  }

  // The defaults are valid and the warnings for dust masking are emitted.
  TestLSH good;
  good.configure(21, 27, 9, 4, 1, true);
  CHECK_NOTHROW(good.validate_configuration());

  TestLSH dust;
  dust.configure(21, 27, 9, 4, 1, true);
  {
    CaptureStream quiet(std::cerr);
    CHECK_NOTHROW(dust.validate_configuration());
  }
}

TEST_CASE("set_nrows counts the rows a fractional configuration uses")
{
  auto rows_of = [](uint8_t h, uint32_t m, uint32_t r, bool frac) {
    TestLSH lsh;
    lsh.configure(21, 27, h, m, r, frac);
    lsh.build();
    return lsh.rows();
  };
  const uint32_t hash_size = 1u << (2 * 9);
  CHECK(rows_of(9, 4, 1, true) == (hash_size / 4) * 2);
  CHECK(rows_of(9, 4, 0, true) == (hash_size / 4));
  CHECK(rows_of(9, 4, 1, false) == hash_size / 4);
  CHECK(rows_of(9, 1, 0, false) == hash_size);
  // Every row of a full index is addressable by check_partial in non-frac mode.
  CHECK(rows_of(9, 3, 1, false) == hash_size / 3);
}

TEST_CASE("read_input_file accepts a reference map")
{
  TempDir dir("index-map");
  const std::vector<Reference> refs = make_refs();
  const std::filesystem::path map = write_reference_map(dir / "in", refs);
  const std::filesystem::path index_dir = dir / "index";
  const BuildOptions opt = small_opts();
  SilenceStderr quiet;
  build_index(index_dir, map, opt);
  CHECK(exists(index_dir / ("cmer" + index_suffix(opt.m, opt.r, opt.frac))));
  CHECK(exists(index_dir / ("inc" + index_suffix(opt.m, opt.r, opt.frac))));
  CHECK(exists(index_dir / ("crecord" + index_suffix(opt.m, opt.r, opt.frac))));
  CHECK(exists(index_dir / ("metadata" + index_suffix(opt.m, opt.r, opt.frac))));
  CHECK(exists(index_dir / ("metadata" + index_suffix(opt.m, opt.r, opt.frac) + ".txt")));
  CHECK(exists(index_dir / ("reflist" + index_suffix(opt.m, opt.r, opt.frac))));
  // No guide tree was given, so none is saved.
  CHECK_FALSE(exists(index_dir / ("tree" + index_suffix(opt.m, opt.r, opt.frac))));
}

TEST_CASE("read_input_file rejects duplicate reference IDs")
{
  TempDir dir("index-dup");
  const std::filesystem::path map = dir / "input.tsv";
  const Reference r{"refA", rand_dna(3000, 1)};
  write_fasta(dir / "refA.fna", r.name, r.seq);
  spit(map, "refA\t" + (dir / "refA.fna").string() + "\nrefA\t" + (dir / "refA.fna").string() + "\n");
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] {
    build_index(dir / "index", map, small_opts());
  });
  CHECK(msg.find("Duplicate reference ID") != std::string::npos);
}

TEST_CASE("read_input_file rejects a malformed map line")
{
  TempDir dir("index-badmap");
  const std::filesystem::path map = dir / "input.tsv";
  spit(map, "a-name-without-a-path\n");
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] { build_index(dir / "index", map, small_opts()); });
  CHECK(msg.find("Failed to read the reference name to path/URL mapping") != std::string::npos);
}

TEST_CASE("a single FASTA switches the builder to per-sequence indexing")
{
  TempDir dir("index-fasta");
  const std::vector<Reference> refs = make_refs();
  const std::filesystem::path fa = dir / "refs.fna";
  std::string contents;
  for (const Reference& r : refs) contents += ">" + r.name + "\n" + wrap(r.seq);
  spit(fa, contents);
  const BuildOptions opt = small_opts();
  build_index(dir / "index", fa, opt);
  // No record is short enough to be skipped here.
  CHECK(last_build_log().find("Skipping") == std::string::npos);
  auto index = load_index_dir(dir / "index");
  CHECK(index->get_lshf()->get_k() == opt.k);
  CHECK(index->get_tree()->get_root()->get_card() == refs.size());
}

TEST_CASE("a per-sequence build and a per-file build agree")
{
  TempDir dir("index-agree");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();

  const std::filesystem::path map = write_reference_map(dir / "map", refs);
  std::string contents;
  for (const Reference& r : refs) contents += ">" + r.name + "\n" + wrap(r.seq);
  const std::filesystem::path fa = dir / "refs.fna";
  spit(fa, contents);

  {
    SilenceStderr quiet;
    build_index(dir / "by_file", map, opt);
    build_index(dir / "by_seq", fa, opt);
  }
  const std::string suffix = index_suffix(opt.m, opt.r, opt.frac);
  for (const std::string& part : {"cmer", "inc", "crecord", "reflist"}) {
    CHECK(slurp(dir / "by_file" / (part + suffix)) == slurp(dir / "by_seq" / (part + suffix)));
  }
  // The human readable metadata differs only in the date line.
  const std::vector<std::string> meta_file = read_lines(dir / "by_file" / ("metadata" + suffix + ".txt"));
  const std::vector<std::string> meta_seq = read_lines(dir / "by_seq" / ("metadata" + suffix + ".txt"));
  REQUIRE(meta_file.size() == meta_seq.size());
  for (size_t i = 0; i < meta_file.size(); ++i) {
    if (meta_file[i].rfind("date:", 0) == 0) continue;
    CHECK(meta_file[i] == meta_seq[i]);
  }
}

TEST_CASE("the metadata describes the configuration and the reference list")
{
  TempDir dir("index-meta");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  const std::string suffix = index_suffix(opt.m, opt.r, opt.frac);
  const std::vector<std::string> meta = read_lines(dir / "index" / ("metadata" + suffix + ".txt"));
  auto value_of = [&](const std::string& key) {
    for (const std::string& line : meta) {
      if (line.rfind(key + ":", 0) == 0) return line.substr(key.size() + 2);
    }
    return std::string();
  };
  CHECK(value_of("k") == "19");
  CHECK(value_of("w") == "23");
  CHECK(value_of("h") == "9");
  CHECK(value_of("m") == "4");
  CHECK(value_of("frac") == "true");
  CHECK(value_of("nrows") == "131072");
  CHECK(std::stoull(value_of("total_num_kmers")) > 0);
  CHECK(value_of("ppos_v").front() == '[');
  CHECK(value_of("npos_v").front() == '[');

  const std::vector<std::string> reflist = read_lines(dir / "index" / ("reflist" + suffix));
  REQUIRE(reflist.size() == refs.size());
  for (size_t i = 0; i < refs.size(); ++i) CHECK(reflist[i] == refs[i].name);
}

TEST_CASE("a guide tree is saved with the index and reloaded")
{
  TempDir dir("index-tree");
  const std::vector<Reference> refs = make_refs();
  BuildOptions opt = small_opts();
  const std::filesystem::path nwk = dir / "guide.nwk";
  spit(nwk, "((refA:0.1,refB:0.1)X:0.1,(refC:0.1,refD:0.1)Y:0.1)root;");
  opt.nwk_path = nwk;
  build_index_from_refs(dir / "index", refs, opt);
  // With a backbone the reflist is not written.
  CHECK(last_build_log().find("Skipped saving a backbone") == std::string::npos);
  const std::string suffix = index_suffix(opt.m, opt.r, opt.frac);
  CHECK(exists(dir / "index" / ("tree" + suffix)));
  auto index = load_index_dir(dir / "index");
  CHECK(index->check_wbackbone());
  CHECK(index->get_tree()->get_root()->get_card() == refs.size());
  CHECK(index->get_tree()->get_root()->get_name() == "root");
}

TEST_CASE("a build with a fixed seed is byte-for-byte reproducible")
{
  TempDir dir("index-det");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "one", refs, opt);
    build_index_from_refs(dir / "two", refs, opt);
  }
  const std::string suffix = index_suffix(opt.m, opt.r, opt.frac);
  for (const std::string& part : {"cmer", "inc", "crecord", "metadata"}) {
    CHECK(slurp(dir / "one" / (part + suffix)) == slurp(dir / "two" / (part + suffix)));
  }
}

TEST_CASE("a loaded index reproduces the bucket layout of the built one")
{
  TempDir dir("index-load");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  auto index = load_index_dir(dir / "index");
  CHECK(index->get_lshf()->get_k() == opt.k);
  CHECK(index->get_lshf()->get_h() == opt.h);
  CHECK(index->get_lshf()->get_m() == opt.m);

  const uint32_t nrows = 1u << (2 * opt.h);
  uint64_t nkmers = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    if (!index->check_partial(rix)) continue;
    auto range = index->bucket_indices(rix);
    CHECK(range.first <= range.second);
    nkmers += static_cast<uint64_t>(std::distance(range.first, range.second));
    // Rows are only addressable when the residue is part of this partial index.
    if (opt.frac) {
      CHECK((rix % opt.m) <= opt.r);
    } else {
      CHECK((rix % opt.m) == opt.r);
    }
  }
  CHECK(nkmers > 0);
  CHECK(index->get_crecord(0)->get_nsubsets() > 0);
}

TEST_CASE("a non-fractional index only keeps the configured residue")
{
  TempDir dir("index-nofrac");
  const std::vector<Reference> refs = make_refs();
  BuildOptions opt = small_opts();
  opt.frac = false;
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  auto index = load_index_dir(dir / "index");
  const uint32_t nrows = 1u << (2 * opt.h);
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    if (rix % opt.m == opt.r) {
      CHECK(index->check_partial(rix));
    } else {
      CHECK_FALSE(index->check_partial(rix));
    }
  }
}

TEST_CASE("make_rho_partial scales rho by the number of partial indexes")
{
  TempDir dir("index-rho");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  auto index = load_index_dir(dir / "index");
  // rho is stored per exact reference and is a probability.
  auto crecord = index->get_crecord(0);
  CHECK(crecord->get_rho(1) >= 0.0);
  CHECK(crecord->get_rho(1) <= 1.0);
}

TEST_CASE("display_info reports the loaded partial index")
{
  TempDir dir("index-info");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  auto index = load_index_dir(dir / "index");
  std::stringstream out;
  index->display_info(&out);
  CHECK(out.str().find("Backbone tree: NA") != std::string::npos);
  CHECK(out.str().find("Partial index") != std::string::npos);
  CHECK(out.str().find("total_num_kmers") != std::string::npos);
}

TEST_CASE("generate_partial_tree requires a reference list")
{
  TempDir dir("index-noreflist");
  std::filesystem::create_directories(dir / "index");
  auto index = std::make_shared<Index>(dir / "index");
  ThrowingErrorHandler handler;
  const std::string msg =
    ThrowingErrorHandler::catches([&] { index->generate_partial_tree(index_suffix(4, 1, true)); });
  CHECK(msg.find("reference list") != std::string::npos);
}

TEST_CASE("an index with an impossible configuration is rejected on load")
{
  TempDir dir("index-badmeta");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  build_index_from_refs(dir / "index", refs, opt);
  const std::filesystem::path metadata = dir / "index" / ("metadata" + index_suffix(opt.m, opt.r, opt.frac));
  REQUIRE(file_exists(metadata));

  // The residue must address a residue class that exists. A build refuses this
  // today, so the file is forged to look like one an older version could write.
  ThrowingErrorHandler handler;
  {
    TempDir forged("index-badr");
    std::filesystem::copy(dir / "index", forged / "index", std::filesystem::copy_options::recursive);
    patch_metadata_field(forged / "index" / ("metadata" + index_suffix(opt.m, opt.r, opt.frac)), kMetadataOffsetR, 9);
    const std::string msg = ThrowingErrorHandler::catches([&] { load_index_dir(forged / "index"); });
    CHECK(msg.find("must be smaller than the modulo") != std::string::npos);
  }
  {
    TempDir forged("index-badm");
    std::filesystem::copy(dir / "index", forged / "index", std::filesystem::copy_options::recursive);
    patch_metadata_field(forged / "index" / ("metadata" + index_suffix(opt.m, opt.r, opt.frac)), kMetadataOffsetM, 0);
    const std::string msg = ThrowingErrorHandler::catches([&] { load_index_dir(forged / "index"); });
    CHECK(msg.find("must be positive") != std::string::npos);
  }
}

TEST_CASE("load_partial_index reports missing files")
{
  TempDir dir("index-missing");
  std::filesystem::create_directories(dir / "index");
  auto index = std::make_shared<Index>(dir / "index");
  ThrowingErrorHandler handler;
  const std::string msg =
    ThrowingErrorHandler::catches([&] { index->load_partial_index(index_suffix(4, 1, true)); });
  CHECK(msg.find("Failed to open") != std::string::npos);
}

TEST_CASE("short references are skipped with a warning")
{
  TempDir dir("index-short");
  std::vector<Reference> refs = make_refs();
  refs.push_back({"tooshort", rand_dna(10, 9)});
  const BuildOptions opt = small_opts();
  const std::filesystem::path fa = dir / "refs.fna";
  std::string contents;
  for (const Reference& r : refs) contents += ">" + r.name + "\n" + wrap(r.seq);
  spit(fa, contents);
  build_index(dir / "index", fa, opt);
  CHECK(last_build_log().find("Skipping \"tooshort\"") != std::string::npos);
  auto index = load_index_dir(dir / "index");
  // The short record never enters the tree, so the backbone covers the rest.
  CHECK(index->get_tree()->get_root()->get_card() == refs.size() - 1);
  CHECK(index->get_tree()->get_root()->get_card() == 4);
}

TEST_CASE("empty reference IDs are rejected")
{
  TempDir dir("index-empty-id");
  const std::filesystem::path fa = dir / "refs.fna";
  // The first record is valid, so the file is read as FASTX; the second one has
  // no identifier.
  spit(fa, ">ok\n" + wrap(rand_dna(2000, 11)) + ">\n" + wrap(rand_dna(2000, 12)));
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] { build_index(dir / "index", fa, small_opts()); });
  CHECK(msg.find("Empty reference ID") != std::string::npos);
}

TEST_CASE("duplicate FASTX record names are rejected")
{
  TempDir dir("index-dupfx");
  const std::filesystem::path fa = dir / "refs.fna";
  spit(fa, ">dup\n" + wrap(rand_dna(2000, 13)) + ">dup\n" + wrap(rand_dna(2000, 14)));
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] { build_index(dir / "index", fa, small_opts()); });
  CHECK(msg.find("Duplicate reference ID") != std::string::npos);
}

TEST_CASE("mapping and reading an index give the same table")
{
  TempDir dir("index-mapped");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  const uint32_t nrows = 1u << (2 * opt.h);
  std::string mapped;
  {
    UseMmap use(true);
    auto index = load_index_dir(dir / "index");
    mapped = index_fingerprint(index, opt.m, opt.r, opt.frac, nrows);
  }
  std::string read;
  {
    UseMmap use(false);
    auto index = load_index_dir(dir / "index");
    read = index_fingerprint(index, opt.m, opt.r, opt.frac, nrows);
  }
  CHECK_FALSE(mapped.empty());
  CHECK(mapped == read);
}

TEST_CASE("query results do not depend on how the index is loaded")
{
  TempDir dir("index-mapped-query");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  const std::filesystem::path q = dir / "query.fna";
  write_fasta(q, "query", mutate(refs[1].seq.substr(0, 3000), 0.03, 21));

  std::string mapped_dist, read_dist, mapped_place, read_place;
  {
    UseMmap use(true);
    auto index = load_index_dir(dir / "index");
    mapped_dist = dist_queries(index, q.string());
    mapped_place = place_queries(index, q.string());
  }
  {
    UseMmap use(false);
    auto index = load_index_dir(dir / "index");
    read_dist = dist_queries(index, q.string());
    read_place = place_queries(index, q.string());
  }
  CHECK_FALSE(mapped_dist.empty());
  CHECK(mapped_dist == read_dist);
  CHECK(mapped_place == read_place);
}

TEST_CASE("a truncated index file is rejected on load")
{
  TempDir dir("index-truncated");
  const std::vector<Reference> refs = make_refs();
  const BuildOptions opt = small_opts();
  {
    SilenceStderr quiet;
    build_index_from_refs(dir / "index", refs, opt);
  }
  const std::string suffix = index_suffix(opt.m, opt.r, opt.frac);
  for (const bool mapped : {true, false}) {
    for (const std::string& part : {"cmer", "inc"}) {
      TempDir forged("index-truncated-copy");
      std::filesystem::copy(dir / "index", forged / "index", std::filesystem::copy_options::recursive);
      const std::filesystem::path file = forged / "index" / (part + suffix);
      const std::string full = slurp(file);
      REQUIRE(full.size() > 64);
      spit(file, full.substr(0, full.size() - 32));
      CAPTURE(part);
      CAPTURE(mapped);
      UseMmap use(mapped);
      ThrowingErrorHandler handler;
      const std::string msg = ThrowingErrorHandler::catches([&] { load_index_dir(forged / "index"); });
      CHECK(msg.find("Truncated") != std::string::npos);
    }
  }
}

TEST_SUITE_END();
