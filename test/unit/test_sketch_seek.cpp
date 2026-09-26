/* Unit tests for building, loading and seeking in a single sketch. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

BuildOptions sketch_opts()
{
  BuildOptions opt;
  opt.k = 19;
  opt.w = 23;
  opt.h = 9;
  opt.m = 4;
  opt.r = 1;
  opt.frac = true;
  opt.seed = 777;
  return opt;
}

struct SeekRow
{
  std::string id;
  double distance;
  bool is_na;
};

std::vector<SeekRow> parse_seek(const std::string& text)
{
  std::vector<SeekRow> rows;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    SeekRow row;
    std::string dist;
    if (!std::getline(ls, row.id, '\t')) continue;
    if (!std::getline(ls, dist)) continue;
    row.is_na = (dist == "NaN");
    row.distance = row.is_na ? std::numeric_limits<double>::quiet_NaN() : std::stod(dist);
    rows.push_back(row);
  }
  return rows;
}

std::string seek_queries(sketch_sptr_t sketch, const std::string& query_path, uint32_t hdist_th = 4)
{
  auto qs = std::make_shared<QSeq>(query_path);
  strstream out;
  out.precision(STRSTREAM_PRECISION);
  out << std::fixed;
  while (qs->read_next_batch() || !qs->is_batch_finished()) {
    SBatch sb(sketch, qs, hdist_th);
    sb.seek_sequences(out);
  }
  return out.str();
}

} // namespace

TEST_SUITE_BEGIN("sketch");

TEST_CASE("a sketch round-trips through save and load")
{
  TempDir dir("sketch-round");
  const std::string seq = rand_dna(30000, 201);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const BuildOptions opt = sketch_opts();
  const std::filesystem::path sketch_path = write_sketch(dir / "ref.sketch", fa.string(), opt);
  REQUIRE(file_exists(sketch_path));

  auto sketch = std::make_shared<Sketch>(sketch_path);
  sketch->load_full_sketch();
  CHECK(sketch->get_lshf()->get_k() == opt.k);
  CHECK(sketch->get_lshf()->get_h() == opt.h);
  CHECK(sketch->get_lshf()->get_m() == opt.m);
  CHECK(sketch->get_sflatht_sptr()->bucket_next(0) >= sketch->get_sflatht_sptr()->bucket_start(0));
  CHECK(sketch->get_rho() > 0.0);
  CHECK(sketch->get_rho() <= 1.0);

  // The bucket rows follow the same residue rule as an index.
  const uint32_t nrows = 1u << (2 * opt.h);
  uint64_t nkmers = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    if (!sketch->check_partial(rix)) continue;
    auto range = sketch->bucket_indices(rix);
    CHECK(range.first <= range.second);
    nkmers += static_cast<uint64_t>(std::distance(range.first, range.second));
  }
  CHECK(nkmers > 0);
}

TEST_CASE("check_partial follows the fraction rule")
{
  TempDir dir("sketch-frac");
  const std::string seq = rand_dna(20000, 202);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);

  for (bool frac : {true, false}) {
    BuildOptions opt = sketch_opts();
    opt.frac = frac;
    const std::filesystem::path path = dir / (frac ? "frac.sketch" : "nofrac.sketch");
    write_sketch(path, fa.string(), opt);
    auto sketch = std::make_shared<Sketch>(path);
    sketch->load_full_sketch();
    const uint32_t nrows = 1u << (2 * opt.h);
    for (uint32_t rix = 0; rix < nrows; rix += 97) {
      if (frac) {
        CHECK(sketch->check_partial(rix) == ((rix % opt.m) <= opt.r));
      } else {
        CHECK(sketch->check_partial(rix) == ((rix % opt.m) == opt.r));
      }
    }
  }
}

TEST_CASE("make_rho_partial scales the subsampling rate")
{
  TempDir dir("sketch-rho");
  const std::string seq = rand_dna(20000, 203);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  BuildOptions opt = sketch_opts();
  opt.m = 4;
  opt.r = 1;
  opt.frac = true;
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), opt);
  auto sketch = std::make_shared<Sketch>(path);
  sketch->load_full_sketch();
  const double before = sketch->get_rho();
  sketch->make_rho_partial();
  // Fractional sketches keep residues 0..r, i.e. (r + 1) / m of the k-mers.
  CHECK(sketch->get_rho() == doctest::Approx(before * 2.0 / 4.0));
  CHECK(sketch->get_rho() <= 1.0);

  BuildOptions no_frac = opt;
  no_frac.frac = false;
  const std::filesystem::path path2 = write_sketch(dir / "ref2.sketch", fa.string(), no_frac);
  auto sketch2 = std::make_shared<Sketch>(path2);
  sketch2->load_full_sketch();
  const double before2 = sketch2->get_rho();
  sketch2->make_rho_partial();
  CHECK(sketch2->get_rho() == doctest::Approx(before2 / 4.0));
}

TEST_CASE("a truncated sketch file is rejected")
{
  TempDir dir("sketch-trunc");
  const std::string seq = rand_dna(5000, 204);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  const std::string full = slurp(path);
  REQUIRE(full.size() > 100);
  const std::filesystem::path truncated = dir / "trunc.sketch";
  spit(truncated, full.substr(0, full.size() / 2));

  auto sketch = std::make_shared<Sketch>(truncated);
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] { sketch->load_full_sketch(); });
  CHECK_FALSE(msg.empty());
}

TEST_CASE("a sketch with an impossible configuration is rejected on load")
{
  TempDir dir("sketch-badcfg");
  const std::string seq = rand_dna(20000, 205);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const BuildOptions opt = sketch_opts();
  const std::filesystem::path good = write_sketch(dir / "ref.sketch", fa.string(), opt);
  const size_t config_offset = sketch_config_offset(good);

  ThrowingErrorHandler handler;
  const std::pair<const char*, uint32_t> cases[] = {
    {"must be smaller than the modulo", 9}, // r = 9 with m = 4
  };
  for (const auto& [needle, value] : cases) {
    const std::filesystem::path bad = dir / "bad.sketch";
    std::filesystem::copy_file(good, bad, std::filesystem::copy_options::overwrite_existing);
    patch_metadata_field(bad, config_offset + 7, value); // r
    auto sketch = std::make_shared<Sketch>(bad);
    const std::string msg = ThrowingErrorHandler::catches([&] { sketch->load_full_sketch(); });
    CHECK(msg.find(needle) != std::string::npos);
  }
  {
    const std::filesystem::path bad = dir / "bad-m.sketch";
    std::filesystem::copy_file(good, bad, std::filesystem::copy_options::overwrite_existing);
    patch_metadata_field(bad, config_offset + 3, 0); // m
    auto sketch = std::make_shared<Sketch>(bad);
    CHECK(ThrowingErrorHandler::catches([&] { sketch->load_full_sketch(); }).find("must be positive") !=
          std::string::npos);
  }
}

TEST_CASE("seek is unaffected by how the sketch is loaded")
{
  TempDir dir("sketch-mapped");
  const std::string seq = rand_dna(30000, 206);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  const std::filesystem::path q = dir / "query.fna";
  write_fasta(q, "query", mutate(seq.substr(0, 15000), 0.05, 207));

  std::string mapped, read;
  {
    UseMmap use(true);
    auto sketch = std::make_shared<Sketch>(path);
    sketch->load_full_sketch();
    mapped = seek_queries(sketch, q.string());
  }
  {
    UseMmap use(false);
    auto sketch = std::make_shared<Sketch>(path);
    sketch->load_full_sketch();
    read = seek_queries(sketch, q.string());
  }
  CHECK_FALSE(mapped.empty());
  CHECK(mapped == read);
}

TEST_CASE("a sketch truncated in its arrays is rejected on load")
{
  TempDir dir("sketch-trunc-arrays");
  const std::string seq = rand_dna(5000, 208);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  // The configuration block follows the arrays, so cut inside the arrays.
  const size_t arrays_end = sketch_config_offset(path);
  REQUIRE(arrays_end > 8);
  const std::filesystem::path truncated = dir / "trunc.sketch";
  spit(truncated, slurp(path).substr(0, arrays_end - 8));
  for (const bool mapped : {true, false}) {
    CAPTURE(mapped);
    UseMmap use(mapped);
    auto sketch = std::make_shared<Sketch>(truncated);
    ThrowingErrorHandler handler;
    const std::string msg = ThrowingErrorHandler::catches([&] { sketch->load_full_sketch(); });
    CHECK(msg.find("Truncated") != std::string::npos);
  }
}

TEST_SUITE_END();

TEST_SUITE_BEGIN("seek");

TEST_CASE("seek places an exact copy at distance zero")
{
  TempDir dir("seek-exact");
  const std::string seq = rand_dna(30000, 301);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  auto sketch = std::make_shared<Sketch>(path);
  sketch->load_full_sketch();

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq,
              {{"exact", seq.substr(1000, 500)},
               {"mutated", mutate(seq.substr(3000, 500), 0.05, 5)},
               {"unrelated", rand_dna(500, 6)}});
  const std::vector<SeekRow> rows = parse_seek(seek_queries(sketch, fq.string()));
  REQUIRE(rows.size() == 3);
  CHECK(rows[0].id == "exact");
  CHECK_FALSE(rows[0].is_na);
  // The estimate is not exactly zero: a sketch only stores the reference's
  // minimizers, so a query k-mer that is not a minimizer can still land in a
  // bucket and be counted as a near match.
  CHECK(rows[0].distance < 0.05);
  CHECK_FALSE(rows[1].is_na);
  CHECK(rows[1].distance > rows[0].distance);
  CHECK(rows[1].distance < 0.3);
  // An unrelated read must not look close, whether or not a stray k-mer
  // happens to collide.
  CHECK(rows[2].id == "unrelated");
  if (!rows[2].is_na) CHECK(rows[2].distance > 0.3);
}

TEST_CASE("seek handles the reverse complement and short queries")
{
  TempDir dir("seek-rc");
  const std::string seq = rand_dna(30000, 302);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  auto sketch = std::make_shared<Sketch>(path);
  sketch->load_full_sketch();

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"reverse", dna_revcomp(seq.substr(5000, 600))}, {"tiny", "ACGTACGTAC"}});
  const std::vector<SeekRow> rows = parse_seek(seek_queries(sketch, fq.string()));
  REQUIRE(rows.size() == 2);
  CHECK(rows[0].id == "reverse");
  CHECK_FALSE(rows[0].is_na);
  CHECK(rows[0].distance < 0.05);
  CHECK(rows[1].id == "tiny");
  CHECK(rows[1].is_na);
}

TEST_CASE("seek batches several queries")
{
  TempDir dir("seek-batch");
  const std::string seq = rand_dna(40000, 303);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  auto sketch = std::make_shared<Sketch>(path);
  sketch->load_full_sketch();

  std::vector<std::pair<std::string, std::string>> recs;
  for (uint32_t i = 0; i < 40; ++i) recs.emplace_back("read" + std::to_string(i), seq.substr(i * 500, 400));
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, recs);
  const std::vector<SeekRow> rows = parse_seek(seek_queries(sketch, fq.string()));
  REQUIRE(rows.size() == recs.size());
  for (uint32_t i = 0; i < recs.size(); ++i) {
    CHECK(rows[i].id == recs[i].first);
    CHECK_FALSE(rows[i].is_na);
    CHECK(rows[i].distance < 0.05);
  }
}

TEST_CASE("seek with hdist-th zero only reports exact k-mer matches")
{
  TempDir dir("seek-hdist0");
  const std::string seq = rand_dna(30000, 304);
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", seq);
  const std::filesystem::path path = write_sketch(dir / "ref.sketch", fa.string(), sketch_opts());
  auto sketch = std::make_shared<Sketch>(path);
  sketch->load_full_sketch();

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"exact", seq.substr(2000, 400)}, {"mutated", mutate(seq.substr(4000, 400), 0.10, 7)}});
  const std::vector<SeekRow> rows = parse_seek(seek_queries(sketch, fq.string(), 0));
  REQUIRE(rows.size() == 2);
  CHECK_FALSE(rows[0].is_na);
  // A heavily mutated read may still match some k-mer exactly, but its
  // estimate cannot be better than the exact copy's.
  if (!rows[1].is_na) CHECK(rows[1].distance >= rows[0].distance);
}

TEST_SUITE_END();
