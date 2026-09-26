/* Unit tests for RSeq (minimizer extraction) and QSeq (query batches). */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

const uint8_t K = 15;
const uint8_t W = 20;
const uint8_t H = 9;

TestLSH make_lsh(uint32_t m = 4, uint32_t r = 1, bool frac = true)
{
  TestLSH lsh;
  lsh.configure(K, W, H, m, r, frac);
  gen.seed(99);
  lsh.build();
  return lsh;
}

/* Independent re-implementation of the extraction loop, used as the oracle for
 * RSeq::extract_mers(). Returns (row, dropped encoding) in emission order. */
std::vector<std::pair<uint32_t, uint32_t>> extract_ref(const std::string& seq, TestLSH& lsh, uint32_t m, uint32_t r, bool frac)
{
  const uint8_t k = K;
  const uint32_t ldiff = (W > k) ? (W - k + 1) : 1;
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));

  std::vector<std::pair<uint32_t, uint32_t>> out;
  std::vector<std::array<uint64_t, 3>> window(ldiff, {0, 0, 0}); // bp, lr, hash
  uint64_t kix = 0;
  uint32_t l = 0;
  for (size_t i = 0; i < seq.size(); ++i) {
    if (seq_nt4_table[static_cast<unsigned char>(seq[i])] >= 4) {
      l = 0;
      continue;
    }
    ++l;
    if (l < k) continue;
    // Encode the k-mer ending at i.
    uint64_t bp = 0, lr = 0;
    for (size_t j = i + 1 - k; j <= i; ++j) {
      const uint8_t code = seq_nt4_table[static_cast<unsigned char>(seq[j])];
      bp = (bp << 2) | nt4_bp_table[code];
      lr = (lr << 1) + nt4_lr_table[code];
    }
    bp &= mask_bp;
    lr &= mask_lr;
    window[kix % ldiff] = {bp, lr, xur64_hash(bp)};
    ++kix;
    if (l < W && (i + 1) != seq.size()) continue;
    const auto min_it = std::min_element(window.begin(), window.end(), [](const std::array<uint64_t, 3>& a, const std::array<uint64_t, 3>& b) {
      return a[2] < b[2];
    });
    const uint64_t min_bp = (*min_it)[0];
    const uint64_t min_lr = (*min_it)[1];
    const uint32_t rix = lsh.lsh()->compute_hash(min_bp);
    const uint32_t rix_res = rix % m;
    if (frac ? rix_res <= r : rix_res == r) {
      const uint32_t row = frac ? rix / m * (r + 1) + rix_res : rix / m;
      out.emplace_back(row, lsh.lsh()->drop_ppos_lr(min_lr));
    }
  }
  return out;
}

std::vector<std::pair<uint32_t, uint32_t>> flatten(const vvec<mer_t>& table)
{
  std::vector<std::pair<uint32_t, uint32_t>> out;
  for (uint32_t row = 0; row < table.size(); ++row) {
    for (const mer_t& mer : table[row]) out.emplace_back(row, mer.encoding);
  }
  return out;
}

} // namespace

TEST_SUITE_BEGIN("rqseq");

TEST_CASE("RSeq walks the records of a FASTA file")
{
  TempDir dir("rseq-fa");
  const std::filesystem::path fa = dir / "refs.fna";
  spit(fa, ">first\n" + wrap(rand_dna(500, 1)) + ">second\n" + wrap(rand_dna(400, 2)));
  TestLSH lsh = make_lsh();
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "first");
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "second");
  CHECK_FALSE(rs.read_next_seq());
}

TEST_CASE("RSeq reads FASTQ as well")
{
  TempDir dir("rseq-fq");
  const std::filesystem::path fq = dir / "refs.fq";
  write_fastq(fq, {{"r1", rand_dna(300, 3)}, {"r2", rand_dna(300, 4)}});
  TestLSH lsh = make_lsh();
  RSeq rs(fq.string(), lsh.lsh(), W, 1, true, 0, 0);
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "r1");
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "r2");
  CHECK_FALSE(rs.read_next_seq());
}

TEST_CASE("RSeq reads gzipped input")
{
  TempDir dir("rseq-gz");
  const std::filesystem::path gz = dir / "refs.fna.gz";
  write_gzip(gz, ">gzref\n" + wrap(rand_dna(600, 5)));
  TestLSH lsh = make_lsh();
  RSeq rs(gz.string(), lsh.lsh(), W, 1, true, 0, 0);
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "gzref");
}

TEST_CASE("RSeq can start reading at a byte offset")
{
  TempDir dir("rseq-offset");
  const std::filesystem::path fa = dir / "refs.fna";
  const std::string rec1 = ">first\n" + wrap(rand_dna(300, 6));
  const std::string rec2 = ">second\n" + wrap(rand_dna(300, 7));
  spit(fa, rec1 + rec2);
  TestLSH lsh = make_lsh();
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0, rec1.size());
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
  CHECK(std::string(rs.get_name()) == "second");
}

TEST_CASE("RSeq refuses a missing file")
{
  ThrowingErrorHandler handler;
  TestLSH lsh = make_lsh();
  const std::string msg =
    ThrowingErrorHandler::catches([&] { RSeq rs("/nonexistent/krepp/refs.fna", lsh.lsh(), W, 1, true, 0, 0); });
  CHECK(msg.find("Failed to open the file at") != std::string::npos);
}

TEST_CASE("read_next_seq reports sequences shorter than the window")
{
  TempDir dir("rseq-short");
  const std::filesystem::path fa = dir / "refs.fna";
  spit(fa, ">tiny\nACGTACGT\n>ok\n" + wrap(rand_dna(400, 8)));
  TestLSH lsh = make_lsh();
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  CHECK(rs.read_next_seq());
  CHECK_FALSE(rs.set_curr_seq()); // shorter than w
  CHECK(rs.read_next_seq());
  CHECK(rs.set_curr_seq());
}

TEST_CASE("extract_mers matches an independent minimizer implementation")
{
  for (uint32_t m : {1u, 4u}) {
    for (bool frac : {true, false}) {
      for (uint64_t seed : {11ull, 12ull, 13ull}) {
        TempDir dir("rseq-mers");
        const std::filesystem::path fa = dir / "ref.fna";
        const std::string seq = rand_dna(900, seed);
        write_fasta(fa, "ref", seq);
        TestLSH lsh = make_lsh(m, 1, frac);
        RSeq rs(fa.string(), lsh.lsh(), W, 1, frac, 0, 0);
        REQUIRE(rs.read_next_seq());
        REQUIRE(rs.set_curr_seq());
        vvec<mer_t> table(lsh.rows());
        rs.extract_mers(table, 42);
        // The table groups entries by row, so both sides are compared as
        // multisets rather than in emission order.
        std::vector<std::pair<uint32_t, uint32_t>> got = flatten(table);
        std::vector<std::pair<uint32_t, uint32_t>> expected = extract_ref(seq, lsh, m, 1, frac);
        REQUIRE(got.size() == expected.size());
        std::sort(got.begin(), got.end());
        std::sort(expected.begin(), expected.end());
        CHECK(got == expected);
        // Every entry carries the shadow hash it was built with.
        for (const auto& row : table) {
          for (const mer_t& mer : row) CHECK(mer.sh == 42);
        }
      }
    }
  }
}

TEST_CASE("extract_mers breaks the k-mer run at ambiguous bases")
{
  TempDir dir("rseq-n");
  const std::filesystem::path fa = dir / "ref.fna";
  const std::string seq = rand_dna(600, 21, 37); // an N every 37 bases
  write_fasta(fa, "ref", seq);
  TestLSH lsh = make_lsh();
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  REQUIRE(rs.read_next_seq());
  REQUIRE(rs.set_curr_seq());
  vvec<mer_t> table(lsh.rows());
  rs.extract_mers(table, 1);
  std::vector<std::pair<uint32_t, uint32_t>> got = flatten(table);
  std::vector<std::pair<uint32_t, uint32_t>> expected = extract_ref(seq, lsh, 4, 1, true);
  std::sort(got.begin(), got.end());
  std::sort(expected.begin(), expected.end());
  CHECK(got == expected);
  CHECK_FALSE(got.empty());
}

TEST_CASE("extract_mers with w <= k falls back to single k-mers")
{
  TempDir dir("rseq-wk");
  const std::filesystem::path fa = dir / "ref.fna";
  const std::string seq = rand_dna(400, 31);
  write_fasta(fa, "ref", seq);
  TestLSH lsh;
  lsh.configure(K, K - 2, H, 4, 1, true); // w < k on purpose
  gen.seed(99);
  lsh.build();
  RSeq rs(fa.string(), lsh.lsh(), lsh.win(), 1, true, 0, 0);
  REQUIRE(rs.read_next_seq());
  REQUIRE(rs.set_curr_seq());
  vvec<mer_t> table(lsh.rows());
  rs.extract_mers(table, 1);
  CHECK_FALSE(flatten(table).empty());
}

TEST_CASE("the rho estimate is a fraction of the k-mer cardinality")
{
  TempDir dir("rseq-rho");
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", rand_dna(30000, 41));
  TestLSH lsh = make_lsh(4, 1, true);
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  REQUIRE(rs.read_next_seq());
  REQUIRE(rs.set_curr_seq());
  vvec<mer_t> table(lsh.rows());
  rs.extract_mers(table, 1);
  rs.compute_rho();
  const double rho = rs.get_rho();
  CHECK(rho > 0.0);
  CHECK(rho <= 1.0);
  // The subsampling rate is m/(r+1) for a fractional configuration, and the
  // estimate should be in the right ballpark.
  CHECK(rho == doctest::Approx(0.25).epsilon(0.35));
}

TEST_CASE("reset_estimates clears the accumulators")
{
  TempDir dir("rseq-reset");
  const std::filesystem::path fa = dir / "ref.fna";
  write_fasta(fa, "ref", rand_dna(5000, 51));
  TestLSH lsh = make_lsh();
  RSeq rs(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  REQUIRE(rs.read_next_seq());
  REQUIRE(rs.set_curr_seq());
  vvec<mer_t> table(lsh.rows());
  rs.extract_mers(table, 1);
  rs.compute_rho();
  CHECK(rs.get_rho() > 0.0);
  // A second pass accumulates on top of the first unless the estimates are
  // reset, which is what per-sequence indexing does.
  rs.reset_estimates();
  vvec<mer_t> table2(lsh.rows());
  rs.extract_mers(table2, 2);
  rs.compute_rho();
  CHECK(rs.get_rho() > 0.0);
}

TEST_CASE("dust masking removes low complexity sequence")
{
  TempDir dir("rseq-dust");
  const std::filesystem::path fa = dir / "ref.fna";
  // A long poly-A stretch followed by random sequence.
  const std::string seq = std::string(400, 'A') + rand_dna(600, 61);
  write_fasta(fa, "ref", seq);
  TestLSH lsh = make_lsh();

  RSeq plain(fa.string(), lsh.lsh(), W, 1, true, 0, 0);
  REQUIRE(plain.read_next_seq());
  REQUIRE(plain.set_curr_seq());
  vvec<mer_t> plain_table(lsh.rows());
  plain.extract_mers(plain_table, 1);

  RSeq masked(fa.string(), lsh.lsh(), W, 1, true, 20, 64);
  REQUIRE(masked.read_next_seq());
  REQUIRE(masked.set_curr_seq());
  vvec<mer_t> masked_table(lsh.rows());
  masked.extract_mers(masked_table, 1);

  CHECK(flatten(masked_table).size() < flatten(plain_table).size());
}

TEST_CASE("QSeq reads batches of queries")
{
  TempDir dir("qseq");
  const std::filesystem::path fq = dir / "reads.fq";
  std::vector<std::pair<std::string, std::string>> recs;
  for (uint32_t i = 0; i < 600; ++i) recs.emplace_back("read" + std::to_string(i), rand_dna(150, i + 1));
  write_fastq(fq, recs);

  QSeq qs(fq.string());
  // Nothing has been read yet, so there is no pending batch.
  CHECK(qs.is_batch_finished());
  uint64_t total = 0;
  uint32_t batches = 0;
  while (qs.read_next_batch() || !qs.is_batch_finished()) {
    const uint64_t n = qs.get_cbatch_size();
    // The batch byte limit is 512 * 150, so each batch holds about 512 reads.
    CHECK(n <= 512);
    total += n;
    ++batches;
    if (qs.is_batch_finished() && n == 0) break;
  }
  CHECK(total == 600);
  CHECK(batches >= 2);
}

TEST_CASE("QSeq keeps a single long record in one batch")
{
  TempDir dir("qseq-long");
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"contig", rand_dna(200000, 71)}});
  QSeq qs(fq.string());
  CHECK(qs.read_next_batch());
  CHECK(qs.get_cbatch_size() == 1);
  CHECK_FALSE(qs.read_next_batch());
  CHECK(qs.is_batch_finished());
}

TEST_CASE("QSeq refuses a missing file")
{
  ThrowingErrorHandler handler;
  CHECK(ThrowingErrorHandler::catches([] { QSeq qs("/nonexistent/krepp/reads.fq"); }).find("Failed to open the file") !=
        std::string::npos);
}

TEST_SUITE_END();
