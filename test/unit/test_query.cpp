/* Integration tests for distance estimation and placement over a small index,
 * driven through IBatch the way the CLI drives it. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

BuildOptions query_opts()
{
  BuildOptions opt;
  opt.k = 19;
  opt.w = 23;
  opt.h = 9;
  opt.m = 4;
  opt.r = 1;
  opt.frac = true;
  opt.seed = 4242;
  return opt;
}

/* Two closely related references, one distant relative and one unrelated
 * sequence, all long enough to produce a few thousand k-mers. */
struct Corpus
{
  std::string base;
  std::vector<Reference> refs;
  Reference query_source; // a piece of refA, the "true" source of the queries
};

Corpus make_corpus()
{
  Corpus c;
  c.base = rand_dna(20000, 101);
  c.refs = {
    {"refA", c.base},
    {"refB", mutate(c.base, 0.01, 102)},
    {"refC", mutate(c.base, 0.25, 103)},
    {"refD", rand_dna(20000, 104)},
  };
  return c;
}

std::string substr_of(const std::string& s, size_t start, size_t len) { return s.substr(start, len); }

/* Fields of one distance row: id, reference, distance, p-value. */
struct DistRow
{
  std::string id;
  std::string reference;
  double distance;
  double p_value;
  bool is_na;
  int nfields = 0;
};

std::vector<DistRow> parse_dist(const std::string& text)
{
  std::vector<DistRow> rows;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    DistRow row;
    std::string dist, pval;
    if (!std::getline(ls, row.id, '\t')) continue;
    if (!std::getline(ls, row.reference, '\t')) continue;
    if (!std::getline(ls, dist, '\t')) continue;
    // The P_VALUE column is optional.
    const bool has_p_value = static_cast<bool>(std::getline(ls, pval));
    row.nfields = static_cast<int>(std::count(line.begin(), line.end(), '\t')) + 1;
    row.is_na = (dist == "NaN");
    row.distance = row.is_na ? std::numeric_limits<double>::quiet_NaN() : std::stod(dist);
    row.p_value = (!has_p_value || pval == "NaN") ? std::numeric_limits<double>::quiet_NaN() : std::stod(pval);
    rows.push_back(row);
  }
  return rows;
}

/* Fields of one tabular placement row. */
struct PlaceRow
{
  std::string id;
  std::string node;
  uint32_t edge;
  double lwr;
  double distance;
};

std::vector<PlaceRow> parse_place(const std::string& text)
{
  std::vector<PlaceRow> rows;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    PlaceRow row;
    std::string edge, lwr, dist;
    if (!std::getline(ls, row.id, '\t')) continue;
    if (!std::getline(ls, row.node, '\t')) continue;
    if (!std::getline(ls, edge, '\t')) continue;
    if (!std::getline(ls, lwr, '\t')) continue;
    if (!std::getline(ls, dist)) continue;
    row.edge = static_cast<uint32_t>(std::stoul(edge));
    row.lwr = std::stod(lwr);
    row.distance = std::stod(dist);
    rows.push_back(row);
  }
  return rows;
}

const DistRow* best_row(const std::vector<DistRow>& rows)
{
  const DistRow* best = nullptr;
  for (const DistRow& row : rows) {
    if (row.is_na) continue;
    if (!best || row.distance < best->distance) best = &row;
  }
  return best;
}

} // namespace

TEST_SUITE_BEGIN("query");

TEST_CASE("dist reports the source reference of an exact substring")
{
  TempDir dir("query-dist");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::string read = substr_of(corpus.base, 5000, 500);
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"exact", read}});

  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string()));
  REQUIRE_FALSE(rows.empty());
  const DistRow* best = best_row(rows);
  REQUIRE(best != nullptr);
  CHECK(best->reference == "refA");
  CHECK(best->distance < 0.01);
}

TEST_CASE("the estimated distance grows with the mutation rate")
{
  TempDir dir("query-mut");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::string read = substr_of(corpus.base, 3000, 2000);
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq,
              {{"m0", read},
               {"m02", mutate(read, 0.02, 1)},
               {"m05", mutate(read, 0.05, 2)},
               {"m10", mutate(read, 0.10, 3)}});

  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string()));
  flat_phmap<std::string, double> best_of;
  for (const DistRow& row : rows) {
    if (row.is_na) continue;
    auto it = best_of.find(row.id);
    if (it == best_of.end() || row.distance < it->second) best_of[row.id] = row.distance;
  }
  REQUIRE(best_of.count("m0") == 1);
  REQUIRE(best_of.count("m02") == 1);
  REQUIRE(best_of.count("m05") == 1);
  REQUIRE(best_of.count("m10") == 1);
  CHECK(best_of["m0"] < best_of["m02"]);
  CHECK(best_of["m02"] < best_of["m05"]);
  CHECK(best_of["m05"] < best_of["m10"]);
  CHECK(best_of["m10"] < 0.4);
}

TEST_CASE("a reverse-complemented query still matches its source")
{
  TempDir dir("query-rc");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::string read = substr_of(corpus.base, 8000, 800);
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"forward", read}, {"reverse", dna_revcomp(read)}});

  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string()));
  flat_phmap<std::string, std::string> best_ref;
  flat_phmap<std::string, double> best_dist;
  for (const DistRow& row : rows) {
    if (row.is_na) continue;
    if (!best_dist.count(row.id) || row.distance < best_dist[row.id]) {
      best_dist[row.id] = row.distance;
      best_ref[row.id] = row.reference;
    }
  }
  REQUIRE(best_ref.count("forward") == 1);
  REQUIRE(best_ref.count("reverse") == 1);
  CHECK(best_ref["forward"] == "refA");
  CHECK(best_ref["reverse"] == "refA");
  CHECK(best_dist["reverse"] == doctest::Approx(best_dist["forward"]).epsilon(0.5));
}

TEST_CASE("dist --no-multi reports exactly one row per query")
{
  TempDir dir("query-nomulti");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q1", substr_of(corpus.base, 100, 500)}, {"q2", substr_of(corpus.base, 9000, 500)}});

  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, false));
  CHECK(rows.size() == 2);
  for (const DistRow& row : rows) CHECK_FALSE(row.is_na);
  CHECK(rows[0].id == "q1");
  CHECK(rows[1].id == "q2");
}

TEST_CASE("dist-max turns distant rows into NA")
{
  TempDir dir("query-distmax");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"unrelated", rand_dna(400, 77)}});
  // Against an unrelated read every hit is far away, so a strict cutoff
  // suppresses the report entirely.
  const std::string strict = dist_queries(index, fq.string(), 4, 2.706, 1e-6, 2, true, true);
  const std::vector<DistRow> rows = parse_dist(strict);
  REQUIRE(rows.size() == 1);
  CHECK(rows[0].is_na);
}

TEST_CASE("the chi-square filter drops rows that are not distinguishable")
{
  TempDir dir("query-filter");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", substr_of(corpus.base, 2000, 1500)}});
  const std::vector<DistRow> unfiltered = parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, true));
  const std::vector<DistRow> filtered = parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, false, true));
  CHECK_FALSE(unfiltered.empty());
  CHECK(filtered.size() <= unfiltered.size());
  CHECK_FALSE(filtered.empty());
  // The best hit always survives the filter.
  const DistRow* best = best_row(unfiltered);
  REQUIRE(best != nullptr);
  bool found = false;
  for (const DistRow& row : filtered) {
    if (row.reference == best->reference) found = true;
  }
  CHECK(found);
}

TEST_CASE("the reported P_VALUE matches the lower tail of the one-sided chi-square test")
{
  // The column is P(chi2_1 < x) = erf(sqrt(x / 2)); the table is allowed 1.6e-4
  // absolute and 4.2e-4 relative deviation, so the test uses looser bounds.
  double worst_abs = 0;
  double worst_rel = 0;
  for (double chisq = 0.0; chisq <= 100.0; chisq += 0.0005) {
    const double exact = std::erf(std::sqrt(chisq / 2.0));
    const double got = chisq_cdf(chisq);
    worst_abs = std::max(worst_abs, std::fabs(got - exact));
    if (exact > 1e-9) {
      worst_rel = std::max(worst_rel, std::fabs(got - exact) / exact);
    }
  }
  CHECK(worst_abs < 2e-4);
  CHECK(worst_rel < 5e-4);

  CHECK(chisq_cdf(0.0) == 0.0);       // tied with the best hit
  CHECK(chisq_cdf(-1.0) == 0.0);      // a statistic cannot be negative
  CHECK(std::isnan(chisq_cdf(std::numeric_limits<double>::quiet_NaN())));
  CHECK(chisq_cdf(2.706) == doctest::Approx(0.9).epsilon(0.01));  // the default --chisq
  CHECK(chisq_cdf(3.841) == doctest::Approx(0.95).epsilon(0.01)); // chi2_1 at 0.05
  CHECK(chisq_cdf(100.0) == 1.0);  // beyond this the two are indistinguishable
  CHECK(chisq_cdf(101.0) == 1.0);  // beyond the table
  double previous = 0.0;
  for (double chisq = 0.0; chisq <= 100.0; chisq += 0.1) {
    const double got = chisq_cdf(chisq);
    CHECK(got >= previous);
    previous = got;
  }
}

TEST_CASE("dist reports a p-value column for every row")
{
  TempDir dir("query-pvalue");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", substr_of(corpus.base, 2000, 1500)}});
  const std::vector<DistRow> rows =
    parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, true, false, true));
  REQUIRE_FALSE(rows.empty());
  for (const DistRow& row : rows) {
    CHECK(row.nfields == 4);
    CHECK_FALSE(row.is_na);
    CHECK(row.p_value >= 0.0);
    CHECK(row.p_value <= 1.0);
  }
  // The closest hit is compared with itself, so it is not distinguishable.
  const DistRow* best = best_row(rows);
  REQUIRE(best != nullptr);
  CHECK(best->p_value == 0.0);
  const DistRow* worst = nullptr;
  for (const DistRow& row : rows) {
    if (!worst || row.distance > worst->distance) worst = &row;
  }
  REQUIRE(worst != nullptr);
  CHECK(worst->p_value > 0.9);
}

TEST_CASE("dist p-values agree with the chi-square filter")
{
  TempDir dir("query-pvalue-filter");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", substr_of(corpus.base, 2000, 1500)}});
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::vector<DistRow> unfiltered =
    parse_dist(dist_queries(index, fq.string(), 4, 2.706, nan, 2, true, true, false, true));
  const std::vector<DistRow> filtered =
    parse_dist(dist_queries(index, fq.string(), 4, 2.706, nan, 2, false, true, false, true));
  REQUIRE_FALSE(unfiltered.empty());
  REQUIRE_FALSE(filtered.empty());
  // A row survives the filter exactly when its statistic is below the cutoff.
  // The default 2.706 is the rounded 90% quantile, so the cutoff P_VALUE is
  // 0.90003.
  const double cutoff = chisq_cdf(2.706);
  for (const DistRow& row : unfiltered) {
    bool present = false;
    for (const DistRow& kept : filtered) {
      if (kept.id == row.id && kept.reference == row.reference) present = true;
    }
    CHECK(present == (row.p_value < cutoff));
  }
}

TEST_CASE("the P_VALUE column is opt-in")
{
  TempDir dir("query-no-pvalue");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", substr_of(corpus.base, 2000, 1500)}});
  const std::vector<DistRow> plain = parse_dist(dist_queries(index, fq.string()));
  REQUIRE_FALSE(plain.empty());
  for (const DistRow& row : plain) {
    CHECK(row.nfields == 3);                          // the old format, by default
    CHECK(std::isnan(row.p_value));                   // and no value read for it
  }
  const std::vector<DistRow> with_column =
    parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, true, false, true));
  REQUIRE(with_column.size() == plain.size());
  for (size_t i = 0; i < plain.size(); ++i) {
    CHECK(with_column[i].reference == plain[i].reference);
    CHECK(with_column[i].distance == plain[i].distance);
    CHECK_FALSE(std::isnan(with_column[i].p_value));
  }
}

TEST_CASE("queries with no usable k-mers are reported as NA")
{
  TempDir dir("query-short");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"short", "ACGTACGTAC"}, {"allN", std::string(200, 'N')}, {"ok", substr_of(corpus.base, 0, 400)}});
  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string()));
  flat_phmap<std::string, bool> seen_na;
  flat_phmap<std::string, bool> seen_ok;
  for (const DistRow& row : rows) {
    if (row.is_na) {
      seen_na[row.id] = true;
    } else {
      seen_ok[row.id] = true;
    }
  }
  CHECK(seen_na.count("short") == 1);
  CHECK(seen_na.count("allN") == 1);
  CHECK(seen_ok.count("ok") == 1);
}

TEST_CASE("dist with hdist-th 0 still matches an exact copy")
{
  TempDir dir("query-hdist0");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"exact", substr_of(corpus.base, 7000, 400)}});
  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string(), 0));
  const DistRow* best = best_row(rows);
  REQUIRE(best != nullptr);
  CHECK(best->reference == "refA");
}

TEST_CASE("compute_branch_lengths splits the edge by the distance ratio")
{
  // d_x = 0.3 and d_y = 0.2 on a b = 0.1 edge: distal = 0.1 * 0.3 / 0.5 = 0.06
  // and pendant = max(0.3 - 0.06, 0.2 - 0.04) = 0.24. A split that admits a
  // pendant must not ask for the leaf fallback.
  bool fallback_used = false;
  const branch_lengths_t split = compute_branch_lengths(0.1, 0.3, 0.2, [&]() {
    fallback_used = true;
    return 0.5;
  });
  CHECK(split.distal == doctest::Approx(0.06));
  CHECK(split.pendant == doctest::Approx(0.24));
  CHECK_FALSE(fallback_used);
}

TEST_CASE("compute_branch_lengths falls back to the closest leaf below x")
{
  // b = 1 is longer than d_x + d_y, so both excesses are negative and the
  // split cannot host a pendant: use the supplied leaf distance.
  bool fallback_used = false;
  const branch_lengths_t split = compute_branch_lengths(1.0, 0.1, 0.1, [&]() {
    fallback_used = true;
    return 0.07;
  });
  CHECK(fallback_used);
  CHECK(split.distal == doctest::Approx(0.5));
  CHECK(split.pendant == doctest::Approx(0.07));

  // No scored leaf below x: the pendant collapses to zero rather than negative.
  const branch_lengths_t empty = compute_branch_lengths(1.0, 0.1, 0.1, []() {
    return std::numeric_limits<double>::quiet_NaN();
  });
  CHECK(empty.distal == doctest::Approx(0.5));
  CHECK(empty.pendant == doctest::Approx(0.0));
}

TEST_CASE("compute_branch_lengths uses 0.33 only for the ratio when the parent is missing")
{
  bool fallback_used = false;
  const branch_lengths_t split = compute_branch_lengths(0.1, 0.03, std::numeric_limits<double>::quiet_NaN(), [&]() {
    fallback_used = true;
    return 0.0;
  });
  // The substitute still positions the split ...
  const double expected_distal = 0.1 * 0.03 / (0.03 + kDefaultParentDistance);
  CHECK(split.distal == doctest::Approx(expected_distal));
  // ... but it must not feed the y-side excess and inflate the pendant.
  const double expected_pendant = 0.03 - expected_distal;
  CHECK(split.pendant == doctest::Approx(expected_pendant));
  CHECK(split.pendant < kDefaultParentDistance - (0.1 - expected_distal));
  CHECK_FALSE(fallback_used);
}

TEST_CASE("compute_branch_lengths takes the larger excess when both distances are real")
{
  // d_x = 0.1 and d_y = 0.5 on a b = 0.1 edge: distal = 0.1 * 0.1 / 0.6.
  // The y-side excess (0.5 - 0.08333 = 0.41667) is larger than the x-side
  // excess (0.1 - 0.01667 = 0.08333), so the pendant comes from d_y.
  bool fallback_used = false;
  const branch_lengths_t split = compute_branch_lengths(0.1, 0.1, 0.5, [&]() {
    fallback_used = true;
    return 0.0;
  });
  CHECK(split.distal == doctest::Approx(0.1 * 0.1 / 0.6));
  CHECK(split.pendant == doctest::Approx(0.5 - (0.1 - 0.1 * 0.1 / 0.6)));
  CHECK_FALSE(fallback_used);
}

TEST_CASE("compute_branch_lengths splits evenly when both distances vanish")
{
  bool fallback_used = false;
  const branch_lengths_t split = compute_branch_lengths(0.4, 0.0, 0.0, [&]() {
    fallback_used = true;
    return 0.0;
  });
  CHECK(split.distal == doctest::Approx(0.2));
  CHECK(split.pendant == doctest::Approx(0.0));
  // Both excesses are negative, so the fallback is what produced the zero.
  CHECK(fallback_used);
}

TEST_CASE("place reports valid edges and normalised weights")
{
  TempDir dir("query-place");
  const Corpus corpus = make_corpus();
  BuildOptions opt = query_opts();
  const std::filesystem::path nwk = dir / "guide.nwk";
  spit(nwk, "((refA:0.05,refB:0.05)AB:0.1,(refC:0.2,refD:0.2)CD:0.1)root;");
  opt.nwk_path = nwk;
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q1", substr_of(corpus.base, 1000, 600)}, {"q2", substr_of(corpus.base, 10000, 600)}});

  const std::vector<PlaceRow> multi = parse_place(place_queries(index, fq.string(), true, 4, 2.706, 2, false, true));
  REQUIRE_FALSE(multi.empty());
  const uint32_t max_edge = index->get_tree()->get_nnodes() - 1;
  flat_phmap<std::string, double> weight_sum;
  flat_phmap<std::string, uint32_t> rows_per_query;
  for (const PlaceRow& row : multi) {
    CHECK(row.edge <= max_edge);
    CHECK(row.lwr >= 0.0);
    CHECK(row.lwr <= 1.0);
    CHECK(row.distance >= 0.0);
    CHECK(row.distance < 1.0);
    weight_sum[row.id] += row.lwr;
    rows_per_query[row.id] += 1;
  }
  CHECK(rows_per_query.count("q1") == 1);
  CHECK(rows_per_query.count("q2") == 1);
  CHECK(weight_sum["q1"] == doctest::Approx(1.0).epsilon(0.02));
  CHECK(weight_sum["q2"] == doctest::Approx(1.0).epsilon(0.02));

  // Single-placement mode reports exactly one edge per query.
  const std::vector<PlaceRow> single = parse_place(place_queries(index, fq.string(), true, 4, 2.706, 2, false, false));
  flat_phmap<std::string, uint32_t> single_rows;
  for (const PlaceRow& row : single) single_rows[row.id] += 1;
  CHECK(single_rows["q1"] == 1);
  CHECK(single_rows["q2"] == 1);
}

TEST_CASE("reported placements keep a non-negative pendant and a distal on the edge")
{
  TempDir dir("query-branch-lengths");
  const Corpus corpus = make_corpus();
  BuildOptions opt = query_opts();
  const std::filesystem::path nwk = dir / "guide.nwk";
  spit(nwk, "((refA:0.05,refB:0.05)AB:0.1,(refC:0.2,refD:0.2)CD:0.1)root;");
  opt.nwk_path = nwk;
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::string read = substr_of(corpus.base, 1000, 600);
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", read}});
  auto qs = std::make_shared<QSeq>(fq.string());
  qs->read_next_batch();
  REQUIRE(qs->get_cbatch_size() == 1);
  IBatch ib(index, qs, 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, true, false, false);
  auto imers_or = std::make_shared<IMers>(index, read.size(), 4);
  auto imers_rc = std::make_shared<IMers>(index, read.size(), 4);
  ib.search_mers(read.data(), read.size(), imers_or, imers_rc);
  ib.summarize_matches(imers_or, imers_rc);

  vec<placement_t> placements;
  REQUIRE(ib.collect_placements(placements));
  REQUIRE_FALSE(placements.empty());
  for (const placement_t& pp : placements) {
    const double b = std::isnan(pp.node->get_blen()) ? 0.0 : pp.node->get_blen();
    INFO("edge " << pp.edge_num << " distal " << pp.distal_length << " pendant " << pp.pendant_length);
    CHECK(pp.pendant_length >= 0.0);
    CHECK(pp.distal_length >= 0.0);
    CHECK(pp.distal_length <= b + 1e-9);
  }
}

TEST_CASE("place emits a jplace-shaped fragment per query")
{
  TempDir dir("query-jplace");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q1", substr_of(corpus.base, 200, 500)}});
  const std::string out = place_queries(index, fq.string(), false, 4, 2.706, 2, false, true);
  INFO("place output: " << out);
  CHECK(out.find("{\"n\" : [\"q1\"], \"p\" : [") != std::string::npos);
  CHECK(out.back() == '}');
  // Every placement array carries the six jplace fields.
  uint32_t placement_lines = 0;
  std::istringstream in(out);
  std::string line;
  while (std::getline(in, line)) {
    const size_t first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] != '[') continue;
    CHECK(std::count(line.begin(), line.end(), ',') == 5); // six fields
    ++placement_lines;
  }
  CHECK(placement_lines >= 1);
}

TEST_CASE("summarize mode accumulates weights without printing placements")
{
  TempDir dir("query-summarize");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q1", substr_of(corpus.base, 300, 500)}, {"q2", substr_of(corpus.base, 4000, 500)}});
  auto qs = std::make_shared<QSeq>(fq.string());
  double total = 0;
  uint32_t entries = 0;
  while (qs->read_next_batch() || !qs->is_batch_finished()) {
    IBatch ib(index, qs, 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, false, true, true, false);
    strstream batch;
    ib.estimate_distances(batch);
    CHECK(batch.str().empty()); // summarize mode prints nothing per batch
    for (const auto& [nd, wcount] : ib.get_summary()) {
      CHECK(nd != nullptr);
      CHECK(wcount > 0.0);
      total += wcount;
      ++entries;
    }
  }
  CHECK(entries > 0);
  CHECK(total > 0.0);
}

TEST_CASE("an empty query file produces no rows")
{
  TempDir dir("query-empty");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");
  const std::filesystem::path fq = dir / "empty.fq";
  spit(fq, "");
  CHECK(parse_dist(dist_queries(index, fq.string())).empty());
}

TEST_CASE("repeated runs report the same results")
{
  TempDir dir("query-determinism");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q1", substr_of(corpus.base, 1500, 700)}, {"q2", mutate(substr_of(corpus.base, 3000, 700), 0.05, 5)}});
  // Distance rows are sorted and carry a per-node estimate, so the report is
  // byte for byte reproducible.
  const std::string first = dist_queries(index, fq.string());
  const std::string second = dist_queries(index, fq.string());
  CHECK(first == second);
  CHECK(sorted_lines(first) == sorted_lines(second));

  // Placement rows are compared as a set: the like-weight-ratio is a ratio of
  // sums over candidates, and which of several equally distant candidates is
  // taken as the reference model still follows the colour map (as it always
  // has), so the last digits of the weights are not fixed.
  const std::string place_first = place_queries(index, fq.string());
  const std::string place_second = place_queries(index, fq.string());
  CHECK(sorted_lines(place_first) == sorted_lines(place_second));
}

TEST_CASE("IBatch::search_mers counts every valid k-mer")
{
  TempDir dir("query-search");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  const std::string read = substr_of(corpus.base, 1000, 500);
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq, {{"q", read}});
  auto qs = std::make_shared<QSeq>(fq.string());
  // The last batch of a file is returned with a false "more to come" flag, so
  // the batch size, not the return value, says whether anything was read.
  qs->read_next_batch();
  REQUIRE(qs->get_cbatch_size() == 1);
  IBatch ib(index, qs, 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, true, false, false);
  auto imers_or = std::make_shared<IMers>(index, read.size(), 4);
  auto imers_rc = std::make_shared<IMers>(index, read.size(), 4);
  CHECK_NOTHROW(ib.search_mers(read.data(), read.size(), imers_or, imers_rc));
  CHECK_NOTHROW(ib.summarize_matches(imers_or, imers_rc));
  // The search found at least one reference to place against.
  vec<placement_t> placements;
  ib.collect_placements(placements);
  CHECK_FALSE(placements.empty());
  for (const placement_t& pp : placements) {
    CHECK(pp.node != nullptr);
    CHECK(pp.edge_num <= index->get_tree()->get_nnodes() - 1);
    CHECK(pp.like_weight_ratio >= 0.0);
    CHECK(pp.like_weight_ratio <= 1.0);
  }
}

TEST_CASE("a query longer than one batch is still reported in order")
{
  TempDir dir("query-long");
  const Corpus corpus = make_corpus();
  const BuildOptions opt = query_opts();
  build_index_from_refs(dir / "index", corpus.refs, opt);
  auto index = load_index_dir(dir / "index");

  // Three records; the batch byte limit forces them into separate batches.
  const std::filesystem::path fq = dir / "reads.fq";
  write_fastq(fq,
              {{"a", substr_of(corpus.base, 0, 900)},
               {"b", substr_of(corpus.base, 2000, 900)},
               {"c", substr_of(corpus.base, 4000, 900)}});
  const std::vector<DistRow> rows = parse_dist(dist_queries(index, fq.string(), 4, 2.706, std::numeric_limits<double>::quiet_NaN(), 2, true, false));
  REQUIRE(rows.size() == 3);
  CHECK(rows[0].id == "a");
  CHECK(rows[1].id == "b");
  CHECK(rows[2].id == "c");
}

TEST_SUITE_END();
