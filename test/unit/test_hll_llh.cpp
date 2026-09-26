/* Unit tests for src/hyperloglog.hpp, src/hdhistllh.hpp and Minfo. */

#include "test_helpers.hpp"
#include <boost/math/tools/minima.hpp>

using namespace ktest;

TEST_SUITE_BEGIN("hll");

TEST_CASE("HyperLogLog validates its register width")
{
  CHECK_NOTHROW(hll::HyperLogLog(4));
  CHECK_NOTHROW(hll::HyperLogLog(30));
  CHECK_THROWS_AS(hll::HyperLogLog(3), std::invalid_argument);
  CHECK_THROWS_AS(hll::HyperLogLog(31), std::invalid_argument);
  CHECK(hll::HyperLogLog(12).registerSize() == 4096);
}

TEST_CASE("HyperLogLog estimates a known cardinality")
{
  const uint32_t b = 14;
  for (uint64_t n : {0ull, 10ull, 1000ull, 50000ull, 500000ull}) {
    hll::HyperLogLog hll(b);
    for (uint64_t i = 0; i < n; ++i) hll.add(xur64_hash(i));
    const double est = hll.estimate();
    if (n == 0) {
      CHECK(est == doctest::Approx(0.0).epsilon(0.05));
    } else {
      // The standard error is 1.04/sqrt(2^b) ~ 0.8% for b = 14.
      CHECK(est == doctest::Approx(static_cast<double>(n)).epsilon(0.10));
    }
  }
}

TEST_CASE("HyperLogLog counts duplicates only once")
{
  hll::HyperLogLog hll(14);
  for (uint32_t rep = 0; rep < 20; ++rep) {
    for (uint64_t i = 0; i < 10000; ++i) hll.add(xur64_hash(i));
  }
  CHECK(hll.estimate() == doctest::Approx(10000.0).epsilon(0.10));
}

TEST_CASE("HyperLogLog merges, clears and swaps registers")
{
  hll::HyperLogLog a(12);
  hll::HyperLogLog b(12);
  for (uint64_t i = 0; i < 20000; ++i) {
    if (i % 2 == 0) {
      a.add(xur64_hash(i));
    } else {
      b.add(xur64_hash(i));
    }
  }
  a.merge(b);
  CHECK(a.estimate() == doctest::Approx(20000.0).epsilon(0.12));

  hll::HyperLogLog c(8);
  CHECK_THROWS_AS(a.merge(c), std::invalid_argument);

  hll::HyperLogLog copy(12);
  copy.swap(a);
  CHECK(copy.estimate() == doctest::Approx(20000.0).epsilon(0.12));
  copy.clear();
  CHECK(copy.estimate() == doctest::Approx(0.0).epsilon(0.05));
}

TEST_CASE("the HIP variant stays within the same tolerance")
{
  hll::HyperLogLogHIP hll(14);
  for (uint64_t i = 0; i < 50000; ++i) hll.add(xur64_hash(i));
  CHECK(hll.estimate() == doctest::Approx(50000.0).epsilon(0.15));
  hll.clear();
  CHECK(hll.estimate() == doctest::Approx(0.0));
}

TEST_CASE("the leading-zero helper matches its contract")
{
  // The rank is clamped to the register width plus one.
  CHECK(_clzll(0, 10) == 11);
  CHECK(_clzll(uint64_t{1} << 63, 10) == 1);
  CHECK(_clzll(uint64_t{1} << 62, 10) == 2);
  CHECK(_clzll(uint64_t{1} << 53, 10) == 11);
  CHECK(_clzll(1, 10) == 11);
}

TEST_SUITE_END();

TEST_SUITE_BEGIN("hdhistllh");

namespace {

/* A histogram with `nmatch` matched k-mers spread over the hdist bins and `uc`
 * unmatched k-mers, which is what Minfo hands to the likelihood. */
struct Hist
{
  std::vector<double> v;
  double uc;
};

Hist make_hist(uint32_t hdist_th, uint32_t nmatch, double uc)
{
  Hist h;
  h.v.assign(hdist_th + 1, 0.0);
  for (uint32_t i = 0; i < nmatch; ++i) h.v[i % (hdist_th + 1)] += 1.0;
  h.uc = uc;
  return h;
}

} // namespace

TEST_CASE("the negative log-likelihood is minimised near the true distance")
{
  const uint32_t hdist_th = 4;
  const uint32_t k = 27;
  const uint32_t h = 11;
  optimize::HDistHistLLH llh(h, k, hdist_th);

  // A query with 30% of its k-mers missing and the rest spread over the bins.
  Hist hist = make_hist(hdist_th, 70, 30);
  llh.set_parameters(hist.v.data(), hist.uc, 1.0);

  const double d_at_min = boost::math::tools::brent_find_minima(llh, 1e-10, 0.5, 24).first;
  const double f_min = llh(d_at_min);
  CHECK(d_at_min > 0.0);
  CHECK(d_at_min < 0.5);
  for (double d : {1e-6, 0.01, 0.05, 0.2, 0.35, 0.49}) {
    CHECK(llh(d) >= f_min - 1e-9);
  }

  // Matches concentrated in the zero-mismatch bin put the estimate near zero.
  Hist near_perfect = make_hist(hdist_th, 100, 0);
  near_perfect.v[0] = 100;
  for (uint32_t x = 1; x <= hdist_th; ++x) near_perfect.v[x] = 0;
  llh.set_parameters(near_perfect.v.data(), near_perfect.uc, 1.0);
  const double d_perfect = boost::math::tools::brent_find_minima(llh, 1e-10, 0.5, 24).first;
  CHECK(d_perfect < 0.05);

  // A query where nothing was found at all is pushed towards a large distance.
  Hist nothing = make_hist(hdist_th, 0, 100);
  llh.set_parameters(nothing.v.data(), nothing.uc, 1.0);
  const double d_nothing = boost::math::tools::brent_find_minima(llh, 1e-10, 0.5, 24).first;
  CHECK(d_nothing > d_perfect);
}

TEST_SUITE_END();

TEST_SUITE_BEGIN("minfo");

TEST_CASE("Minfo counts each query position once and keeps the best distance")
{
  const uint32_t hdist_th = 4;
  Minfo mi(hdist_th, 100, 1.0);
  // No matches yet: the tau count is zero.
  CHECK(mi.get_leq_tau(2) == doctest::Approx(0.0));

  // First match at position 0 with hamming distance 3.
  mi.update_match(0x1234, 0, 3);
  CHECK(mi.get_leq_tau(2) == doctest::Approx(0.0));
  CHECK(mi.get_leq_tau(3) == doctest::Approx(1.0));

  // A better match at the same position replaces the previous one.
  mi.update_match(0x1235, 0, 1);
  CHECK(mi.get_leq_tau(0) == doctest::Approx(0.0));
  CHECK(mi.get_leq_tau(1) == doctest::Approx(1.0));
  CHECK(mi.get_leq_tau(3) == doctest::Approx(1.0));

  // A worse match at the same position is ignored.
  mi.update_match(0x1236, 0, 4);
  CHECK(mi.get_leq_tau(1) == doctest::Approx(1.0));
  CHECK(mi.get_leq_tau(3) == doctest::Approx(1.0));

  // A new position counts separately.
  mi.update_match(0x1237, 1, 2);
  CHECK(mi.get_leq_tau(1) == doctest::Approx(1.0));
  CHECK(mi.get_leq_tau(2) == doctest::Approx(2.0));
}

TEST_CASE("Minfo optimises the likelihood and reports a finite distance")
{
  const uint32_t hdist_th = 4;
  const uint32_t k = 27;
  const uint32_t h = 11;
  optimize::HDistHistLLH llh(h, k, hdist_th);
  Minfo mi(hdist_th, 100, 1.0);
  for (uint32_t pos = 0; pos < 60; ++pos) {
    mi.update_match(0x1000 + pos, pos, pos % 3);
  }
  mi.optimize_likelihood(llh);
  const double d = mi.jukes_cantor_dist();
  CHECK(std::isfinite(d));
  CHECK(d > 0.0);
  CHECK(d < 1.0);
  // jukes_cantor_dist() transforms the estimate, so it is not the argmin
  // itself; the ratio is still non-negative everywhere and positive away from
  // the optimum.
  CHECK(mi.likelihood_ratio(d, llh) >= 0.0);
  CHECK(mi.likelihood_ratio(0.45, llh) > 0.0);
}

TEST_CASE("Minfo accumulates children with add")
{
  const uint32_t hdist_th = 4;
  Minfo parent(hdist_th);
  auto child = std::make_shared<Minfo>(hdist_th, 100, 1.0);
  child->update_match(1, 0, 1);
  child->update_match(2, 1, 2);
  parent.add(child, 1.0);
  CHECK(parent.get_leq_tau(2) == doctest::Approx(2.0));

  auto sibling = std::make_shared<Minfo>(hdist_th, 100, 1.0);
  sibling->update_match(3, 5, 0);
  parent.add(sibling, 0.5);
  CHECK(parent.get_leq_tau(0) == doctest::Approx(0.5));
  CHECK(parent.get_leq_tau(2) == doctest::Approx(2.5));

  // join() halves the counts when the destination already tracks nmers.
}

TEST_SUITE_END();
