/* Unit tests for the LSH. The hash values are part of the on-disk format, so
 * a few golden values are pinned alongside the structural checks. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

/* Independent, string-level references for the two encodings the query loop
 * computes for every k-mer. */

uint8_t base_code(char c)
{
  switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 4;
  }
}

/* LSH positions are offsets from the low-order end of the encoding: position p
 * of a k-mer is the base that sits at bits 2p..2p+1, i.e. the (k-1-p)-th
 * character of the string. */
uint8_t code_at(const std::string& mer, uint8_t p) { return base_code(mer[mer.size() - 1 - p]); }

/* LSH(x): the 2-bit codes of the selected positions, lowest bit position first,
 * which is what a pext with mask_hash_bp produces. */
uint32_t hash_ref(LSHF& lshf, const std::string& mer)
{
  const vec<uint8_t> ppos = lshf.get_ppos(); // sorted descending
  uint32_t out = 0;
  uint32_t rank = 0;
  for (auto it = ppos.rbegin(); it != ppos.rend(); ++it) {
    out |= static_cast<uint32_t>(code_at(mer, *it)) << (2 * rank);
    ++rank;
  }
  return out;
}

/* drop_ppos(x): one bit per non-selected position in each half of the word, so
 * that folding the halves counts differing positions. */
uint32_t drop_ref(LSHF& lshf, const std::string& mer)
{
  const vec<uint8_t> npos = lshf.get_npos(); // sorted ascending
  uint32_t out = 0;
  for (size_t i = 0; i < npos.size() && i < 16; ++i) {
    const uint8_t code = code_at(mer, npos[i]);
    out |= static_cast<uint32_t>(code & 1u) << i;
    out |= static_cast<uint32_t>(code >> 1) << (16 + i);
  }
  return out;
}

uint64_t encode_lr(const std::string& mer)
{
  uint64_t low = 0, high = 0;
  for (char c : mer) {
    low = (low << 1) | (base_code(c) & 1u);
    high = (high << 1) | (base_code(c) >> 1);
  }
  return (high << 32) | low;
}

uint64_t encode_bp(const std::string& mer)
{
  uint64_t enc = 0;
  for (char c : mer) enc = (enc << 2) | base_code(c);
  return enc;
}

/* Bit-by-bit bit compaction: collect the set bits of `mask` lowest first.
 * This is the definition the pext/compress paths have to reproduce. */
uint64_t compact_ref(uint64_t x, uint64_t mask)
{
  uint64_t out = 0;
  uint32_t bit = 0;
  for (uint32_t i = 0; i < 64; ++i) {
    if (mask & (uint64_t{1} << i)) {
      if (x & (uint64_t{1} << i)) out |= (uint64_t{1} << bit);
      ++bit;
    }
  }
  return out;
}

uint32_t positions_differing(const std::string& a, const std::string& b, const vec<uint8_t>& pos)
{
  uint32_t n = 0;
  for (uint8_t p : pos) n += (code_at(a, p) != code_at(b, p));
  return n;
}

/* k/h pairs that exercise the compact (k - h > 0) and the padded
 * (k - h < 16) shapes of the dropped encoding. */
const std::vector<std::pair<uint8_t, uint8_t>> configs = {{21, 9}, {23, 9}, {27, 11}, {31, 15}, {19, 9}, {31, 16}};

} // namespace

TEST_SUITE_BEGIN("lshf");

TEST_CASE("get_random_positions partitions 0..k-1 into ppos and npos")
{
  gen.seed(1);
  for (const auto& [k, h] : configs) {
    LSHF lshf(k, h, 4, 1, true);
    const vec<uint8_t> ppos = lshf.get_ppos();
    const vec<uint8_t> npos = lshf.get_npos();
    CHECK(lshf.get_k() == k);
    CHECK(lshf.get_h() == h);
    CHECK(lshf.get_m() == 4);
    REQUIRE(ppos.size() == h);
    REQUIRE(npos.size() == k - h);

    // ppos is sorted descending, npos ascending; together they cover 0..k-1.
    CHECK(std::is_sorted(ppos.begin(), ppos.end(), std::greater<uint8_t>()));
    CHECK(std::is_sorted(npos.begin(), npos.end()));
    std::vector<bool> seen(k, false);
    for (uint8_t p : ppos) {
      REQUIRE(p < k);
      CHECK_FALSE(seen[p]);
      seen[p] = true;
    }
    for (uint8_t p : npos) {
      REQUIRE(p < k);
      CHECK_FALSE(seen[p]);
      seen[p] = true;
    }
    for (uint8_t p = 0; p < k; ++p) CHECK(seen[p]);
    CHECK(std::set<uint8_t>(ppos.begin(), ppos.end()).size() == h);
  }
}

TEST_CASE("the LSH positions are seeded deterministically")
{
  gen.seed(7);
  LSHF a(27, 11, 4, 1, true);
  gen.seed(7);
  LSHF b(27, 11, 4, 1, true);
  CHECK(a.get_ppos() == b.get_ppos());
  CHECK(a.get_npos() == b.get_npos());
  gen.seed(8);
  LSHF c(27, 11, 4, 1, true);
  CHECK(a.get_ppos() != c.get_ppos());
}

TEST_CASE("the reconstructing constructor recovers k and h")
{
  gen.seed(3);
  LSHF original(27, 11, 4, 2, false);
  LSHF rebuilt(4, original.get_ppos(), original.get_npos(), 2, false);
  CHECK(rebuilt.get_k() == original.get_k());
  CHECK(rebuilt.get_h() == original.get_h());
  CHECK(rebuilt.get_m() == original.get_m());
  CHECK(rebuilt.get_ppos() == original.get_ppos());
  CHECK(rebuilt.get_npos() == original.get_npos());
  CHECK(rebuilt.check_compatible(std::make_shared<LSHF>(27, 11, 4, 2, false)) == false); // different positions
  CHECK(rebuilt.check_compatible(std::make_shared<LSHF>(4, original.get_ppos(), original.get_npos(), 2, false)));
}

TEST_CASE("compute_hash selects the ppos bits in pext order")
{
  for (const auto& [k, h] : configs) {
    gen.seed(100 + k);
    LSHF lshf(k, h, 4, 1, true);
    for (uint64_t seed = 0; seed < 40; ++seed) {
      const std::string mer = rand_dna(k, seed * 7919 + k);
      const uint32_t from_bits = lshf.compute_hash(encode_bp(mer) & (std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2)));
      CHECK(from_bits == hash_ref(lshf, mer));
      // The hash fits in 2h bits (h == 16 fills the whole word).
      if (2 * h < 32) CHECK(from_bits < (uint32_t{1} << (2 * h)));
    }
  }
}

TEST_CASE("drop_ppos_lr keeps one bit of every non-selected position per half")
{
  for (const auto& [k, h] : configs) {
    gen.seed(200 + k);
    LSHF lshf(k, h, 4, 1, true);
    const uint64_t mask_lr = ((std::numeric_limits<uint64_t>::max() >> (64 - k)) << 32) +
                             ((std::numeric_limits<uint64_t>::max() << 32) >> (64 - k));
    for (uint64_t seed = 0; seed < 40; ++seed) {
      const std::string mer = rand_dna(k, seed * 104729 + k);
      const uint64_t lr = encode_lr(mer) & mask_lr;
      CHECK(lshf.drop_ppos_lr(lr) == drop_ref(lshf, mer));
    }
  }
}

TEST_CASE("the dropped encoding measures the Hamming distance over the non-selected positions")
{
  for (const auto& [k, h] : configs) {
    gen.seed(300 + k);
    LSHF lshf(k, h, 4, 1, true);
    const uint64_t mask_lr = ((std::numeric_limits<uint64_t>::max() >> (64 - k)) << 32) +
                             ((std::numeric_limits<uint64_t>::max() << 32) >> (64 - k));
    const vec<uint8_t> npos = lshf.get_npos();
    for (uint64_t seed = 0; seed < 60; ++seed) {
      const std::string a = rand_dna(k, seed * 31 + 1);
      // b differs from a in at most 5 positions, which is the interesting range.
      std::string b = a;
      std::mt19937_64 rng(seed);
      for (uint32_t m = 0; m < 5; ++m) {
        const uint32_t p = static_cast<uint32_t>(rng() % k);
        const char alt = "ACGT"[rng() & 3];
        b[p] = alt;
      }
      // Positions are measured from the low-order end, so the comparison has
      // to happen on the encoded form, not on the string index.
      const uint32_t expected = positions_differing(a, b, npos);
      CHECK(hdist_lr32(lshf.drop_ppos_lr(encode_lr(a) & mask_lr), lshf.drop_ppos_lr(encode_lr(b) & mask_lr)) == expected);
      CHECK(hdist_lr32(lshf.drop_ppos_lr(encode_lr(a) & mask_lr), lshf.drop_ppos_lr(encode_lr(a) & mask_lr)) == 0);
    }
  }
}

TEST_CASE("drop_ppos_bp keeps the two bits of every non-selected position")
{
  for (const auto& [k, h] : configs) {
    gen.seed(400 + k);
    LSHF lshf(k, h, 4, 1, true);
    const vec<uint8_t> npos = lshf.get_npos();
    const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
    for (uint64_t seed = 0; seed < 40; ++seed) {
      const std::string mer = rand_dna(k, seed * 6151 + 3);
      uint32_t expected = 0;
      for (uint32_t i = 0; i < npos.size(); ++i) {
        expected |= static_cast<uint32_t>(code_at(mer, npos[i])) << (2 * i);
      }
      CHECK(lshf.drop_ppos_bp(encode_bp(mer) & mask_bp) == expected);
    }
  }
}

TEST_CASE("the three hash extractions agree with a bit-by-bit compaction")
{
  // extract_bits() walks the set bits of the mask one at a time; it is the
  // straightforward definition of what the pext/compress paths have to do.
  // This is the guard that keeps the (fast) extraction faithful.
  for (const auto& [k, h] : configs) {
    gen.seed(500 + k);
    LSHF lshf(k, h, 4, 1, true);
    const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
    const uint64_t mask_lr = ((std::numeric_limits<uint64_t>::max() >> (64 - k)) << 32) +
                             ((std::numeric_limits<uint64_t>::max() << 32) >> (64 - k));
    // Rebuild the LSH and dropped masks from the position lists.
    uint64_t hash_mask_bp = 0, drop_mask_lr = 0;
    for (uint8_t p : lshf.get_ppos()) hash_mask_bp += (uint64_t{3} << (p * 2));
    for (uint8_t p : lshf.get_npos()) drop_mask_lr += (0x0000000100000001ull << p);
    for (uint32_t i = 0; i < 16 - (k - h); ++i) drop_mask_lr += uint64_t{1} << (i + k);

    for (uint64_t seed = 0; seed < 200; ++seed) {
      const std::string mer = rand_dna(k, seed * 7919 + k);
      const uint64_t bp = encode_bp(mer) & mask_bp;
      const uint64_t lr = encode_lr(mer) & mask_lr;
      CHECK(lshf.compute_hash(bp) == static_cast<uint32_t>(compact_ref(bp, hash_mask_bp)));
      CHECK(lshf.drop_ppos_lr(lr) == static_cast<uint32_t>(compact_ref(lr, drop_mask_lr)));
    }
  }
}

TEST_CASE("golden hash values pin the on-disk hash configuration")
{
  // ppos/npos taken from a real index (k=27, h=11), so these values describe
  // exactly what an existing database on disk expects.
  const vec<uint8_t> ppos = {26, 25, 22, 16, 14, 12, 11, 5, 4, 3, 1};
  const vec<uint8_t> npos = {0, 2, 6, 7, 8, 9, 10, 13, 15, 17, 18, 19, 20, 21, 23, 24};
  LSHF lshf(4u, ppos, npos, 1u, true);
  const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - 27) * 2);
  const uint64_t mask_lr =
    ((std::numeric_limits<uint64_t>::max() >> (64 - 27)) << 32) + ((std::numeric_limits<uint64_t>::max() << 32) >> (64 - 27));

  struct golden_t
  {
    const char* mer;
    uint32_t hash;
    uint32_t dropped;
  };
  const std::vector<golden_t> golden = {
    {"AAAAAAAAAAAAAAAAAAAAAAAAAAA", 0x00000000u, 0x00000000u},
    {"ACGTACGTACGTACGTACGTACGTACG", 0x00048b6du, 0xd9196ba8u},
    {"TTTTTTTTTTTTTTTTTTTTTTTTTTT", 0x003fffffu, 0xffffffffu},
    {"GATTACAGATTACAGATTACAGATTAC", 0x0020d88cu, 0xca62e26bu},
  };
  for (const golden_t& g : golden) {
    const std::string mer = g.mer;
    REQUIRE(mer.size() == 27);
    CHECK(lshf.compute_hash(encode_bp(mer) & mask_bp) == g.hash);
    CHECK(lshf.drop_ppos_lr(encode_lr(mer) & mask_lr) == g.dropped);
    // The golden values also agree with the independent references.
    CHECK(g.hash == hash_ref(lshf, mer));
    CHECK(g.dropped == drop_ref(lshf, mer));
  }
}

TEST_CASE("check_compatible rejects mismatched configurations")
{
  gen.seed(11);
  LSHF base(27, 11, 4, 1, true);
  auto same = std::make_shared<LSHF>(4u, base.get_ppos(), base.get_npos(), 1u, true);
  {
    // check_compatible() prints a diagnostic on mismatch; keep it out of the
    // test log (it writes to stdout, not stderr).
    CaptureStream quiet_out(std::cout);
    CaptureStream quiet_err(std::cerr);
    CHECK(base.check_compatible(same));
    CHECK(base.check_compatible(nullptr)); // a null lshf is always compatible
    // Different m, k/h, frac or r are all incompatible.
    CHECK_FALSE(base.check_compatible(std::make_shared<LSHF>(27, 11, 8, 1, true)));
    CHECK_FALSE(base.check_compatible(std::make_shared<LSHF>(27, 10, 4, 1, true)));
    CHECK_FALSE(base.check_compatible(std::make_shared<LSHF>(4u, base.get_ppos(), base.get_npos(), 1u, false)));
    CHECK_FALSE(base.check_compatible(std::make_shared<LSHF>(4u, base.get_ppos(), base.get_npos(), 2u, true)));
    // r only matters when the configuration is fractional.
    LSHF whole(27, 11, 4, 1, false);
    std::vector<uint8_t> p = whole.get_ppos(), n = whole.get_npos();
    CHECK(whole.check_compatible(std::make_shared<LSHF>(4u, p, n, 3u, false)));
  }
}

TEST_CASE("data pointers expose the raw position arrays")
{
  const vec<uint8_t> ppos = {26, 25, 22, 16, 14, 12, 11, 5, 4, 3, 1};
  const vec<uint8_t> npos = {0, 2, 6, 7, 8, 9, 10, 13, 15, 17, 18, 19, 20, 21, 23, 24};
  LSHF lshf(4u, ppos, npos, 1u, true);
  REQUIRE(lshf.ppos_data() != nullptr);
  REQUIRE(lshf.npos_data() != nullptr);
  for (size_t i = 0; i < ppos.size(); ++i) CHECK(static_cast<uint8_t>(lshf.ppos_data()[i]) == ppos[i]);
  for (size_t i = 0; i < npos.size(); ++i) CHECK(static_cast<uint8_t>(lshf.npos_data()[i]) == npos[i]);
}

TEST_SUITE_END();
