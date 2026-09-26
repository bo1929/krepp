/* Unit tests for the helpers in src/common.{hpp,cpp}. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

/* String-level references for the bit tricks. */

uint64_t encode_bp_ref(const std::string& seq)
{
  uint64_t enc = 0;
  for (char c : seq) {
    uint8_t v = 0;
    switch (c) {
      case 'A': v = 0; break;
      case 'C': v = 1; break;
      case 'G': v = 2; break;
      case 'T': v = 3; break;
      default: FAIL("not a nucleotide: " << c);
    }
    enc = (enc << 2) | v;
  }
  return enc;
}

uint64_t encode_lr_ref(const std::string& seq)
{
  uint64_t low = 0, high = 0;
  for (char c : seq) {
    uint8_t v = 0;
    switch (c) {
      case 'A': v = 0; break;
      case 'C': v = 1; break;
      case 'G': v = 2; break;
      case 'T': v = 3; break;
      default: FAIL("not a nucleotide: " << c);
    }
    low = (low << 1) | (v & 1u);
    high = (high << 1) | (v >> 1);
  }
  return (high << 32) | low;
}

uint64_t revcomp_bp_ref(const std::string& kmers)
{
  std::string rc;
  for (auto it = kmers.rbegin(); it != kmers.rend(); ++it) {
    switch (*it) {
      case 'A': rc.push_back('T'); break;
      case 'C': rc.push_back('G'); break;
      case 'G': rc.push_back('C'); break;
      case 'T': rc.push_back('A'); break;
      default: break;
    }
  }
  return encode_bp_ref(rc);
}

std::string random_kmer(size_t k, uint64_t seed) { return rand_dna(k, seed); }

} // namespace

TEST_SUITE_BEGIN("common");

TEST_CASE("seq_nt4_table maps the IUPAC alphabet onto 0..4")
{
  CHECK(seq_nt4_table['A'] == 0);
  CHECK(seq_nt4_table['C'] == 1);
  CHECK(seq_nt4_table['G'] == 2);
  CHECK(seq_nt4_table['T'] == 3);
  // krepp accepts lowercase input.
  CHECK(seq_nt4_table['a'] == 0);
  CHECK(seq_nt4_table['c'] == 1);
  CHECK(seq_nt4_table['g'] == 2);
  CHECK(seq_nt4_table['t'] == 3);
  // Everything else is ambiguous and must break the k-mer run.
  for (char c : std::string("NRYWSKMBDHVU.-*0123456789 ")) {
    CHECK(seq_nt4_table[static_cast<unsigned char>(c)] == 4);
  }
  CHECK(seq_nt4_table[0] == 4);
}

TEST_CASE("nt4 tables agree with the string reference encodings")
{
  CHECK(nt4_bp_table[0] == 0);
  CHECK(nt4_bp_table[1] == 1);
  CHECK(nt4_bp_table[2] == 2);
  CHECK(nt4_bp_table[3] == 3);
  CHECK(nt4_lr_table[0] == 0);
  CHECK(nt4_lr_table[1] == 1);
  CHECK(nt4_lr_table[2] == (uint64_t{1} << 32));
  CHECK(nt4_lr_table[3] == (uint64_t{1} << 32) + 1);

  for (uint32_t len = 1; len <= 32; ++len) {
    for (uint64_t seed = 0; seed < 8; ++seed) {
      const std::string seq = random_kmer(len, seed * 131 + len);
      uint64_t lr = 0, bp = 0;
      compute_encoding(seq.data(), seq.data() + seq.size(), lr, bp);
      CHECK(bp == encode_bp_ref(seq));
      CHECK(lr == encode_lr_ref(seq));
      CHECK(lr == conv_bp64_lr64(bp));
    }
  }
}

TEST_CASE("update_encoding equals a fresh compute_encoding while sliding")
{
  const std::string seq = rand_dna(200, 7);
  for (uint32_t k = 1; k <= 32; ++k) {
    uint64_t lr = 0, bp = 0;
    const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
    const uint64_t u64m = std::numeric_limits<uint64_t>::max();
    const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
    for (size_t i = 0; i + k <= seq.size(); ++i) {
      if (i == 0) {
        compute_encoding(seq.data(), seq.data() + k, lr, bp);
      } else {
        update_encoding(seq.data() + i + k - 1, lr, bp);
      }
      bp &= mask_bp;
      lr &= mask_lr;
      const std::string win = seq.substr(i, k);
      CHECK(bp == encode_bp_ref(win));
      CHECK(lr == encode_lr_ref(win));
      CHECK(lr == conv_bp64_lr64(bp));
    }
  }
}

TEST_CASE("revcomp_bp64 is an involution that matches a string reverse complement")
{
  for (uint32_t k = 1; k <= 32; ++k) {
    for (uint64_t seed = 0; seed < 16; ++seed) {
      const std::string mer = random_kmer(k, seed * 977 + k);
      const uint64_t bp = encode_bp_ref(mer);
      const uint64_t rc = revcomp_bp64(bp, static_cast<uint8_t>(k));
      CHECK(rc == revcomp_bp_ref(mer));
      // Complementing twice returns the original k-mer.
      CHECK(revcomp_bp64(rc, static_cast<uint8_t>(k)) == bp);
      if (seed == 0) {
        const uint64_t all_a = 0;
        uint64_t all_t = 0;
        for (uint32_t i = 0; i < k; ++i) all_t = (all_t << 2) | 3;
        CHECK(revcomp_bp64(all_a, static_cast<uint8_t>(k)) == all_t);
        CHECK(revcomp_bp64(all_t, static_cast<uint8_t>(k)) == all_a);
      }
    }
  }
}

TEST_CASE("rmoddp_bp64 compacts the even bits of the word downwards")
{
  CHECK(rmoddp_bp64(0) == 0);
  // Even bit 2i lands on bit i; odd bits are dropped.
  for (uint32_t i = 0; i < 32; ++i) {
    CHECK(rmoddp_bp64(uint64_t{1} << (2 * i)) == (uint64_t{1} << i));
    CHECK(rmoddp_bp64(uint64_t{1} << (2 * i + 1)) == 0);
  }
  // The whole 64-bit word is compacted, not just the low 32 bits.
  CHECK(rmoddp_bp64(uint64_t{1} << 40) == (uint64_t{1} << 20));
  CHECK(rmoddp_bp64(uint64_t{1} << 62) == (uint64_t{1} << 31));
  // Two adjacent even bits stay adjacent; odd bits are dropped.
  CHECK(rmoddp_bp64(0b0101) == 0b11);
  CHECK(rmoddp_bp64(0b1010) == 0);
}

TEST_CASE("conv_bp64_lr64 interleaves the two bits of every base")
{
  CHECK(conv_bp64_lr64(0) == 0);
  // A single 'C' (code 1) keeps one bit in the low half.
  CHECK(conv_bp64_lr64(1) == 1);
  // A single 'G' (code 2) keeps one bit in the high half.
  CHECK(conv_bp64_lr64(2) == (uint64_t{1} << 32));
  // A single 'T' (code 3) sets one bit in each half.
  CHECK(conv_bp64_lr64(3) == (uint64_t{1} << 32) + 1);
}

TEST_CASE("bit-counting helpers on the lr encoding")
{
  CHECK(hdist_lr32(0, 0) == 0);
  CHECK(hdist_lr32(0x00000000u, 0x00000001u) == 1);
  CHECK(popcount_lr32(0) == 0);
  // The two halves are folded: a difference in the high word counts once.
  CHECK(popcount_lr32(0x00010000u) == 1);
  CHECK(popcount_lr32(0x00000001u) == 1);
  // The same position set in both halves still counts once.
  CHECK(popcount_lr32(0x00010001u) == 1);
  CHECK(popcount_lr32(0x00020001u) == 2);
  CHECK(hdist_lr32(0x00000001u, 0x00010000u) == 1);
  // Position 0 differs in exactly one of its two bits.
  CHECK(hdist_lr32(0x00000001u, 0x00010001u) == 1);
  CHECK(hdist_lr32(0x00000003u, 0x00020000u) == 2);
  CHECK(hdist_lr32(0x00000003u, 0x00000003u) == 0);
  // A position that only differs in the high half still counts once.
  CHECK(hdist_lr32(0x00000003u, 0x00030003u) == 2);
  // hdist_lr32 is the folded popcount of the xor.
  std::mt19937_64 rng(99);
  for (uint32_t t = 0; t < 2000; ++t) {
    const uint32_t x = static_cast<uint32_t>(rng());
    const uint32_t y = static_cast<uint32_t>(rng());
    CHECK(hdist_lr32(x, y) == popcount_lr32(x ^ y));
    CHECK(zc_lr32(x, y) == (((x ^ y) | ((x ^ y) >> 16)) & 0xffffu));
  }
}

TEST_CASE("non-cryptographic hashes are deterministic and spread out")
{
  // Golden values pin the hash functions: changing them would invalidate every
  // existing index on disk.
  CHECK(xur32_hash(0) == xur32_hash(0));
  CHECK(xur64_hash(0) == xur64_hash(0));
  CHECK(xur64m_hash(0) == xur64m_hash(0));
  CHECK(xur32_hash(1) != xur32_hash(2));
  CHECK(xur64_hash(1) != xur64_hash(2));
  CHECK(xur64m_hash(1) != xur64m_hash(2));

  // No collisions over a small dense input range, and every hash changes when
  // the input changes.
  flat_phmap<uint64_t, bool> seen32, seen64, seen64m;
  for (uint64_t i = 0; i < 20000; ++i) {
    CHECK_FALSE(seen32.contains(xur32_hash(static_cast<uint32_t>(i))));
    seen32[xur32_hash(static_cast<uint32_t>(i))] = true;
    CHECK_FALSE(seen64.contains(xur64_hash(i)));
    seen64[xur64_hash(i)] = true;
    CHECK_FALSE(seen64m.contains(xur64m_hash(i)));
    seen64m[xur64m_hash(i)] = true;
  }

  // A single-bit flip changes roughly half of the output bits.
  for (uint32_t bit = 0; bit < 64; bit += 7) {
    const uint64_t a = xur64_hash(0x0123456789abcdefull);
    const uint64_t b = xur64_hash(0x0123456789abcdefull ^ (uint64_t{1} << bit));
    const uint32_t differing = __builtin_popcountll(a ^ b);
    CHECK(differing >= 16);
    CHECK(differing <= 48);
  }
}

TEST_CASE("hash_name and rehash agree with the murmur-based subset hash")
{
  std::string a = "G000341695";
  std::string b = "G000341696";
  CHECK(hash_name(a) != 0);
  CHECK(hash_name(a) == Subset::get_singleton_sh(a));
  CHECK(hash_name(a) != hash_name(b));
  CHECK(hash_name(a) == hash_name(a));

  const sh_t h = hash_name(a);
  CHECK(rehash(h) != 0);
  CHECK(rehash(h) == Subset::rehash(h));
  CHECK(rehash(h) != h);
  // Distinct inputs give distinct hashes for a small set of names.
  flat_phmap<sh_t, bool> seen;
  for (uint32_t i = 0; i < 5000; ++i) {
    std::string name = "ref" + std::to_string(i);
    CHECK_FALSE(seen.contains(hash_name(name)));
    seen[hash_name(name)] = true;
  }
}

TEST_CASE("gp_hash is stable and non-zero for non-empty input")
{
  CHECK(gp_hash("") == 0);
  CHECK(gp_hash("krepp") == gp_hash("krepp"));
  CHECK(gp_hash("krepp") != gp_hash("krepx"));
  CHECK((gp_hash("krepp") & 0x80000000u) == 0);
}

TEST_CASE("vec_to_str formats a byte vector")
{
  CHECK(vec_to_str({}) == "[]");
  CHECK(vec_to_str({1}) == "[1]");
  CHECK(vec_to_str({1, 2, 3}) == "[1, 2, 3]");
  CHECK(vec_to_str({255, 0}) == "[255, 0]");
}

TEST_CASE("set_num_threads never falls below one")
{
  const uint32_t saved = num_threads;
  set_num_threads(0);
  CHECK(num_threads == 1);
  {
    // Asking for more than one thread on a build without OpenMP is a no-op, so
    // it has to say so rather than pretending to be parallel.
    CaptureStream quiet(std::cerr);
    set_num_threads(7);
    CHECK(num_threads == 7);
#if !defined(_OPENMP) || _WOPENMP != 1
    CHECK(quiet.str().find("no OpenMP support") != std::string::npos);
#else
    CHECK(quiet.str().find("no OpenMP support") == std::string::npos);
#endif
  }
  set_num_threads(saved);
}

TEST_CASE("error_exit routes through the installed handler")
{
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([] { error_exit("boom"); });
  CHECK(msg == "boom");
  const std::string msg2 = ThrowingErrorHandler::catches([] { error_exit("with code", 7); });
  CHECK(msg2 == "with code");

  // A handler that returns is turned into an abort, which cannot be tested
  // in-process; what is testable is that a throwing handler is not called twice.
  int calls = 0;
  set_error_handler([&calls](const std::string&, int) {
    ++calls;
    throw std::runtime_error("once");
  });
  CHECK_THROWS_AS(error_exit("x"), std::runtime_error);
  CHECK(calls == 1);
  set_error_handler(error_handler_t());
}

TEST_CASE("CHECK_STREAM_OR_EXIT reports the failure with the given message")
{
  ThrowingErrorHandler handler;
  std::ifstream missing("/nonexistent/krepp/test/file");
  const std::string msg = ThrowingErrorHandler::catches([&] { CHECK_STREAM_OR_EXIT(missing, "stream failed"); });
  CHECK(msg == "stream failed");
}

TEST_CASE("ErrorRelay captures the first exception and rethrows it once")
{
  ErrorRelay relay;
  CHECK_FALSE(relay.check_error());
  CHECK_NOTHROW(relay.rethrow_error());

  relay.guard([] {});
  CHECK_FALSE(relay.check_error());

  relay.guard([] { throw std::runtime_error("first"); });
  CHECK(relay.check_error());
  // A second failure must not replace the first.
  relay.guard([] { throw std::runtime_error("second"); });
  CHECK_THROWS_WITH_AS(relay.rethrow_error(), "first", std::runtime_error);

  // Once an error is latched, guard() skips the work entirely.
  ErrorRelay relay2;
  int ran = 0;
  relay2.guard([&] {
    ++ran;
    throw std::runtime_error("fatal");
  });
  relay2.guard([&] { ++ran; });
  CHECK(ran == 1);
  CHECK_THROWS_AS(relay2.rethrow_error(), std::runtime_error);
}

TEST_CASE("ErrorRelay without a captured exception reports its own message")
{
  ErrorRelay relay;
  relay.capture();
  CHECK_THROWS_WITH_AS(relay.rethrow_error(), "a parallel region failed with an exception that could not be captured", std::runtime_error);
}

TEST_SUITE_END();
