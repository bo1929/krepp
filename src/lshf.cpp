#include "lshf.hpp"

LSHF::LSHF(uint8_t k, uint8_t h, uint32_t m, uint32_t r, bool frac)
  : k(k)
  , h(h)
  , m(m)
  , r(r)
  , frac(frac)
{
  get_random_positions();
  set_lshf();
}

bool LSHF::check_configuration(uint8_t k, uint8_t w, uint8_t h, uint32_t m, uint32_t r, bool frac)
{
  (void)frac;
  bool is_invalid = true;
  if ((is_invalid = (w < k))) {
    error_exit("The minimizer window (-w) must be at least the k-mer length (-k).");
  }
  if ((is_invalid = (h < 9))) {
    error_exit("The number of LSH positions (-h) must be at least 9.");
  }
  if ((is_invalid = (h > 15))) {
    error_exit("The number of LSH positions (-h) must be at most 15.");
  }
  if ((is_invalid = (k > 31))) {
    error_exit("The k-mer length (-k) must be at most 31.");
  }
  if ((is_invalid = (k < 19))) {
    error_exit("The k-mer length (-k) must be at least 19.");
  }
  if ((is_invalid = ((k - h) > 16))) {
    error_exit("For compact k-mer encodings, h must be >= k-16.");
  }
  if ((is_invalid = (m == 0))) {
    error_exit("The modulo value (-m) must be positive.");
  }
  if ((is_invalid = (r >= m))) {
    error_exit("The LSH residue (-r) must be smaller than the modulo value (-m).");
  }
  return !is_invalid;
}

static void gen_compress_mv(uint64_t m, uint64_t* mv)
{
  uint64_t mk = ~m << 1;
  for (uint32_t i = 0; i < 6; ++i) {
    uint64_t mp = mk ^ (mk << 1);
    mp ^= mp << 2;
    mp ^= mp << 4;
    mp ^= mp << 8;
    mp ^= mp << 16;
    mp ^= mp << 32;
    mv[i] = mp & m;
    m = (m ^ mv[i]) | (mv[i] >> (1u << i));
    mk &= ~mp;
  }
}

void LSHF::set_lshf()
{
  for (int i = npos_v.size() - 1; i >= 0; --i) {
    mask_drop_lr += (0x0000000100000001ull << npos_v[i]);
    mask_drop_bp += (0x0000000000000003ull << (npos_v[i] * 2));
  }
  // Pad the dropped lr encoding up to 16 bits per half. The bound is signed:
  // for (invalid) configurations with k-h > 16 it is negative and the loop is
  // skipped; the old unsigned comparison wrapped it to ~2^32 iterations and
  // shifted past bit 63 (UB). check_configuration rejects such configs at the
  // CLI level; this keeps direct constructions safe too.
  const int pad = 16 - (static_cast<int>(k) - static_cast<int>(h));
  for (int i = 0; i < pad; ++i) {
    mask_drop_lr += 0x0000000000000001ull << (i + static_cast<int>(k));
  }
  for (int i = ppos_v.size() - 1; i >= 0; --i) {
    mask_hash_bp += (0x0000000000000003ull << (ppos_v[i] * 2));
  }
  gen_compress_mv(mask_hash_bp, mv_hash_bp);
  gen_compress_mv(mask_drop_lr, mv_drop_lr);
  gen_compress_mv(mask_drop_bp, mv_drop_bp);
}

#if defined(__BMI2__)
uint32_t LSHF::compute_hash(uint64_t enc64_bp) { return static_cast<uint32_t>(_pext_u64(enc64_bp, mask_hash_bp)); }

uint32_t LSHF::drop_ppos_lr(uint64_t enc64_lr) { return static_cast<uint32_t>(_pext_u64(enc64_lr, mask_drop_lr)); }

uint32_t LSHF::drop_ppos_bp(uint64_t enc64_bp) { return static_cast<uint32_t>(_pext_u64(enc64_bp, mask_drop_bp)); }
#else
uint32_t LSHF::compute_hash(uint64_t enc64_bp)
{
  return static_cast<uint32_t>(compress_mv(enc64_bp & mask_hash_bp, mv_hash_bp));
}

uint32_t LSHF::drop_ppos_lr(uint64_t enc64_lr)
{
  return static_cast<uint32_t>(compress_mv(enc64_lr & mask_drop_lr, mv_drop_lr));
}

uint32_t LSHF::drop_ppos_bp(uint64_t enc64_bp)
{
  return static_cast<uint32_t>(compress_mv(enc64_bp & mask_drop_bp, mv_drop_bp));
}
#endif

void LSHF::get_random_positions()
{
  uint8_t n;
  assert(h <= 16);
  assert(h < k);
  std::uniform_int_distribution<uint8_t> distrib(0, k - 1);
  while (ppos_v.size() < h) {
    n = distrib(gen);
    if (!std::count(ppos_v.begin(), ppos_v.end(), n)) {
      ppos_v.push_back(n);
    }
  }
  std::sort(ppos_v.begin(), ppos_v.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < k; ++i) {
    if (ix_pos < h && i == ppos_v[ix_pos])
      ix_pos++;
    else
      npos_v.push_back(i);
  }
  std::sort(ppos_v.begin(), ppos_v.end(), std::greater<uint8_t>());
}

LSHF::LSHF(uint32_t m, vec<uint8_t> ppos_v, vec<uint8_t> npos_v, uint32_t r, bool frac)
  : m(m)
  , r(r)
  , frac(frac)
  , npos_v(npos_v)
  , ppos_v(ppos_v)
{
  const size_t kk = npos_v.size() + ppos_v.size();
  const size_t hh = ppos_v.size();
  // Metadata/sketch files are untrusted input: mirror check_configuration's
  // bounds so a corrupt position vector is a clean load error instead of
  // shift UB in set_lshf().
  if (kk < 19 || kk > 31 || hh < 9 || hh > 15 || (kk - hh) > 16) {
    error_exit("Invalid LSH configuration in the index/sketch metadata (k=" + std::to_string(kk) +
               ", h=" + std::to_string(hh) + "); the file may be corrupt.");
  }
  for (uint8_t p : ppos_v) {
    if (p >= kk) {
      error_exit("Invalid LSH position in the index/sketch metadata; the file may be corrupt.");
    }
  }
  for (uint8_t p : npos_v) {
    if (p >= kk) {
      error_exit("Invalid LSH position in the index/sketch metadata; the file may be corrupt.");
    }
  }
  k = static_cast<uint8_t>(kk);
  h = static_cast<uint8_t>(hh);
  set_lshf();
}

bool LSHF::check_compatible(lshf_sptr_t lshf)
{
  if (!lshf) return true;
  bool is_compatible = (lshf->m == m) && (lshf->h == h) && (lshf->k == k) && (lshf->frac == frac) &&
                       (lshf->npos_v == npos_v) && (lshf->ppos_v == ppos_v);
  if (is_compatible && frac && (lshf->r != r)) {
    is_compatible = false;
  }
  if (!is_compatible) {
    std::cout << "m: " << static_cast<uint32_t>(m) << "/" << static_cast<uint32_t>(lshf->m) << std::endl;
    std::cout << "h: " << static_cast<uint32_t>(h) << "/" << static_cast<uint32_t>(lshf->h) << std::endl;
    std::cout << "k: " << static_cast<uint32_t>(k) << "/" << static_cast<uint32_t>(lshf->k) << std::endl;
    std::cout << "frac: " << frac << "/" << lshf->frac << std::endl;
    std::cout << "r: " << r << "/" << lshf->r << std::endl;
    // Bound the debug dumps by BOTH configurations' sizes: on a mismatch the
    // other object's vectors can be shorter, and indexing them with this
    // object's h / k-h read out of bounds.
    const uint32_t h_both = std::min<uint32_t>(h, lshf->h);
    std::cout << "ppos_v:";
    for (uint32_t i = 0; i < h_both; ++i) {
      std::cout << " " << static_cast<uint32_t>(ppos_v[i]) << "/" << static_cast<uint32_t>(lshf->ppos_v[i]);
    }
    std::cout << std::endl;
    const uint32_t n_both = std::min<uint32_t>(k - h, lshf->k - lshf->h);
    std::cout << "npos_v:";
    for (uint32_t i = 0; i < n_both; ++i) {
      std::cout << " " << static_cast<uint32_t>(npos_v[i]) << "/" << static_cast<uint32_t>(lshf->npos_v[i]);
    }
    std::cout << std::endl;
    return false;
  } else {
    return true;
  }
}
