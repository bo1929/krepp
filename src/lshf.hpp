#ifndef _LSHF_H
#define _LSHF_H

#include "common.hpp"

inline uint64_t compress_mv(uint64_t x, const uint64_t* mv)
{
  uint64_t t;
  t = x & mv[0];
  x = (x ^ t) | (t >> 1);
  t = x & mv[1];
  x = (x ^ t) | (t >> 2);
  t = x & mv[2];
  x = (x ^ t) | (t >> 4);
  t = x & mv[3];
  x = (x ^ t) | (t >> 8);
  t = x & mv[4];
  x = (x ^ t) | (t >> 16);
  t = x & mv[5];
  x = (x ^ t) | (t >> 32);
  return x;
}

class LSHF
{
public:
  LSHF(uint8_t k, uint8_t h, uint32_t m, uint32_t r, bool frac);
  LSHF(uint32_t m, vec<uint8_t> ppos_v, vec<uint8_t> npos_v, uint32_t r, bool frac);
  static bool check_configuration(uint8_t k, uint8_t w, uint8_t h, uint32_t m, uint32_t r, bool frac);
  void set_lshf();
  void get_random_positions();
  uint8_t get_k() { return k; }
  uint8_t get_h() { return h; }
  uint32_t get_m() { return m; }
  vec<uint8_t> get_npos() { return npos_v; }
  vec<uint8_t> get_ppos() { return ppos_v; }
  bool check_compatible(lshf_sptr_t lshf);
  uint32_t compute_hash(uint64_t enc_bp);
  uint32_t drop_ppos_lr(uint64_t enc64_lr);
  uint32_t drop_ppos_bp(uint64_t enc64_bp);
  char* npos_data() { return reinterpret_cast<char*>(npos_v.data()); }
  char* ppos_data() { return reinterpret_cast<char*>(ppos_v.data()); }

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  uint32_t r;
  bool frac;
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  uint64_t mask_drop_lr = 0;
  uint64_t mask_drop_bp = 0;
  uint64_t mask_hash_bp = 0;
  // Compaction masks for compress_mv(), one set per mask above.
  uint64_t mv_hash_bp[6] = {0, 0, 0, 0, 0, 0};
  uint64_t mv_drop_lr[6] = {0, 0, 0, 0, 0, 0};
  uint64_t mv_drop_bp[6] = {0, 0, 0, 0, 0, 0};
};

#endif
