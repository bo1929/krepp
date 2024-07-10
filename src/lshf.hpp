#ifndef _LSHF_H
#define _LSHF_H

#include "common.hpp"
#include "keremet.hpp"
#include "rqseq.hpp"

class LSHF
{
public:
  LSHF(uint8_t k, uint8_t h, uint32_t m);
  LSHF(uint32_t m, vec<uint8_t> ppos_v, vec<uint8_t> npos_v);
  void set_lshashf();
  void get_random_positions();
  uint32_t compute_hash(uint64_t enc_bp);
  uint32_t drop_ppos_lr(uint64_t enc64_lr);
  uint32_t drop_ppos_bp(uint64_t enc64_bp);
  bool check_compatible(lshf_sptr_t lshashf);
  char* npos_data() { return reinterpret_cast<char*>(npos_v.data()); }
  char* ppos_data() { return reinterpret_cast<char*>(ppos_v.data()); }
  vec<uint8_t> get_npos() { return npos_v; }
  vec<uint8_t> get_ppos() { return ppos_v; }
  uint8_t get_k() { return k; }
  uint8_t get_h() { return h; }
  uint32_t get_m() { return m; }

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  vec<std::pair<int8_t, int8_t>> glsh_v;
};

#endif
