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
  void get_random_positions();
  void set_lshf();
  uint8_t get_k();
  uint8_t get_h();
  uint32_t get_m();
  char* npos_data();
  char* ppos_data();
  vec<uint8_t> get_npos();
  vec<uint8_t> get_ppos();
  uint32_t compute_hash(uint64_t enc_bp);
  uint32_t drop_ppos_lr(uint64_t enc64_lr);
  uint32_t drop_ppos_bp(uint64_t enc64_bp);
  bool check_compatible(lshf_sptr_t lshf);

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  vec<std::pair<int8_t, int8_t>> glsh_v;
};

#endif
