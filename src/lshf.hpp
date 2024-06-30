#ifndef _LSHF_H
#define _LSHF_H

#include "common.hpp"
#include "refseq.hpp"

class LSHF
{
  friend class RefSeq;

public:
  LSHF(uint8_t k, uint8_t h, uint32_t m, uint32_t r, bool frac);
  void get_random_positions();
  void drop_ppos_enc(uint64_t enc64_bp, uint64_t enc64_lr, uint32_t& enc32_bp, uint32_t& enc32_lr);
  uint32_t compute_hash(uint64_t enc_bp);

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  uint32_t r;
  bool frac;
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  vec<std::pair<int8_t, int8_t>> glsh_v;
};

#endif
