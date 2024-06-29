#ifndef _LSHF_H
#define _LSHF_H

#include "common.hpp"

class LSHF {
public:
  LSHF(vec<uint8_t> &ppos_v, vec<uint8_t> &npos_v);
  void drop_ppos_enc(uint64_t enc64_bp, uint64_t enc64_lr,
                          uint32_t &enc32_bp, uint32_t &enc32_lr);
  uint32_t compute_hash(uint64_t enc_bp);

private:
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  vec<std::pair<int8_t, int8_t>> glsh_v;
};

#endif
