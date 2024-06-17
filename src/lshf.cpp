#include "lshf.hpp"

LSHF::LSHF(vec<uint8_t> &ppos_v, vec<uint8_t> &npos_v)
    : ppos_v(ppos_v), npos_v(npos_v) {
  std::sort(ppos_v.begin(), ppos_v.end(), std::greater<uint8_t>());
  std::vector<int8_t> v;
  std::vector<int8_t> g;
  int8_t lp = 31;
  int8_t jp = 0;
  for (int8_t j = 0; j < ppos_v.size(); j++) {
    if (j == 0) {
      v.push_back((lp - ppos_v[j]) * 2);
      lp = ppos_v[j];
      jp += 2;
    } else if ((ppos_v[j - 1] - ppos_v[j]) != 1) {
      v.push_back((lp - ppos_v[j]) * 2);
      lp = ppos_v[j];
      g.push_back(jp);
      jp = 2;
    } else {
      jp += 2;
    }
  }
  g.push_back(jp);
  v.push_back(-1);
  g.push_back(-1);
  glsh_v.resize(v.size() < g.size() ? v.size() : g.size());
  for (unsigned int i = 0; i < glsh_v.size(); i++) {
    glsh_v[i] = std::make_pair(v[i], g[i]);
  }
}

uint32_t LSHF::compute_hash(uint64_t enc_bp) {
  uint64_t res = 0;
  unsigned int i = 0;
  while (glsh_v[i].first != -1) {
    enc_bp = enc_bp << glsh_v[i].first;
    asm("shld %b3, %2, %0"
        : "=rm"(res)
        : "0"(res), "r"(enc_bp), "ic"(glsh_v[i].second)
        : "cc");
    i++;
  }
  return static_cast<uint32_t>(res);
}

void LSHF::drop_ppos_encoding(uint64_t enc64_bp, uint64_t enc64_lr,
                              uint32_t &enc32_bp, uint32_t &enc32_lr) {
  enc32_bp = 0;
  enc32_lr = 0;
  for (uint i = npos_v.size() - 1; i >= 0; --i) {
    enc32_lr += static_cast<uint32_t>((enc64_lr >> npos_v[i]) & 1);
    enc32_lr += static_cast<uint32_t>((enc64_lr >> (npos_v[i] + 32)) & 1) << 16;
    enc32_bp += static_cast<uint32_t>((enc64_bp >> (npos_v[i] * 2)) & 3);
    enc32_lr <<= 1 * (i & 0x00000001);
    enc32_bp <<= 2 * (i & 0x00000001);
  }
}
