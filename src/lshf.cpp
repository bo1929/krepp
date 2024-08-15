#include "lshf.hpp"

LSHF::LSHF(uint8_t k, uint8_t h, uint32_t m)
  : k(k)
  , h(h)
  , m(m)
{
  get_random_positions();
  set_lshf();
}

void LSHF::set_lshf()
{
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

uint32_t LSHF::compute_hash(uint64_t enc_bp)
{
  uint64_t res = 0;
  unsigned int i = 0;
  while (glsh_v[i].first != -1) {
    enc_bp = enc_bp << glsh_v[i].first;
    asm("shld %b3, %2, %0" : "=rm"(res) : "0"(res), "r"(enc_bp), "ic"(glsh_v[i].second) : "cc");
    i++;
  }
  return static_cast<uint32_t>(res);
}

uint32_t LSHF::drop_ppos_lr(uint64_t enc64_lr)
{
  uint32_t enc32_lr = 0;
  for (int i = npos_v.size() - 1; i >= 0; --i) {
    enc32_lr += static_cast<uint32_t>((enc64_lr >> npos_v[i]) & 1);
    enc32_lr += static_cast<uint32_t>((enc64_lr >> (npos_v[i] + 32)) & 1) << 16;
    enc32_lr <<= 1 * (i & 0x00000001);
  }
  return enc32_lr;
}

uint32_t LSHF::drop_ppos_bp(uint64_t enc64_bp)
{
  uint32_t enc32_bp = 0;
  for (int i = npos_v.size() - 1; i >= 0; --i) {
    enc32_bp += static_cast<uint32_t>((enc64_bp >> (npos_v[i] * 2)) & 3);
    enc32_bp <<= 2 * (i & 0x00000001);
  }
  return enc32_bp;
}

void LSHF::get_random_positions()
{
  uint8_t n;
  assert(h <= 16);
  assert(h < k);
  std::uniform_int_distribution<uint8_t> distrib(0, k - 1);
  for (uint8_t m = 0; m < h; m++) {
    n = distrib(gen);
    if (std::count(ppos_v.begin(), ppos_v.end(), n)) {
      m -= 1;
    } else {
      ppos_v.push_back(n);
    }
  }
  std::sort(ppos_v.begin(), ppos_v.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < k; ++i) {
    if (i != ppos_v[ix_pos])
      npos_v.push_back(i);
    else
      ix_pos++;
  }
  std::sort(ppos_v.begin(), ppos_v.end(), std::greater<uint8_t>());
}

LSHF::LSHF(uint32_t m, vec<uint8_t> ppos_v, vec<uint8_t> npos_v)
  : m(m)
  , ppos_v(ppos_v)
  , npos_v(npos_v)
{
  k = npos_v.size() + ppos_v.size();
  h = ppos_v.size();
  set_lshf();
}

bool LSHF::check_compatible(lshf_sptr_t lshf)
{
  if (!lshf)
    return true;
  if (!((lshf->m == m) && (lshf->h == h) && (lshf->k == k) && (lshf->npos_v == npos_v) &&
        (lshf->ppos_v == ppos_v))) {
    std::cout << "m: " << m << "/" << lshf->m << std::endl;
    std::cout << "h: " << h << "/" << lshf->h << std::endl;
    std::cout << "k: " << k << "/" << lshf->k << std::endl;
    std::cout << "ppos_v:";
    for (uint8_t i = 0; i < h; ++i) {
      std::cout << " " << static_cast<uint32_t>(ppos_v[i]) << "/"
                << static_cast<uint32_t>(lshf->ppos_v[i]);
    }
    std::cout << std::endl;
    std::cout << "npos_v:";
    for (uint8_t i = 0; i < k - h; ++i) {
      std::cout << " " << static_cast<uint32_t>(npos_v[i]) << "/"
                << static_cast<uint32_t>(lshf->npos_v[i]);
    }
    std::cout << std::endl;
    return false;
  } else {
    return true;
  }
}
