#ifndef _SKETCH_H
#define _SKETCH_H

#include "common.hpp"
#include "lshf.hpp"
#include "table.hpp"

typedef std::vector<enc_t>::const_iterator vec_enc_it;

class Sketch
{
public:
  Sketch(std::filesystem::path sketch_path)
    : sketch_path(sketch_path) {};
  void load_full_sketch();
  void make_rho_partial();
  std::pair<vec_enc_it, vec_enc_it> bucket_indices(uint32_t rix);
  bool check_partial(uint32_t rix)
  {
    uint32_t rix_res = rix % m;
    return (frac && (rix_res <= r)) || (rix_res == r);
  }
  sflatht_sptr_t get_sflatht_sptr() { return sflatht; };
  lshf_sptr_t get_lshf() { return lshf; }
  double get_rho() { return rho; }

private:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  bool frac;
  uint32_t r;
  uint32_t m;
  double rho;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
  sflatht_sptr_t sflatht = nullptr;
  std::filesystem::path sketch_path;
};

#endif
