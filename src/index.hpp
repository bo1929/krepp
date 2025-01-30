#ifndef _INDEX_H
#define _INDEX_H

#include "common.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "record.hpp"
#include "table.hpp"

class Index
{
public:
  Index(std::filesystem::path index_dir)
    : index_dir(index_dir){};
  void add_partial_index(std::string suffix);
  void make_rho_partial();
  crecord_sptr_t get_crecord();
  std::vector<cmer_t>::const_iterator bucket_start();
  std::vector<cmer_t>::const_iterator bucket_next();
  bool set_partial(uint32_t rix)
  {
    rix_res = rix % m;
    offset = rix / m;
    if (r_to_numerator[rix_res] > 1) {
      offset = offset * r_to_numerator[rix_res] + rix_res;
    }
    return r_to_flatht.contains(rix_res);
  }
  flatht_sptr_t get_flatht_sptr(uint32_t rix) { return r_to_flatht[rix % m]; };
  lshf_sptr_t get_lshf() { return lshf; }
  tree_sptr_t get_tree() { return tree; }
  void display_info()
  {
    for (auto const& [key, val] : r_to_flatht) {
      val->display_info(key);
    }
  }

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  uint32_t nrows;
  uint32_t offset;
  uint32_t rix_res;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  std::filesystem::path index_dir;
  parallel_flat_phmap<uint32_t, flatht_sptr_t> r_to_flatht;
  parallel_flat_phmap<uint32_t, uint32_t> r_to_numerator;
};

#endif
