#ifndef _INDEX_H
#define _INDEX_H

#include "common.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "record.hpp"
#include "table.hpp"

typedef std::vector<cmer_t>::const_iterator vec_cmer_it;

class Index
{
public:
  Index(std::filesystem::path index_dir)
    : index_dir(index_dir){};
  void make_rho_partial();
  crecord_sptr_t get_crecord(uint32_t rix);
  void load_partial_tree(std::string suffix);
  void load_partial_index(std::string suffix);
  void generate_partial_tree(std::string suffix);
  void display_info(std::ostream* output_stream);
  std::pair<vec_cmer_it, vec_cmer_it> bucket_indices(uint32_t rix);
  lshf_sptr_t get_lshf() { return lshf; }
  tree_sptr_t get_tree() { return tree; }
  bool check_wbackbone() { return wbackbone; }
  bool check_partial(uint32_t rix) { return r_to_flatht.contains(rix % m); }
  flatht_sptr_t get_flatht_sptr(uint32_t rix) { return r_to_flatht[rix % m]; };

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  uint32_t nrows;
  bool wbackbone = false;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  std::filesystem::path index_dir;
  fparallel_flat_phmap<uint32_t, flatht_sptr_t> r_to_flatht;
  fparallel_flat_phmap<uint32_t, uint32_t> r_to_numerator;
  fparallel_flat_phmap<uint32_t, std::string> r_to_info;
};

#endif
