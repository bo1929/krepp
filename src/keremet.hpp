#ifndef _KEREMET_H
#define _KEREMET_H

#include "CLI11.hpp"
#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "table.hpp"

class Bkrmt
{
public:
  Bkrmt(CLI::App& sub_build);
  void set_lshf();
  void save_metadata();
  void read_input_file();
  void parse_newick_tree();
  void save_library();
  void build_library();
  void build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht);

private:
  uint8_t k = 29;
  uint8_t w = k + 3;
  uint8_t h = 13;
  uint32_t m = 2;
  uint32_t r = 1;
  bool frac = false;
  tuint_t build_count = 0;
  uint32_t nrows = pow(2, 2 * h - 1);
  std::string suffix;
  std::filesystem::path nwk_path;
  std::filesystem::path input_path;
  std::filesystem::path library_dir;
  lshf_sptr_t lshf = nullptr;
  tree_sptr_t tree = nullptr;
  dynht_sptr_t root_dynht = nullptr;
  parallel_flat_phmap<std::string, std::string> name_to_input;
};

class Pkrmt
{
public:
  Pkrmt(CLI::App& sub_query);
  void load_library();
  void place_sequences();

private:
  std::filesystem::path output_dir = "./";
  std::filesystem::path library_dir;
  std::filesystem::path query_path;
  library_sptr_t library = nullptr;
  uint32_t max_hdist = 5;
  float min_covpos = 0.5;
};

#endif
