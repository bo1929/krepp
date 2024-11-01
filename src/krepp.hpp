#ifndef _KREPP_H
#define _KREPP_H

#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include <CLI.hpp>

class Bkrepp
{
public:
  Bkrepp(CLI::App& sub_build);
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
  dynht_sptr_t root_dynht = nullptr;
  lshf_sptr_t lshf = nullptr;
  tree_sptr_t tree = nullptr;
  std::string suffix;
  std::filesystem::path nwk_path;
  std::filesystem::path input_path;
  std::filesystem::path library_dir;
  parallel_flat_phmap<std::string, std::string> name_to_path;
};

class Qkrepp
{
public:
  Qkrepp(CLI::App& sub_query);
  void load_library();
  void place_sequences();
  void estimate_distances();
  static void simulate_hdhistllh();

private:
  std::filesystem::path output_dir = "./";
  double min_gamma = 0.5; // TODO: Perhaps remove this?
  uint32_t hdist_th = 5;
  library_sptr_t library = nullptr;
  std::filesystem::path query_path;
  std::filesystem::path library_dir;
};

#endif
