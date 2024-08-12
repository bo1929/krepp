#ifndef _KEREMET_H
#define _KEREMET_H

#include "CLI11.hpp"
#include "common.hpp"
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
  void set_lshashf();
  void save_metadata();
  void read_input_file();
  void initialize_record();
  void parse_newick_tree();
  void save_library(DynHT& root_dynht);
  void build_library(DynHT& root_dynht);
  void build_for_subtree(node_sptr_t nd, DynHT& dynht);

private:
  uint32_t nrows = pow(2, 2 * h - 1);
  uint8_t k = 29;
  uint8_t w = k + 3;
  uint8_t h = 13;
  uint32_t m = 2;
  uint32_t r = 1;
  bool frac = false;
  std::string suffix;
  tuint_t build_count = 0;
  std::filesystem::path nwk_path;
  std::filesystem::path input_path;
  std::filesystem::path library_dir;
  tree_sptr_t tree = nullptr;
  record_sptr_t record = nullptr;
  lshf_sptr_t lshashf = nullptr;
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

class Library
{
public:
  Library(std::filesystem::path library_dir)
    : library_dir(library_dir){};
  void add_partial_tree(std::string suffix);
  void add_partial_flatht(std::string suffix);
  void add_partial_crecord(std::string suffix);
  std::vector<cmer_t>::const_iterator get_first(uint32_t rix);
  std::vector<cmer_t>::const_iterator get_next(uint32_t rix);
  bool check_partial(uint32_t rix) { return r_to_flatht.contains(rix % m); }
  flatht_sptr_t get_flatht_sptr(uint32_t rix) { return r_to_flatht[rix % m]; };
  crecord_sptr_t get_crecord() { return crecord; }
  lshf_sptr_t get_lshashf() { return lshashf; }
  tree_sptr_t get_tree() { return tree; }

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  uint32_t nrows;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshashf = nullptr;
  crecord_sptr_t crecord = nullptr;
  std::filesystem::path library_dir;
  parallel_flat_phmap<uint32_t, flatht_sptr_t> r_to_flatht;
};

#endif
