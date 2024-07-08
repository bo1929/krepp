#ifndef _KEREMET_H
#define _KEREMET_H

#include "CLI11.hpp"
#include "common.hpp"
#include "lshf.hpp"
#include "query.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "tree.hpp"
#include <cstdint>

class Bkrmt
{
public:
  Bkrmt(CLI::App& sub_build);
  void set_lshashf();
  void build_library();
  void save_metadata();
  void read_input_file();
  void parse_newick_tree();
  void build_for_subtree(node_sptr_t nd, DynHT& dt);

private:
  uint8_t k = 29;
  uint8_t w = k + 3;
  uint8_t h = 13;
  uint32_t m = 2;
  uint32_t r = 1;
  bool frac = false;
  uint32_t nrows = pow(2, 2 * h - 1);
  tree_sptr_t ref_tree = nullptr;
  record_sptr_t record = nullptr;
  lshf_sptr_t lshashf = nullptr;
  std::filesystem::path nwk_path;
  std::filesystem::path input_path;
  std::filesystem::path library_dir;
  std::string suffix;
  std::unordered_map<std::string, std::string> name_to_input;
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
};

class Library
{
public:
  Library(std::filesystem::path library_dir)
    : library_dir(library_dir){};
  void add_partial_tree(std::string suffix);
  void add_partial_flatht(std::string suffix);
  void add_partial_crecord(std::string suffix);
  lshf_sptr_t get_lshashf() { return lshashf; }
  std::vector<cmer_t>::const_iterator begin(uint32_t rix);
  std::vector<cmer_t>::const_iterator end(uint32_t rix);
  bool partial_exists(uint32_t rix) { return r_to_flatht.find(rix % m) != r_to_flatht.end(); }

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshashf = nullptr;
  crecord_sptr_t crecord = nullptr;
  std::filesystem::path library_dir;
  std::unordered_map<uint32_t, flatht_sptr_t> r_to_flatht;
};

#endif
