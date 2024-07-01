#ifndef _KEREMET_H
#define _KEREMET_H

#include "CLI11.hpp"
#include "common.hpp"
#include "lshf.hpp"
#include "record.hpp"
#include "refseq.hpp"
#include "table.hpp"
#include "tree.hpp"

class Bkrmt
{

public:
  Bkrmt(CLI::App& sub_build);
  void set_hash_func();
  void read_input_file();
  void parse_newick_tree();
  void build_for_subtree(node_sptr_t nd, DynTable& dt);
  void build_library();

private:
  uint8_t k = 29;
  uint8_t w = k + 3;
  uint8_t h = 13;
  uint32_t m = 2;
  uint32_t r = 1;
  uint32_t nrows = pow(2, 2 * h - 1);
  bool frac = false;
  std::string library_dir;
  std::string nwk_filepath;
  std::string input_filepath;
  std::unordered_map<std::string, std::string> name_to_gpath;
  std::string suffix;
  tree_sptr_t ref_tree = nullptr;
  lshf_sptr_t hash_func = nullptr;
};

class Qkrmt
{

public:
  Qkrmt(CLI::App& sub_query);

private:
  std::vector<std::string> library_dir_v;
  std::string output_dir = "./";
  std::string query_file;
};

#endif
