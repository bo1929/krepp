#ifndef _INDEX_H
#define _INDEX_H

#include "common.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "record.hpp"
#include "table.hpp"
#include <optional>

typedef std::vector<cmer_t>::const_iterator vec_cmer_it;

class Index
{
public:
  Index(std::filesystem::path index_dir)
    : index_dir(index_dir) {};
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

struct IndexConfig
{
  std::filesystem::path input;     // TSV of reference ID -> path/URL, or one FASTA/FASTQ
  std::filesystem::path index_dir; // directory to write the index into
  std::filesystem::path nwk_path;  // guide tree; empty means generate one
  uint8_t k = 29;
  std::optional<uint8_t> w; // unset: k + 6
  std::optional<uint8_t> h; // unset: k - 16
  uint32_t m = 4;
  uint32_t r = 1;
  bool frac = true;
  uint32_t sdust_t = 0;
  uint32_t sdust_w = 0;
};

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void save_configuration(std::ofstream& cfg_stream);
  void set_sketch_defaults()
  {
    k = 25;
    w = k + 6;
    h = 10;
    m = 4;
    r = 1;
    frac = true;
    nrows = pow(2, 2 * h - 1);
    sdust_t = 0;
    sdust_w = 0;
  }
  bool validate_configuration()
  {
    bool is_invalid = true;
    if ((is_invalid = (w < k))) {
      error_exit("The minimum minimizer window size (-w) is k (-k).");
    }
    if ((is_invalid = (h < 9))) {
      error_exit("The minimum number of LSH positions (-h) is 9.");
    }
    if ((is_invalid = (h > 15))) {
      error_exit("The maximum number of LSH positions (-h) is 15.");
    }
    if ((is_invalid = (k > 31))) {
      error_exit("The maximum allowed k-mer length (-k) is 31.");
    }
    if ((is_invalid = (k < 19))) {
      error_exit("The minimum allowed k-mer length (-k) is 19.");
    }
    if ((is_invalid = ((k - h) > 16))) {
      error_exit("For compact k-mer encodings, h must be >= k-16.");
    }
    if ((sdust_t != 0) && (sdust_w != 0)) {
      std::cerr << "Setting --sdust-w and --sdust-t to >0 will enable dustmasker." << std::endl;
      std::cerr << "With dustmasker, krepp might fail to model subsampling and be slightly inaccurate." << std::endl;
    }
    return !is_invalid;
  }

protected:
  uint8_t w;
  uint8_t k;
  uint8_t h;
  bool frac;
  uint32_t m;
  uint32_t r;
  uint32_t nrows;
  uint32_t sdust_t;
  uint32_t sdust_w;
  lshf_sptr_t lshf = nullptr;
};

class IndexMultiple : public BaseLSH
{
public:
  IndexMultiple(const IndexConfig& config);
  void obtain_build_tree();
  void read_input_file();
  void save_index();
  void build_index();
  void index_sequences();
  void index_files();
  void build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht, ErrorRelay& relay);
  void save_info(std::ofstream& info_stream);

private:
  std::string suffix;
  tuint_t build_count = 0;
  vec<std::string> names_v;
  vec<std::string> fastx_names;
  vec<uint64_t> fastx_offsets;
  std::filesystem::path input;
  std::filesystem::path index_dir;
  std::filesystem::path nwk_path;
  bool per_sequence = false;
  tree_sptr_t tree = nullptr;
  flatht_sptr_t root_flatht = nullptr;
  parallel_flat_phmap<std::string, std::string> name_to_path;
};

#endif
