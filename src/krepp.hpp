#ifndef _KREPP_H
#define _KREPP_H

#include "common.hpp"
#include "index.hpp"
#include "sketch.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "seek.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include <CLI.hpp>

const auto url_validator = CLI::Validator(
  [](std::string& input) {
    const std::regex url_regexp = std::regex(
      R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
    if (std::regex_match(input, url_regexp)) {
      return std::string("");
    } else {
      return "Given URL is not valid: " + input;
    }
  },
  "URL",
  "URL validator");

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void save_configuration(std::ofstream& cfg_stream);
  void set_sketch_defaults()
  {
    k = 26;
    w = k + 6;
    h = 10;
    m = 5;
    r = 1;
    frac = true;
    nrows = pow(2, 2 * h - 1);
  }
  void set_index_defaults()
  {
    k = 29;
    w = k + 6;
    h = 13;
    m = 5;
    r = 1;
    frac = true;
    nrows = pow(2, 2 * h - 1);
  }
  bool validate_configuration()
  {
    bool is_valid = true;
    if (is_valid = w < k) {
      std::cerr << "The minimum minimizer window size (-w) is k (-k)." << std::endl;
    }
    if (is_valid = h < 3) {
      std::cerr << "The minimum number of LSH positions (-h) is 3." << std::endl;
    }
    if (is_valid = h > 15) {
      std::cerr << "The maximum number of LSH positions (-h) is 15." << std::endl;
    }
    if (is_valid = k > 31) {
      std::cerr << "The maximum allowed k-mer length (-k) is 31." << std::endl;
    }
    if (is_valid = k < 19) {
      std::cerr << "The minimum allowed k-mer length (-k) is 19." << std::endl;
    }
    if (is_valid = (k - h) > 16) {
      std::cerr << "For compact k-mer encodings, h must be >= k-16." << std::endl;
    }
    return !is_valid;
  }

protected:
  uint8_t w;
  uint8_t k;
  uint8_t h;
  bool frac;
  uint32_t m;
  uint32_t r;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
};

class TargetSketch
{
public:
  void load_sketch();

protected:
  sketch_sptr_t sketch = nullptr;
  std::filesystem::path sketch_path;
};

class TargetIndex
{
public:
  void load_index();
  void ensure_wbackbone();

protected:
  index_sptr_t index = nullptr;
  std::filesystem::path index_dir;
};

class SketchSingle : public BaseLSH
{
public:
  SketchSingle(CLI::App& sc);
  void create_sketch();
  void save_sketch();

private:
  double rho;
  std::string input;
  std::filesystem::path sketch_path;
  sflatht_sptr_t sketch_sflatht = nullptr;
};

class IndexMultiple : public BaseLSH
{
public:
  IndexMultiple(CLI::App& sc);
  void obtain_build_tree();
  void read_input_file();
  void save_index();
  void build_index();
  void build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht);

private:
  std::string suffix;
  tuint_t build_count = 0;
  vec<std::string> names_v;
  std::filesystem::path input;
  std::filesystem::path index_dir;
  std::filesystem::path nwk_path;
  tree_sptr_t tree = nullptr;
  flatht_sptr_t root_flatht = nullptr;
  parallel_flat_phmap<std::string, std::string> name_to_path;
};

class QuerySketch : public TargetSketch
{
public:
  QuerySketch(CLI::App& sc);
  void seek_sequences();
  void header_dreport(strstream& dreport_stream);
  uint32_t get_total_qseq() { return total_qseq; }

private:
  std::string query;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
  uint64_t total_qseq = 0;
};

class QueryIndex : public TargetIndex
{
public:
  QueryIndex(CLI::App& sc);
  void init_sc_place(CLI::App& sc);
  void init_sc_dist(CLI::App& sc);
  void estimate_distances();
  void place_sequences();
  void header_dreport(strstream& dreport_stream);
  void begin_jplace(strstream& jplace_stream);
  void end_jplace(strstream& jplace_stream);
  uint32_t get_total_qseq() { return total_qseq; }

private:
  std::string query;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t tau = 2;
  uint32_t hdist_th = 4;
  double dist_max = 0.2;
  bool no_filter = true;
  bool multi = true;
  uint64_t total_qseq = 0;
};

class InfoIndex : public TargetIndex
{
public:
  InfoIndex(CLI::App& sc);
  void display_info();
};

#endif
