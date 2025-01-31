#ifndef _KREPP_H
#define _KREPP_H

#include "common.hpp"
#include "index.hpp"
#include "sketch.hpp"
#include "lshf.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "compare.hpp"
#include "record.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include <CLI.hpp>

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void save_configuration(std::ofstream& cfg_stream);

protected: // TODO: adaptive defaults for sketch versus index.
  uint8_t w = k + 3;
  uint8_t k = 29;
  uint8_t h = 14;
  uint32_t m = 2;
  uint32_t r = 1;
  bool frac = false;
  uint32_t nrows = pow(2, 2 * h - 1);
  lshf_sptr_t lshf = nullptr;
};

class SketchSingle : public BaseLSH
{
public:
  SketchSingle(CLI::App& sub_ss);
  void create_sketch();
  void save_sketch();

private:
  double rho;
  sflatht_sptr_t sketch_sflatht = nullptr;
  std::filesystem::path input_path;
  std::filesystem::path output_path;
};

class TargetSketch
{
public:
  void load_sketch();

protected:
  sketch_sptr_t sketch = nullptr;
  std::filesystem::path sketch_path;
};

class CompareSketch : public TargetSketch
{
public:
  CompareSketch(CLI::App& sub_sscomp);
  void estimate_distances();

private:
  uint32_t hdist_th = 4;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  std::filesystem::path query_path;
};

class IndexMultiple : public BaseLSH
{
public:
  IndexMultiple(CLI::App& sub_im);
  void obtain_build_tree();
  void save_metadata();
  void save_index();
  void build_index();
  void build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht);
  void read_input_file();

private:
  std::string suffix;
  tuint_t build_count = 0;
  std::filesystem::path input_path;
  std::filesystem::path index_dir;
  std::filesystem::path nwk_path;
  tree_sptr_t tree = nullptr;
  flatht_sptr_t root_flatht = nullptr;
  parallel_flat_phmap<std::string, std::string> name_to_path;
};

class TargetIndex
{
public:
  void load_index();
  bool check_wtree() { return index->check_wtree(); }

protected:
  index_sptr_t index = nullptr;
  std::filesystem::path index_dir;
};

class QueryIndex : public TargetIndex
{
public:
  QueryIndex(CLI::App& sub_imquery);
  void estimate_distances();
  void place_sequences();
  void header_dreport(strstream& dreport_stream);
  void begin_jplace(strstream& jplace_stream);
  void end_jplace(strstream& jplace_stream);

private:
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
  uint32_t tau = 3;
  bool no_filter = false;
  std::filesystem::path query_path;
};

class InfoIndex : public TargetIndex
{
public:
  InfoIndex(CLI::App& sub_iminfo);
  void display_info();
};

#endif
