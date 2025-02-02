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
    frac = false;
    nrows = pow(2, 2 * h - 1);
  }
  void set_index_defaults()
  {
    k = 29;
    w = k + 6;
    h = 13;
    m = 5;
    r = 1;
    frac = false;
    nrows = pow(2, 2 * h - 1);
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
  std::string input_file;
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
  std::filesystem::path input_file;
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

private:
  std::string query_file;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t hdist_th = 4;
};

class QueryIndex : public TargetIndex
{
public:
  QueryIndex(CLI::App& sc);
  void estimate_distances();
  void place_sequences();
  void header_dreport(strstream& dreport_stream);
  void begin_jplace(strstream& jplace_stream);
  void end_jplace(strstream& jplace_stream);

private:
  std::string query_file;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  uint32_t tau = 3;
  uint32_t hdist_th = 4;
  bool no_filter = false;
};

class InfoIndex : public TargetIndex
{
public:
  InfoIndex(CLI::App& sc);
  void display_info();
};

#endif
