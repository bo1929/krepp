#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "record.hpp"
#include "rqseq.hpp"

class SDynHT
{
  friend class SFlatHT;

public:
  void make_unique();
  void sort_columns();
  void fill_table(uint32_t nrows, rseq_sptr_t rs);
  uint64_t get_nkmers() { return nkmers; }

protected:
  uint64_t nkmers = 0;
  vvec<enc_t> enc_vvec;
};

class DynHT
{
  friend class FlatHT;

public:
  DynHT(uint32_t nrows, tree_sptr_t tree, record_sptr_t record)
    : nrows(nrows)
    , tree(tree)
    , record(record)
  {
  }
  DynHT()
    : tree(nullptr)
    , record(nullptr)
  {
  }
  void print_info();
  void clear_rows();
  void make_unique();
  void sort_columns();
  void update_nkmers();
  void update_size_hist();
  void ensure_sorted_columns();
  void fill_table(sh_t sh, rseq_sptr_t rqseq);
  void prune_columns(size_t max_size);
  void union_table(dynht_sptr_t source);
  void reserve() { mer_vvec.reserve(nrows); }
  uint64_t get_nkmers() { return nkmers; }
  tree_sptr_t get_tree() { return tree; }
  record_sptr_t get_record() { return record; }
  void set_tree(tree_sptr_t source) { tree = source; }
  void set_record(record_sptr_t source) { record = source; }
  void union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v);
  cmer_t conv_mer_cmer(mer_t x) { return std::make_pair(x.encoding, record->map_compact(x.sh)); }
  static bool comp_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding < right.encoding;
  }
  static bool eq_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding == right.encoding;
  }

private:
  uint64_t nkmers = 0;
  uint32_t nrows = 0;
  vvec<mer_t> mer_vvec;
  tree_sptr_t tree = nullptr;
  record_sptr_t record = nullptr;
  flat_phmap<uint64_t, uint32_t> size_hist;
};

class SFlatHT
{
  friend class SDynHT;

public:
  SFlatHT(sdynht_sptr_t source);
  SFlatHT() {};
  ~SFlatHT()
  {
    inc_v.clear();
    enc_v.clear();
  }
  void save(std::ofstream& sketch_stream);
  void load(std::ifstream& sketch_stream);
  std::vector<enc_t>::const_iterator bucket_start(uint32_t rix)
  {
    if (rix) {
      return std::next(enc_v.begin(), inc_v[rix - 1]);
    } else {
      return enc_v.begin();
    }
  }
  std::vector<enc_t>::const_iterator bucket_next(uint32_t rix)
  {
    if (rix < inc_v.size()) {
      return std::next(enc_v.begin(), inc_v[rix]);
    } else {
      return enc_v.end();
    }
  }

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<enc_t> enc_v;
};

class FlatHT
{
  friend class DynHT;

public:
  FlatHT(dynht_sptr_t source);
  FlatHT(tree_sptr_t tree, crecord_sptr_t crecord)
    : tree(tree)
    , crecord(crecord) {};
  ~FlatHT()
  {
    inc_v.clear();
    cmer_v.clear();
    crecord.reset();
    tree.reset();
  }
  void load(std::ifstream& mer_stream, std::ifstream& inc_stream);
  void save(std::ofstream& mer_stream, std::ofstream& inc_stream);
  void set_crecord(crecord_sptr_t source) { crecord = source; }
  void set_tree(tree_sptr_t source) { tree = source; }
  uint64_t get_nkmers() { return nkmers; }
  tree_sptr_t get_tree() { return tree; }
  crecord_sptr_t get_crecord() { return crecord; }
  inc_t get_inc(uint32_t rix) { return inc_v[rix]; }
  std::vector<cmer_t>::const_iterator bucket_start(uint32_t rix)
  {
    if (rix) {
      return std::next(cmer_v.begin(), inc_v[rix - 1]);
    } else {
      return cmer_v.begin();
    }
  }
  std::vector<cmer_t>::const_iterator bucket_next(uint32_t rix)
  {
    if (rix < inc_v.size()) {
      return std::next(cmer_v.begin(), inc_v[rix]);
    } else {
      return cmer_v.end();
    }
  }
  void display_info(uint32_t r);

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<cmer_t> cmer_v;
  tree_sptr_t tree = nullptr;
  crecord_sptr_t crecord = nullptr;
};

#endif
