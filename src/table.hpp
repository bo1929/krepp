#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "record.hpp"
#include "rqseq.hpp"

#define NONSTD_UNION
// TODO: Choose one, one of the options has a bug.

class DynHT
{
  friend class FlatHT;

public:
  DynHT(uint32_t nrows, tree_sptr_t tree, record_sptr_t record)
    : nrows(nrows)
    , tree(tree)
    , record(record)
    , nkmers(0)
  {}
  void print_info();
  void clear_rows();
  void make_unique();
  void sort_columns();
  void update_nkmers();
  void update_size_hist();
  void ensure_sorted_columns();
  void union_table(DynHT& source);
  void fill_table(rseq_sptr_t rqseq);
  void prune_columns(size_t max_size);
  void reserve() { mer_vvec.reserve(nrows); }
#ifdef NONSTD_UNION // TODO: Benchmark for static vs non-static, and std non-std, in-place etc.
  void union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v);
#else
  void union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v, bool in_place);
#endif
  static bool comp_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding < right.encoding;
  }
  static bool eq_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding == right.encoding;
  }
  uint64_t get_nkmers() { return nkmers; }
  cmer_t conv_mer_cmer(mer_t x) { return std::make_pair(x.encoding, record->map_compact(x.shash)); }

private:
  uint32_t nrows;
  uint64_t nkmers;
  tree_sptr_t tree;
  record_sptr_t record;
  vvec<mer_t> mer_vvec;
  std::map<size_t, uint32_t> size_hist;
};

class FlatHT
{
  friend class DynHT;

public:
  FlatHT(DynHT& source);
  FlatHT(tree_sptr_t tree, crecord_sptr_t crecord)
    : tree(tree)
    , crecord(crecord){};
  void load(std::filesystem::path library_dir, std::string suffix);
  void save(std::filesystem::path library_dir, std::string suffix);
  void set_crecord(crecord_sptr_t source) { crecord = source; }
  void set_tree(tree_sptr_t source) { tree = source; }
  crecord_sptr_t get_crecord() { return crecord; }
  tree_sptr_t get_tree() { return tree; }
  std::vector<cmer_t>::const_iterator begin() { return cmer_v.begin(); }
  std::vector<cmer_t>::const_iterator end() { return cmer_v.end(); }
  std::vector<cmer_t>::const_iterator at(uint32_t rix)
  {
    return std::next(cmer_v.begin(), inc_v[rix]);
  }
  inc_t get_inc(uint32_t rix) { return inc_v[rix]; }

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<cmer_t> cmer_v;
  tree_sptr_t tree = nullptr;
  crecord_sptr_t crecord = nullptr;
};

#endif
