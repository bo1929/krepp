#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "refseq.hpp"
#include "subset.hpp"
#include <cstdint>

#define NONSTD_UNION

class DynTable {
public:
  DynTable(uint32_t nrows, record_sptr_t record) : nrows(nrows), record(record), nkmers(0) {}
  void clear_rows();
  void make_unique();
  void sort_columns();
  void ensure_sorted_columns();
  void prune_columns(size_t max_size);
  void union_table(dyntable_sptr_t source);
  void update_size_hist();
#ifdef NONSTD_UNION
  static void union_row(vec<mer_t> &dest_v, vec<mer_t> &source_v,
                        record_sptr_t record);
#else
  static void union_row(vec<mer_t> &dest_v, vec<mer_t> &source_v,
                        record_sptr_t record, bool in_place);
#endif
  void fill_table(refseq_sptr_t refseq);
  // void convert_table(FlatTable &dest);
  static bool comp_encoding(const mer_t &left, const mer_t &right) {
    return left.encoding < right.encoding;
  }
  static bool eq_encoding(const mer_t &left, const mer_t &right) {
    return left.encoding == right.encoding;
  }

private:
  uint32_t nrows;
  uint64_t nkmers;
  vvec<mer_t> table;
  record_sptr_t record;
  std::map<size_t, uint32_t> size_hist;
};

class FlatTable {
public:
  FlatTable(){};

private:
};

#endif
