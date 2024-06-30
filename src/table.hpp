#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "refseq.hpp"
#include "record.hpp"

#define NONSTD_UNION
// TODO: Choose one, one of the options has a bug.

class DynTable
{
  friend class FlatTable;

public:
  DynTable(uint32_t nrows, record_sptr_t record)
    : nrows(nrows)
    , record(record)
    , nkmers(0)
  {}
  void print_info();
  void clear_rows();
  void make_unique();
  void sort_columns();
  void ensure_sorted_columns();
  void prune_columns(size_t max_size);
  void union_table(DynTable& source);
  void update_size_hist();
  void update_nkmers();
  void reserve() { mer_vvec.reserve(nrows); }
  void fill_table(refseq_sptr_t refseq);
  uint64_t get_nkmers() { return nkmers; }
#ifdef NONSTD_UNION
  static void union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v, record_sptr_t record);
#else
  static void
  union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v, record_sptr_t record, bool in_place);
#endif
  static bool comp_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding < right.encoding;
  }
  static bool eq_encoding(const mer_t& left, const mer_t& right)
  {
    return left.encoding == right.encoding;
  }

private:
  uint32_t nrows;
  uint64_t nkmers;
  vvec<mer_t> mer_vvec;
  record_sptr_t record;
  std::map<size_t, uint32_t> size_hist;
};

class FlatTable
{
  friend class DynTable;

public:
  FlatTable(DynTable& source);
  FlatTable(std::filesystem::path library_dir, std::string suffix);
  void save(std::filesystem::path library_dir, std::string suffix);

private:
  uint32_t nrows;
  uint64_t nkmers;
  vec<mer_t> mer_v;
  vec<inc_t> inc_v;
};

#endif
