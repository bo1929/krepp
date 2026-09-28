#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "filemap.hpp"
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
  void update_size_hist();
  void fill_table(sh_t sh, rseq_sptr_t rqseq, bool curr = false);
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
  static bool comp_encoding(const mer_t& left, const mer_t& right) { return left.encoding < right.encoding; }
  static bool eq_encoding(const mer_t& left, const mer_t& right) { return left.encoding == right.encoding; }

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
  ~SFlatHT() = default;
  void save(std::ofstream& sketch_stream);
  size_t load(std::ifstream& sketch_stream);
  size_t load(std::ifstream& sketch_stream, const std::filesystem::path& path);
  const enc_t* bucket_start(uint32_t rix) const
  {
    if (rix) {
      return enc_v + inc_at(rix - 1);
    } else {
      return enc_v;
    }
  }
  const enc_t* bucket_next(uint32_t rix) const
  {
    if (rix < nrows) {
      return enc_v + inc_at(rix);
    } else {
      return enc_v + nkmers;
    }
  }

private:
  inc_t inc_at(uint32_t rix) const
  {
    inc_t value = 0;
    std::memcpy(&value, inc_bytes + static_cast<size_t>(rix) * sizeof(inc_t), sizeof(inc_t));
    return value;
  }

  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  const enc_t* enc_v = nullptr;
  const char* inc_bytes = nullptr;
  vec<inc_t> inc_owned;
  vec<enc_t> enc_owned;
  krepp::FileMap map;
};

class FlatHT
{
  friend class DynHT;

public:
  FlatHT(dynht_sptr_t source);
  FlatHT(tree_sptr_t tree, crecord_sptr_t crecord)
    : tree(tree)
    , crecord(crecord) {};
  ~FlatHT() = default;
  void load(const std::filesystem::path& mer_path, const std::filesystem::path& inc_path);
  void load(std::ifstream& mer_stream, std::ifstream& inc_stream);
  void save(std::ofstream& mer_stream, std::ofstream& inc_stream);
  void set_crecord(crecord_sptr_t source) { crecord = source; }
  void set_tree(tree_sptr_t source) { tree = source; }
  uint64_t get_nkmers() { return nkmers; }
  uint32_t get_nrows() { return nrows; }
  bool is_mapped() const { return mer_map.is_open() || inc_map.is_open(); }
  tree_sptr_t get_tree() { return tree; }
  crecord_sptr_t get_crecord() { return crecord; }
  inc_t get_inc(uint32_t rix) const { return inc_at(rix); }
  const cmer_t* bucket_data(uint32_t rix) const
  {
    if (rix == 0) return cmer_v;
    if (rix <= nrows) return cmer_v + inc_at(rix - 1);
    return cmer_v + nkmers;
  }
  const cmer_t* bucket_start(uint32_t rix) const { return bucket_data(rix); }
  const cmer_t* bucket_next(uint32_t rix) const { return bucket_data(rix + 1); }
  void display_info(std::ostream* output_stream, uint32_t r);

private:
  inc_t inc_at(uint32_t rix) const
  {
    inc_t value = 0;
    std::memcpy(&value, inc_bytes + static_cast<size_t>(rix) * sizeof(inc_t), sizeof(inc_t));
    return value;
  }
  void bind();

  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  const cmer_t* cmer_v = nullptr;
  const char* inc_bytes = nullptr;
  vec<cmer_t> cmer_owned;
  vec<inc_t> inc_owned;
  krepp::FileMap mer_map;
  krepp::FileMap inc_map;
  tree_sptr_t tree = nullptr;
  crecord_sptr_t crecord = nullptr;
};

#endif
