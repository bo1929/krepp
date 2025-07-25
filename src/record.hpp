#ifndef _RECORD_H
#define _RECORD_H

#include "MurmurHash3.hpp"
#include "common.hpp"
#include "phytree.hpp"

#define MMHSEED1 1
#define MMHSEED0 0

class Subset : public std::enable_shared_from_this<Subset>
{
  friend class Record;
  friend class CRecord;

private:
  sh_t sh = 0;
  sh_t ch = 0;
  tuint_t card = 0;
  sh_t nonce = 0;

public:
  Subset(sh_t sh, sh_t ch, tuint_t card, sh_t nonce = 0)
    : sh(sh)
    , ch(ch)
    , card(card)
    , nonce(nonce)
  {
  }
  subset_sptr_t getptr() { return shared_from_this(); }
  static inline sh_t get_singleton_sh(std::string& name)
  {
    sh_t sh = 0;
    uint32_t a1, a2;
    MurmurHash3_x86_32(name.c_str(), name.length(), MMHSEED0, &a1);
    MurmurHash3_x86_32(name.c_str(), name.length(), MMHSEED1, &a2);
    sh = sh | static_cast<uint64_t>(a1);
    sh = sh << 32;
    sh = sh | static_cast<uint64_t>(a2);
    return sh;
  }
  static inline sh_t rehash(sh_t sh)
  {
    uint32_t a1, a2;
    MurmurHash3_x86_32(&sh, sizeof(sh), MMHSEED0, &a1);
    MurmurHash3_x86_32(&sh, sizeof(sh), MMHSEED1, &a2);
    sh = 0;
    sh = sh | static_cast<uint64_t>(a1);
    sh = sh << 32;
    sh = sh | static_cast<uint64_t>(a2);
    return sh;
  }
};

class Record : public std::enable_shared_from_this<Record>
{
  friend class CRecord;
  friend class Subset;

public:
  Record(tree_sptr_t tree);
  Record(record_sptr_t source1, record_sptr_t source2);
  bool check_tree_collision();
  void rehash_tree();
  void make_compact();
  sh_t add_subset(sh_t sh1, sh_t sh2);
  void union_record(record_sptr_t source);
  void decode_sh(sh_t sh, vec<node_sptr_t>& subset_v);
  bool check_subset_collision(subset_sptr_t s, subset_sptr_t s1, subset_sptr_t s2);
  uint64_t get_size() { return sh_to_subset.size(); }
  record_sptr_t getptr() { return shared_from_this(); }
  void insert_rho(sh_t sh, double rho) { sh_to_rho[sh] = rho; }
  se_t map_compact(sh_t sh) { return sh_to_se[sh]; }
  tree_sptr_t get_tree() { return tree; }

private:
  tree_sptr_t tree = nullptr;
  parallel_flat_phmap<sh_t, se_t> sh_to_se = {};
  parallel_flat_phmap<sh_t, node_sptr_t> sh_to_node = {};
  parallel_flat_phmap<sh_t, subset_sptr_t> sh_to_subset = {};
  parallel_flat_phmap<sh_t, double> sh_to_rho = {};
};

class CRecord : public std::enable_shared_from_this<CRecord>
{
public:
  CRecord(tree_sptr_t tree);
  CRecord(record_sptr_t record);
  ~CRecord()
  {
    se_to_pse.clear();
    se_to_rho.clear();
    tree.reset();
  }
  void print_info();
  void apply_rho_coef(double coef);
  void decode_se(se_t se, vec<node_sptr_t>& subset_v);
  void load(std::ifstream& crecord_stream);
  void save(std::ofstream& crecord_stream);
  se_t get_nsubsets() { return nsubsets; }
  crecord_sptr_t getptr() { return shared_from_this(); }
  pse_t get_pse(se_t se) const { return se_to_pse[se]; }
  double get_rho(se_t se) const { return se_to_rho[se]; }
  void display_info(uint32_t r, vec<uint64_t>& se_to_count);
  // TODO: void merge(crecord_sptr_t rhs);

private:
  se_t nnodes = 0;
  se_t nsubsets = 0;
  uint32_t num_partials = 1;
  // Remove tree from this completely or make it a star tree.
  tree_sptr_t tree = nullptr;
  std::vector<pse_t> se_to_pse = {};
  std::vector<double> se_to_rho = {};
};

#endif
