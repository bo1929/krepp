#ifndef _RECORD_H
#define _RECORD_H

#define MMHSEED1 1
#define MMHSEED0 0

#include "MurmurHash3.hpp"
#include "common.hpp"
#include "phytree.hpp"

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
  {}
  subset_sptr_t getptr() { return shared_from_this(); }
  static inline sh_t get_singleton_shash(std::string& name)
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
  Record(node_sptr_t nd);
  Record(record_sptr_t source1, record_sptr_t source2);
  se_t map_compact(sh_t sh) { return shash_to_senc[sh]; }
  record_sptr_t getptr() { return shared_from_this(); }
  tree_sptr_t get_tree() { return tree; }
  sh_t add_subset(sh_t sh1, sh_t sh2);
  void union_record(record_sptr_t source);
  bool check_subset_collision(sh_t sh, subset_sptr_t subset1, subset_sptr_t subset2);
  // bool check_record_conflict(record_sptr_t x);
  vec<node_sptr_t> decode_shash(sh_t sh);
  bool check_tree_collision();
  void rehash_tree();
  void make_compact();

private:
  tree_sptr_t tree = nullptr;
  node_sptr_t subtree_root = nullptr;
  parallel_flat_phmap<sh_t, subset_sptr_t> shash_to_subset = {};
  parallel_flat_phmap<sh_t, node_sptr_t> shash_to_node = {};
  parallel_flat_phmap<sh_t, se_t> shash_to_senc = {};
};

class CRecord : public std::enable_shared_from_this<CRecord>
{
public:
  CRecord(tree_sptr_t tree);
  CRecord(record_sptr_t record);
  std::pair<se_t, se_t> get_psenc(se_t se) { return senc_to_psenc[se]; }
  void load(std::filesystem::path library_dir, std::string suffix);
  void save(std::filesystem::path library_dir, std::string suffix);
  bool check_node(se_t se) { return senc_to_node.contains(se); }
  node_sptr_t get_node(se_t se) { return senc_to_node[se]; }
  crecord_sptr_t getptr() { return shared_from_this(); }
  se_t get_nnodes() { return senc_to_node.size(); }
  bool check_compatible(crecord_sptr_t crecord);
  vec<node_sptr_t> decode_senc(se_t se);
  void merge(crecord_sptr_t crecord);
  void print_info();

private:
  tree_sptr_t tree;
  se_t nsubsets = 0;
  parallel_flat_phmap<se_t, std::pair<se_t, se_t>> senc_to_psenc = {};
  parallel_flat_phmap<se_t, node_sptr_t> senc_to_node = {};
};

#endif
