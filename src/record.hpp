#ifndef _RECORD_H
#define _RECORD_H

#define MMHSEED1 1
#define MMHSEED0 0

#include "MurmurHash3.hpp"
#include "common.hpp"
#include "tree.hpp"

class Subset : public std::enable_shared_from_this<Subset>
{
  friend class Record;
  friend class CRecord;

private:
  sh_t shash = 0;
  sh_t chash = 0;
  tuint_t card = 0;

public:
  Subset(sh_t shash1, sh_t shash2, record_sptr_t record);
  Subset(sh_t shash, sh_t chash, tuint_t card)
    : shash(shash)
    , chash(chash)
    , card(card)
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
  bool tree_with_collision();
  void rehash_tree();
  void union_record(record_sptr_t source);
  void add_subset(subset_sptr_t new_subset);
  record_sptr_t getptr() { return shared_from_this(); }
  tree_sptr_t get_tree() { return tree; }
  // bool check_collision(subset_sptr_t x);
  // bool check_conflict(record_sptr_t x);
  void make_compact();
  se_t map_compact(sh_t shash) { return sh_to_se[shash]; }

private:
  std::unordered_map<sh_t, subset_sptr_t> sh_to_subset = {};
  std::unordered_map<sh_t, node_sptr_t> sh_to_node = {};
  std::unordered_map<sh_t, se_t> sh_to_se = {};
  node_sptr_t subtree_root = nullptr;
  tree_sptr_t tree = nullptr;
};

class CRecord : public std::enable_shared_from_this<CRecord>
{
public:
  CRecord(tree_sptr_t tree);
  CRecord(record_sptr_t record);
  vec<se_t> decode_se(se_t se);
  void load(std::filesystem::path library_dir, std::string suffix);
  void save(std::filesystem::path library_dir, std::string suffix);
  crecord_sptr_t getptr() { return shared_from_this(); }
  bool check_compatible(crecord_sptr_t crecord);
  void merge(crecord_sptr_t crecord);

private:
  tree_sptr_t tree;
  se_t nsubsets = 0;
  std::unordered_map<se_t, std::pair<se_t, se_t>> se_to_pse = {};
  std::unordered_map<se_t, node_sptr_t> se_to_node = {};
};

#endif