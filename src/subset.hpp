#ifndef _SUBSETS_H
#define _SUBSETS_H

#define MMHSEED1 1
#define MMHSEED0 0

#include "MurmurHash3.hpp"
#include "common.hpp"
#include "tree.hpp"

class Subset : public std::enable_shared_from_this<Subset>
{
  friend class Record;

private:
  sh_t shash = 0;
  sh_t chash = 0;
  tuint card = 0;

public:
  Subset(sh_t shash1, sh_t shash2, record_sptr_t record);
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
  friend class Subset;

public:
  Record(node_sptr_t nd);
  Record(record_sptr_t source1, record_sptr_t source2);
  bool tree_with_collision();
  void rehash_tree();
  void union_record(record_sptr_t source);
  void add_subset(subset_sptr_t new_subset);
  record_sptr_t getptr() { return shared_from_this(); }
  // bool check_collision(subset_sptr_t x);
  // bool check_conflict(record_sptr_t x);

private:
  std::unordered_map<sh_t, subset_sptr_t> sh_to_subset = {};
  std::unordered_map<sh_t, node_sptr_t> sh_to_node = {};
  node_sptr_t subtree_root = nullptr;
  tree_sptr_t tree = nullptr;
  omp_lock_t main_lock;
};

#endif
