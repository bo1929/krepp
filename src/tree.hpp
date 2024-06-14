#ifndef _TREE_H
#define _TREE_H

#include "common.hpp"
#include "subsets.hpp"

class Tree : public std::enable_shared_from_this<Tree> {
  friend class Node;

public:
  Tree(std::string nwk_filepath) : nwk_filepath(nwk_filepath) {}
  void parse();
  void print_info();
  void split_nwk(vec<std::string> &n_vec);
  node_sptr_t next_post_order();
  static node_sptr_t compute_lca(node_sptr_t x, node_sptr_t y);
  node_sptr_t get_root() { return root; }
  tree_sptr_t getptr() { return shared_from_this(); }
  void set_subtree(node_sptr_t nd) { subtree_root = nd; }
  void reset_traversal() {
    curr = nullptr;
    subtree_root = root;
  }

private:
  std::string nwk_filepath;
  tuint atter = 0;
  tuint nnodes = 0;
  double total_len_branch = 0;
  node_sptr_t root = nullptr;
  node_sptr_t subtree_root = nullptr;
  node_sptr_t curr = nullptr;
};

class Node : public std::enable_shared_from_this<Node> {
  friend class Tree;
  friend class Subset;

public:
  Node(tree_sptr_t tree) : tree(tree) {}
  void parse(vec<std::string> &n_vec);
  void print_info();
  inline sh_t get_shash() { return shash; }
  inline bool check_leaf() { return is_leaf; }
  inline void set_shash(sh_t sh) { shash = sh; }
  inline tree_sptr_t get_tree() { return tree; }
  inline void set_parent(node_sptr_t nd) { parent = nd; }
  inline node_sptr_t getptr() { return shared_from_this(); }
  sh_t sum_children_shash() {
    sh_t sh = 0;
    std::for_each(children.begin(), children.end(),
                  [&](node_sptr_t nd) { sh += nd->shash; });
    return sh;
  }

private:
  std::string name = "";
  node_sptr_t parent = nullptr;
  tree_sptr_t tree = nullptr;
  vec<node_sptr_t> children;
  bool is_leaf = true;
  double len_branch = 0;
  double depth_level = 0;
  double depth_branch = 0;
  tuint nchildren = 0;
  tuint ix_child = 0;
  tuint card = 1;
  sh_t shash = 0;
};

#endif
