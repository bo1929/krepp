#ifndef _TREE_H
#define _TREE_H

#include "common.hpp"
#include "record.hpp"

class Tree : public std::enable_shared_from_this<Tree>
{
  friend class Node;

public:
  Tree() {}
  void set_subtree(node_sptr_t source) { subtree_root = source; }
  tree_sptr_t getptr() { return shared_from_this(); }
  node_sptr_t get_root() { return root; }
  void reset_traversal()
  {
    curr = nullptr;
    subtree_root = root;
  }
  void print_info();
  node_sptr_t next_post_order();
  bool check_compatible(tree_sptr_t tree);
  static node_sptr_t compute_lca(node_sptr_t x, node_sptr_t y);
  void parse(std::filesystem::path nwk_path);
  void split_nwk(vec<std::string>& n_vec);
  void save(std::filesystem::path library_dir, std::string suffix);
  void load(std::filesystem::path library_dir, std::string suffix);

private:
  tuint_t atter = 0;
  tuint_t nnodes = 0;
  double total_len_branch = 0;
  node_sptr_t subtree_root = nullptr;
  node_sptr_t root = nullptr;
  node_sptr_t curr = nullptr;
  std::string nwk_str;
};

class Node : public std::enable_shared_from_this<Node>
{
  friend class Tree;
  friend class Subset;

public:
  Node(tree_sptr_t tree)
    : tree(tree)
  {}
  void print_info();
  void parse(vec<std::string>& n_vec);
  void set_shash(sh_t source) { shash = source; }
  void set_parent(node_sptr_t source) { parent = source; }
  node_sptr_t getptr() { return shared_from_this(); }
  vec<node_sptr_t> get_children() { return children; }
  tuint_t get_nchildren() { return nchildren; }
  tree_sptr_t get_tree() { return tree; }
  std::string get_name() { return name; }
  node_sptr_t get_parent() { return parent; }
  bool check_leaf() { return is_leaf; }
  tuint_t get_card() { return card; }
  sh_t get_shash() { return shash; }
  sh_t sum_children_shash()
  {
    sh_t sh = 0;
    std::for_each(children.begin(), children.end(), [&sh](node_sptr_t nd) { sh += nd->shash; });
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
  tuint_t nchildren = 0;
  tuint_t ix_child = 0;
  tuint_t card = 0;
  sh_t shash = 0;
};

#endif
