#ifndef _PHYTREE_H
#define _PHYTREE_H

#include "common.hpp"
#include "record.hpp"

typedef std::vector<std::string>::const_iterator vec_str_iter;

class Tree : public std::enable_shared_from_this<Tree>
{
  friend class Node;

public:
  Tree() {}
  void print_info();
  node_sptr_t next_post_order();
  node_sptr_t next_post_order(node_sptr_t nd_curr);
  void compute_bdepth();
  bool check_compatible(tree_sptr_t tree);
  void split_nwk(vec<std::string>& n_vec);
  void parse(std::filesystem::path nwk_path);
  void save(std::ofstream& tree_stream);
  void load(std::ifstream& tree_stream);
  void generate_tree(vec<std::string>& names_v);
  void map_to_qtree(tree_sptr_t qtree);
  static double compute_distance(node_sptr_t a, node_sptr_t b);
  static node_sptr_t compute_lca(node_sptr_t x, node_sptr_t y);
  void stream_newick_str(std::stringstream& nwk_strstream, node_sptr_t nd);
  void set_subtree(node_sptr_t source) { subtree_root = source; }
  node_sptr_t get_node(se_t se) const { return se_to_node[se]; }
  bool check_node(se_t se) const { return se <= nnodes; }
  tree_sptr_t getptr() { return shared_from_this(); }
  node_sptr_t get_root() { return root; }
  node_sptr_t get_subtree() { return subtree_root; }
  tuint_t get_nnodes() { return nnodes; }
  void reset_traversal()
  {
    curr = nullptr;
    subtree_root = root;
  }

private:
  std::string nwk_str;
  tuint_t atter = 0;
  tuint_t nnodes = 0;
  double total_blen = 0;
  node_sptr_t root = nullptr;
  node_sptr_t curr = nullptr;
  node_sptr_t subtree_root = nullptr;
  std::vector<node_sptr_t> se_to_node = {nullptr};
};

class Node : public std::enable_shared_from_this<Node>
{
  friend class Tree;
  friend class Subset;

public:
  Node(tree_sptr_t tree)
    : tree(tree)
  {
  }
  void print_info();
  void parse(vec<std::string>& n_vec);
  void generate_tree(vec_str_iter name_first, vec_str_iter name_last);
  void set_sh(sh_t source) { sh = source; }
  void set_parent(node_sptr_t source) { parent = source; }
  node_sptr_t getptr() { return shared_from_this(); }
  node_sptr_t* get_children() { return children.data(); }
  tuint_t get_nchildren() { return nchildren; }
  node_sptr_t get_parent() { return parent; }
  tree_sptr_t get_tree() { return tree; }
  double get_bdepth() { return bdepth; }
  double get_blen() { return blen; }
  double get_ldepth() { return ldepth; }
  std::string get_name() { return name; }
  tuint_t get_card() { return card; }
  sh_t get_sh() { return sh; }
  se_t get_se() { return se; }
  bool check_leaf() { return is_leaf; }
  sh_t sum_children_sh()
  {
    sh_t sh = 0;
    std::for_each(children.begin(), children.end(), [&sh](node_sptr_t nd) { sh += nd->sh; });
    return sh;
  }
  double total_blen = 0;

private:
  vec<node_sptr_t> children;
  std::string name = "";
  node_sptr_t parent = nullptr;
  tree_sptr_t tree = nullptr;
  double blen = 0;
  double bdepth = 0;
  uint32_t ldepth = 0;
  bool is_leaf = true;
  tuint_t nchildren = 0;
  tuint_t ix_child = 0;
  tuint_t card = 0;
  sh_t sh = 0;
  se_t se = 0;
};

#endif
