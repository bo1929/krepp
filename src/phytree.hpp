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
  void parse_lineages(std::ifstream& tree_stream);
  void parse(std::filesystem::path nwk_path);
  void save(std::ofstream& tree_stream);
  void load(std::ifstream& tree_stream);
  void generate_tree(vec<std::string>& names_v);
  void map_to_qtree(tree_sptr_t qtree);
  static double compute_distance(node_sptr_t a, node_sptr_t b);
  static node_sptr_t compute_lca(node_sptr_t x, node_sptr_t y);
  void stream_nwk_basic(std::stringstream& nwk_strstream, node_sptr_t nd);
  void stream_nwk_jplace(std::stringstream& nwk_strstream, node_sptr_t nd);
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
  double tblen = 0;
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
  Node(tree_sptr_t tree, std::string name, node_sptr_t parent, bool is_leaf = false)
    : tree(tree)
    , name(name)
    , parent(parent)
    , is_leaf(is_leaf)
  {
    if (is_leaf) {
      card = 1;
    }
    ix_child = parent ? parent->get_nchildren() : -1;
    sh = hash_name(name);
    while (!sh) {
      sh = rehash(reinterpret_cast<uint64_t>(&sh));
    }
    blen = std::numeric_limits<double>::quiet_NaN();
    bdepth = std::numeric_limits<double>::quiet_NaN();
    ldepth = parent ? parent->get_ldepth() + 1 : 0;
    tblen = std::numeric_limits<double>::quiet_NaN();
  }
  void print_info();
  void parse(vec<std::string>& n_vec);
  void generate_tree(vec_str_iter name_first, vec_str_iter name_last);
  void set_sh(sh_t source) { sh = source; }
  void set_se(se_t source) { se = source; }
  void set_rank(std::string r)
  {
    rank = r;
    is_taxon = true;
  }
  void set_parent(node_sptr_t source)
  {
    if (!source) {
      return;
    }
    parent = source;
    ix_child = parent->get_nchildren();
    parent->add_children(getptr());
    ldepth = parent ? parent->get_ldepth() + 1 : 0;
  }
  node_sptr_t getptr() { return shared_from_this(); }
  node_sptr_t* get_children() { return children.data(); }
  void add_children(node_sptr_t child)
  {
    nchildren++;
    tblen += child->get_blen() + child->get_tblen();
    // children.push_back(child);
    card = card + child->get_card();
    sh = sh + child->get_sh();
    is_leaf = false;
  }
  tuint_t get_nchildren() { return nchildren; }
  node_sptr_t get_parent() { return parent; }
  tree_sptr_t get_tree() { return tree; }
  double get_bdepth() { return bdepth; }
  double get_blen() { return blen; }
  double get_tblen() { return tblen; }
  double get_midpoint_pendant()
  {
    if (!std::isnan(blen)) {
      return blen / 2.0;
    } else {
      return 0;
    }
  }
  double get_ldepth() { return ldepth; }
  std::string const get_name(bool return_na = false)
  {
    if (!name.empty()) {
      return name;
    } else {
      if (return_na) {
        return "NA";
      } else {
        return std::to_string(se - 1);
      }
    }
  }
  void stream_nwk_entry(strstream& nwk_strstream)
  {
    if (std::isnan(blen)) {
      nwk_strstream << name;
    } else {
      nwk_strstream << name << ":" << blen;
    }
  }
  tuint_t get_card() { return card; }
  sh_t get_sh() { return sh; }
  se_t get_se() { return se; }
  se_t get_en() { return se - 1; }
  bool check_leaf() { return is_leaf; }
  bool check_taxon() { return is_taxon; }
  bool is_labeled() { return !name.empty(); }
  sh_t sum_children_sh()
  {
    sh_t sh = 0;
    std::for_each(children.begin(), children.end(), [&sh](node_sptr_t nd) { sh += nd->sh; });
    return sh;
  }

private:
  vec<node_sptr_t> children;
  std::string name = "";
  std::string rank = "";
  node_sptr_t parent = nullptr;
  tree_sptr_t tree = nullptr;
  double blen = 0;
  double bdepth = 0;
  uint32_t ldepth = 0;
  double tblen = 0;
  bool is_leaf = true;
  bool is_taxon = false;
  tuint_t nchildren = 0;
  tuint_t ix_child = 0;
  tuint_t card = 0;
  sh_t sh = 0;
  se_t se = 0;
};

#endif
