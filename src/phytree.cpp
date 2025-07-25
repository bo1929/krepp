#include "phytree.hpp"

bool Tree::check_compatible(tree_sptr_t tree)
{
  if (!tree) return true;
  bool is_compatible = true;
  node_sptr_t ndx, ndy;
  tree->reset_traversal();
  reset_traversal();
  while (true) {
    ndx = next_post_order();
    ndy = tree->next_post_order();
    if (!ndx && !ndy) {
      is_compatible = true;
      break;
    } else if (!ndx || !ndy) {
      is_compatible = false;
      break;
    } else if (ndx->name != ndy->name) {
      is_compatible = false;
      break;
    } else {
      continue;
    }
  }
  tree->reset_traversal();
  reset_traversal();
  return is_compatible;
}

void Tree::generate_tree(vec<std::string>& names_v)
{
  root = std::make_shared<Node>(getptr());
  root->generate_tree(names_v.begin(), names_v.end());
  se_to_node.push_back(root);
  subtree_root = root;
  compute_bdepth();
}

void Tree::stream_newick_str(strstream& nwk_strstream, node_sptr_t nd)
{
  if (!nd->check_leaf()) {
    nwk_strstream << "(";
    for (uint32_t nix = 0; nix < nd->get_nchildren(); ++nix) {
      stream_newick_str(nwk_strstream, *std::next(nd->get_children(), nix));
      if (nix < (nd->get_nchildren() - 1)) {
        nwk_strstream << ",";
      }
    }
    nwk_strstream << ")";
  }
  nwk_strstream << nd->get_name() << ":" << nd->get_blen() << "{" << nd->get_se() - 1 << "}";
  if (nd == get_root()) {
    nwk_strstream << ";";
  }
}

void Tree::split_nwk(vec<std::string>& nd_v)
{
  int i = 0;
  std::string buf = "";
  bool is_label = false, quote = false;
  buf.reserve(nd_v.size());
  for (; i < nwk_str.length(); i++) {
    quote = (nwk_str[i] == '\'' || nwk_str[i] == '"');
    is_label = (is_label != quote);
    if (is_label && !quote) {
      buf += nwk_str[i];
    } else if (quote || nwk_str[i] == '\n' || nwk_str[i] == ';') {
      continue;
    } else if (nwk_str[i] == '(' || nwk_str[i] == ')' || nwk_str[i] == ':' || nwk_str[i] == ',') {
      if (nwk_str[i] != '(' && nwk_str[i - 1] != '(') {
        if (buf.empty()) buf = "''";
        nd_v.push_back(buf);
        buf = "";
      }
      nd_v.push_back(std::string() + nwk_str[i]);
    } else {
      buf += nwk_str[i];
    }
  }
  if (buf.length() > 0) {
    nd_v.push_back(buf);
  }
}

void Node::parse(vec<std::string>& nd_v)
{
  ldepth = parent ? parent->ldepth + 1 : 0;
  if (tree->atter >= nd_v.size()) return;
  if (nd_v[tree->atter] != "(") {
    if (nd_v[tree->atter] != ":") {
      name = nd_v[tree->atter];
      tree->atter++;
    }
    if (nd_v[tree->atter] == ":") {
      blen = std::atof(nd_v[tree->atter + 1].c_str());
      tree->total_blen += blen;
      tree->atter += 2;
    }
    is_leaf = true;
    card = 1;
    sh = Subset::get_singleton_sh(name);
    while (!sh) {
      sh = Subset::rehash(reinterpret_cast<uint64_t>(&sh));
    }
    tree->nnodes++;
    se = tree->nnodes;
    tree->se_to_node.push_back(getptr());
    return;
  }
  if (nd_v[tree->atter] == "(") {
    while (1) {
      tree->atter++;
      children.emplace_back(std::make_shared<Node>(tree));
      (children.back())->set_parent(getptr());
      (children.back())->ix_child = nchildren;
      (children.back())->parse(nd_v);
      card += (children.back())->card;
      sh += (children.back())->sh;
      total_blen += (children.back())->blen;
      total_blen += (children.back())->total_blen;
      nchildren++;
      if (nd_v[tree->atter] == ",")
        continue;
      else
        break;
    }
    is_leaf = false;
    tree->nnodes++;
    se = tree->nnodes;
    tree->se_to_node.push_back(getptr());
  }
  if (nd_v[tree->atter] == ")") {
    tree->atter++;
    if (nd_v[tree->atter] == ")") return;
  }
  if (nd_v[tree->atter] != ":") {
    name = nd_v[tree->atter];
    tree->atter++;
  }
  if (nd_v[tree->atter] == ":") {
    blen = std::atof(nd_v[tree->atter + 1].c_str());
    tree->total_blen += blen;
    tree->atter += 2;
  }
  // if (!name.empty() && name[name.length() - 1] == '\n') {
  //   name.erase(name.length() - 1);
  // }
  return;
}

void Node::generate_tree(vec_str_iter name_first, vec_str_iter name_last)
{
  std::size_t diff_size = (name_last - name_first);
  if (diff_size == 1) {
    name = *(name_first);
    blen = 1.0;
    tree->total_blen += blen;
    is_leaf = true;
    card = 1;
    sh = Subset::get_singleton_sh(name);
    while (!sh) {
      sh = Subset::rehash(reinterpret_cast<uint64_t>(&sh));
    }
    tree->nnodes++;
    se = tree->nnodes;
    tree->se_to_node.push_back(getptr());
  } else {
    vec_str_iter name_half = std::next(name_first, diff_size / 2);
    for (uint32_t pix = 0; pix < 2; ++pix) {
      children.emplace_back(std::make_shared<Node>(tree));
      (children.back())->set_parent(getptr());
      (children.back())->ix_child = nchildren;
      if (pix) {
        (children.back())->generate_tree(name_first, name_half);
      } else {
        (children.back())->generate_tree(name_half, name_last);
      }
      card += (children.back())->card;
      sh += (children.back())->sh;
      total_blen += (children.back())->blen;
      total_blen += (children.back())->total_blen;
      nchildren++;
    }
    blen = 1.0;
    is_leaf = false;
    tree->nnodes++;
    se = tree->nnodes;
    name = "N" + std::to_string(se);
    tree->se_to_node.push_back(getptr());
    tree->total_blen += blen;
  }
}

void Node::print_info()
{
  std::string pname = parent ? parent->name : "";
  std::cout << name << "\t" << sh << "\t" << blen << "\t" << pname << std::endl;
}

node_sptr_t Tree::next_post_order(node_sptr_t nd_curr)
{
  if (nd_curr == subtree_root) return nullptr;
  if (!nd_curr) {
    nd_curr = subtree_root;
    while (!nd_curr->is_leaf)
      nd_curr = nd_curr->children.front();
    return nd_curr;
  } else {
    if (nd_curr->is_leaf) {
      if (nd_curr == nd_curr->parent->children.back()) {
        nd_curr = nd_curr->parent;
        return nd_curr;
      } else {
        nd_curr = nd_curr->parent->children[nd_curr->ix_child + 1];
        while (!nd_curr->is_leaf)
          nd_curr = nd_curr->children.front();
        return nd_curr;
      }
    } else {
      if (nd_curr == nd_curr->parent->children.back()) {
        nd_curr = nd_curr->parent;
        return nd_curr;
      } else {
        nd_curr = nd_curr->parent->children[nd_curr->ix_child + 1];
        while (!nd_curr->is_leaf)
          nd_curr = nd_curr->children.front();
        return nd_curr;
      }
    }
  }
}

node_sptr_t Tree::next_post_order()
{
  curr = next_post_order(curr);
  return curr;
}

void Tree::print_info()
{
  std::cout << "The tree rooted at : " << root->name << std::endl;
  std::cout << "\tNumber of nodes : " << nnodes << std::endl;
  std::cout << "\tTotal branch length : " << total_blen << std::endl;
}

node_sptr_t Tree::compute_lca(node_sptr_t a, node_sptr_t b)
{
  if (!a || !b) // LCA(x,0) = LCA(0,x) = x
    return a ? a : b;
  while (a->sh != b->sh) {
    if (a->ldepth < b->ldepth)
      b = b->parent;
    else
      a = a->parent;
  }
  return a;
}

double Tree::compute_distance(node_sptr_t a, node_sptr_t b)
{
  double distance = 0;
  if (!a || !b) return std::numeric_limits<double>::max();
  while (a->sh != b->sh) {
    if (a->ldepth < b->ldepth) {
      distance += b->get_blen();
      b = b->parent;
    } else {
      distance += a->get_blen();
      a = a->parent;
    }
  }
  return distance;
}

void Tree::save(std::ofstream& tree_stream)
{
  std::ostream_iterator<char> output_iterator(tree_stream);
  std::copy(nwk_str.begin(), nwk_str.end(), output_iterator);
}

void Tree::load(std::ifstream& tree_stream)
{
  nwk_str =
    std::string((std::istreambuf_iterator<char>(tree_stream)), std::istreambuf_iterator<char>());
  vec<std::string> nd_v;
  split_nwk(nd_v);
  root = std::make_shared<Node>(getptr());
  atter = 0, nnodes = 0, total_blen = 0;
  root->parse(nd_v);
  subtree_root = root;
  compute_bdepth();
}

void Tree::compute_bdepth()
{
  std::queue<node_sptr_t> nd_q;
  node_sptr_t nd_curr;
  nd_q.push(root);
  while (!nd_q.empty()) {
    nd_curr = nd_q.front();
    nd_q.pop();
    for (tuint_t i = 0; i < nd_curr->get_nchildren(); ++i) {
      nd_q.push((*std::next(nd_curr->get_children(), i)));
      (nd_q.front())->bdepth = nd_curr->bdepth + (nd_q.front())->blen;
    }
  }
}

void Tree::map_to_qtree(tree_sptr_t qtree)
{
  parallel_flat_phmap<std::string, se_t> name_to_se;
  reset_traversal();
  while (curr = next_post_order()) {
    if (curr->check_leaf()) {
      name_to_se[curr->get_name()] = curr->get_se();
      se_to_node[curr->get_se()] = nullptr;
    }
  }
  root = qtree->get_root();
  subtree_root = qtree->get_subtree();
  reset_traversal();
  while (curr = next_post_order()) {
    if (curr->check_leaf()) {
      if (name_to_se.contains(curr->get_name())) {
        se_to_node[name_to_se[curr->get_name()]] = curr;
      } else {
        error_exit("Given placement tree contains a reference that does not appear in the index.");
      }
    }
  }
  reset_traversal();
}
