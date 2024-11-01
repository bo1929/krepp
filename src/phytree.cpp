#include "phytree.hpp"

bool Tree::check_compatible(tree_sptr_t tree)
{
  if (!tree)
    return true;
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

void Tree::parse(std::filesystem::path nwk_path)
{
  std::ifstream tree_stream(nwk_path);
  nwk_str =
    std::string((std::istreambuf_iterator<char>(tree_stream)), std::istreambuf_iterator<char>());
  tree_stream.close();
  vec<std::string> nd_v;
  split_nwk(nd_v);
  root = std::make_shared<Node>(getptr());
  root->parse(nd_v);
  subtree_root = root;
  compute_bdepth();
}

void Tree::split_nwk(vec<std::string>& nd_v)
{
  int i = 0;
  std::string buf = "";
  buf.reserve(nd_v.size());
  for (; i < nwk_str.length(); i++) {
    if (nwk_str[i] == '\t' || nwk_str[i] == ' ' || nwk_str[i] == ';')
      continue;
    if (nwk_str[i] == '(' || nwk_str[i] == ')' || nwk_str[i] == ':' || nwk_str[i] == ',') {
      if (buf.length() > 0) {
        nd_v.push_back(buf);
        buf = "";
      }
      nd_v.push_back(std::string() + nwk_str[i]);
      continue;
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
  if (tree->atter >= nd_v.size())
    return;
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
    if (nd_v[tree->atter] == ")")
      return;
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
  if (!name.empty() && name[name.length() - 1] == '\n') {
    name.erase(name.length() - 1);
  }
  return;
}

void Node::print_info()
{
  std::string pname = parent ? parent->name : "";
  std::cout << name << "\t" << sh << "\t" << blen << "\t" << pname << std::endl;
}

node_sptr_t Tree::next_post_order()
{
  if (curr == subtree_root)
    return nullptr;
  if (!curr) {
    curr = subtree_root;
    while (!curr->is_leaf)
      curr = curr->children.front();
    return curr;
  } else {
    if (curr->is_leaf) {
      if (curr == curr->parent->children.back()) {
        curr = curr->parent;
        return curr;
      } else {
        curr = curr->parent->children[curr->ix_child + 1];
        while (!curr->is_leaf)
          curr = curr->children.front();
        return curr;
      }
    } else {
      if (curr == curr->parent->children.back()) {
        curr = curr->parent;
        return curr;
      } else {
        curr = curr->parent->children[curr->ix_child + 1];
        while (!curr->is_leaf)
          curr = curr->children.front();
        return curr;
      }
    }
  }
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
  if (!a || !b)
    return std::numeric_limits<double>::max();
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

void Tree::save(std::filesystem::path library_dir, std::string suffix)
{
  std::ofstream tree_stream(library_dir / ("tree" + suffix));
  std::ostream_iterator<char> output_iterator(tree_stream);
  std::copy(nwk_str.begin(), nwk_str.end(), output_iterator);
  if (!tree_stream.good()) {
    std::cerr << "Writing the parsed reference tree has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  tree_stream.close();
}

void Tree::load(std::filesystem::path library_dir, std::string suffix)
{
  std::ifstream tree_stream(library_dir / ("tree" + suffix));
  nwk_str =
    std::string((std::istreambuf_iterator<char>(tree_stream)), std::istreambuf_iterator<char>());
  tree_stream.close();
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
