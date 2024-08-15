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
  vec<std::string> n_vec;
  split_nwk(n_vec);
  root = std::make_shared<Node>(getptr());
  root->parse(n_vec);
  subtree_root = root;
}

void Tree::split_nwk(vec<std::string>& n_vec)
{
  int i = 0;
  int at = 0;
  char buf[nwk_str.length()];
  std::memset(buf, 0, nwk_str.length());
  for (; i < nwk_str.length(); i++) {
    if (nwk_str[i] == '\t' || nwk_str[i] == ' ' || nwk_str[i] == ';')
      continue;
    if (nwk_str[i] == '(' || nwk_str[i] == ')' || nwk_str[i] == ':' || nwk_str[i] == ',') {
      if (strlen(buf) > 0) {
        n_vec.push_back(buf);
        std::memset(buf, 0, nwk_str.length());
        at = 0;
      }
      n_vec.push_back(std::string() + nwk_str[i]);
      continue;
    } else {
      buf[at++] = nwk_str[i];
    }
  }
  if (strlen(buf) > 0) {
    n_vec.push_back(buf);
    std::memset(buf, 0, nwk_str.length());
    at = 0;
  }
}

void Node::parse(vec<std::string>& n_vec)
{
  if (tree->atter >= n_vec.size())
    return;
  if (n_vec[tree->atter] != "(") {
    if (n_vec[tree->atter] != ":") {
      name = n_vec[tree->atter];
      tree->atter++;
    }
    if (n_vec[tree->atter] == ":") {
      branch_len = std::atof(n_vec[tree->atter + 1].c_str());
      tree->total_len_branch += branch_len;
      ldepth = parent->ldepth + 1;
      bdepth = parent->bdepth + branch_len;
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
    return;
  }
  if (n_vec[tree->atter] == "(") {
    while (1) {
      tree->atter++;
      node_sptr_t child = std::make_shared<Node>(tree);
      child->set_parent(getptr());
      child->ix_child = nchildren;
      child->parse(n_vec);
      card += child->card;
      sh += child->sh;
      children.push_back(std::move(child));
      nchildren++;
      if (n_vec[tree->atter] == ",")
        continue;
      else
        break;
    }
    is_leaf = false;
    tree->nnodes++;
    se = tree->nnodes;
  }
  if (n_vec[tree->atter] == ")") {
    tree->atter++;
    if (n_vec[tree->atter] == ")")
      return;
  }
  if (n_vec[tree->atter] != ":") {
    name = n_vec[tree->atter];
    tree->atter++;
  }
  if (n_vec[tree->atter] == ":") {
    branch_len = std::atof(n_vec[tree->atter + 1].c_str());
    tree->total_len_branch += branch_len;
    ldepth = parent ? parent->ldepth + 1 : 0;
    bdepth = parent ? parent->bdepth + branch_len : 0;
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
  std::cout << name << "\t" << sh << "\t" << branch_len << "\t" << pname << std::endl;
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
  std::cout << "\tTotal branch length : " << total_len_branch << std::endl;
}

node_sptr_t Tree::compute_lca(node_sptr_t a, node_sptr_t b)
{
  if (!a || !b) // LCA(x,0) = LCA(0,x) = x
    return a ? a : b;
  while (a != b) {
    if (a->ldepth > b->ldepth)
      b = b->parent;
    else
      a = a->parent;
  }
  return a;
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
  vec<std::string> n_vec;
  split_nwk(n_vec);
  root = std::make_shared<Node>(getptr());
  atter = 0, nnodes = 0, total_len_branch = 0;
  root->parse(n_vec);
  subtree_root = root;
}
