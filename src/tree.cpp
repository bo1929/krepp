#include "tree.hpp"
#include "common.hpp"

bool Tree::check_compatible(tree_sptr_t tree)
{
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
  tree->nnodes++;
  if (tree->atter >= n_vec.size())
    return;
  if (n_vec[tree->atter] != "(") {
    if (n_vec[tree->atter] != ":") {
      name = n_vec[tree->atter];
      tree->atter++;
    }
    if (n_vec[tree->atter] == ":") {
      len_branch = std::atof(n_vec[tree->atter + 1].c_str());
      tree->total_len_branch += len_branch;
      depth_level = parent->depth_level + 1;
      depth_branch = parent->depth_branch + len_branch;
      tree->atter += 2;
    }
    is_leaf = true;
    card = 1;
    shash = Subset::get_singleton_shash(name);
    while (!shash)
      shash = Subset::rehash(reinterpret_cast<uint64_t>(&shash));
    return;
  }
  if (n_vec[tree->atter] == "(") {
    while (1) {
      tree->atter++;
      node_sptr_t child = std::make_shared<Node>(tree);
      child->set_parent(getptr());
      child->ix_child = nchildren;
      child->parse(n_vec);
      shash += child->shash;
      card += child->card;
      children.push_back(std::move(child));
      nchildren++;
      if (n_vec[tree->atter] == ",")
        continue;
      else
        break;
    }
    is_leaf = false;
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
    len_branch = std::atof(n_vec[tree->atter + 1].c_str());
    tree->total_len_branch += len_branch;
    depth_level = parent ? parent->depth_level + 1 : 0;
    depth_branch = parent ? parent->depth_branch + len_branch : 0;
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
  std::cout << name << "\t" << shash << "\t" << len_branch << "\t" << pname << std::endl;
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
    if (a->depth_level > b->depth_level)
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
