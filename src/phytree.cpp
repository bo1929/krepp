#include "phytree.hpp"

bool is_number(const std::string& s)
{
  std::istringstream iss(s);
  double d;
  return iss >> std::noskipws >> d && iss.eof();
}

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
  // compute_bdepth();
}

void Tree::stream_nwk_jplace(strstream& nwk_strstream, node_sptr_t nd)
{
  if (!nd->check_leaf()) {
    nwk_strstream << "(";
    for (uint32_t nix = 0; nix < nd->get_nchildren(); ++nix) {
      stream_nwk_jplace(nwk_strstream, *std::next(nd->get_children(), nix));
      if (nix < (nd->get_nchildren() - 1)) {
        nwk_strstream << ",";
      }
    }
    nwk_strstream << ")";
  }
  nd->stream_nwk_entry(nwk_strstream);
  nwk_strstream << "{" << nd->get_en() << "}";
  if (nd == get_root()) {
    nwk_strstream << ";";
  }
}

void Tree::stream_nwk_basic(strstream& nwk_strstream, node_sptr_t nd)
{
  if (!nd->check_leaf()) {
    nwk_strstream << "(";
    for (uint32_t nix = 0; nix < nd->get_nchildren(); ++nix) {
      stream_nwk_basic(nwk_strstream, *std::next(nd->get_children(), nix));
      if (nix < (nd->get_nchildren() - 1)) {
        nwk_strstream << ",";
      }
    }
    nwk_strstream << ")";
  }
  nd->stream_nwk_entry(nwk_strstream);
  if (nd == get_root()) {
    nwk_strstream << ";";
  }
}

void Tree::split_nwk(vec<std::string>& el_v)
{ // TODO: This doesn't parse The Rich Newick format. Mention in the docs.
  std::string buf = "";
  bool is_quoted = false, quote = false, quote_p = false, is_comment = false;
  if (nwk_str.empty()) {
    error_exit("Given Newick tree seems to be empty?!?.");
  }
  if (nwk_str.back() == '\n') {
    nwk_str.pop_back();
  }
  buf.reserve(el_v.size());
  if (nwk_str.back() != ';') {
    error_exit("Given Newick tree ends with a character other than ';'.");
  }
  for (uint32_t i = 0; i < nwk_str.length(); i++) {
    if (is_comment) {
      is_comment = is_comment != (nwk_str[i] == ']');
      continue;
    }
    quote = (nwk_str[i] == '\'' || nwk_str[i] == '"');
    if (quote & quote_p) {
      is_quoted = false;
      buf += "'";
      continue;
    }
    quote_p = quote;
    if (quote) {
      is_quoted = (is_quoted != quote);
      continue;
    } else if (is_quoted) {
      is_comment = is_comment != (nwk_str[i] == '[');
      if (!is_comment) {
        buf += nwk_str[i];
      }
    } else if (nwk_str[i] == '(' || nwk_str[i] == ')' || nwk_str[i] == ':' || nwk_str[i] == ',') {
      if (nwk_str[i] != '(' && nwk_str[i - 1] != '(') {
        if (buf.empty()) buf = "";
        el_v.push_back(buf);
        buf = "";
      }
      el_v.push_back(std::string() + nwk_str[i]);
    } else {
      if (nwk_str[i] == '[' || nwk_str[i] == ']') {
        error_exit("Given Newick tree contains an unquoted label or length with '[' or ']'.");
      }
      if (nwk_str[i] == ';' && i != (nwk_str.length() - 1)) {
        if (nwk_str.length() > (i + 1) && nwk_str[i + 1] == '\n') {
          error_exit("Given Newick file may contain multiple trees, encountered unexpected ';'.");
        } else {
          error_exit("Given Newick tree contains an unquoted label or length with ';'.");
        }
      }
      if ((nwk_str[i] == ' ' || nwk_str[i] == '\n') && !buf.empty()) {
        error_exit("Given Newick tree contains an unquoted label or length with ' ' or newline.");
      }
      buf += nwk_str[i];
    }
  }
  if (buf.length() > 0) {
    el_v.push_back(buf);
  }
}

void Node::parse(vec<std::string>& el_v)
{
  ldepth = parent ? parent->ldepth + 1 : 0;
  if (tree->atter >= el_v.size()) return;
  if (el_v[tree->atter] == "(") {
    while (true) {
      tree->atter++;
      children.emplace_back(std::make_shared<Node>(tree));
      (children.back())->parse(el_v);
      (children.back())->set_parent(getptr());
      if (el_v[tree->atter] == ",")
        continue;
      else
        break;
    }
    if (nchildren == 1) {
      error_exit("A node has a single child in the backbone tree! Please suppress unifurcations.");
    }
    tree->nnodes++;
    se = tree->nnodes;
    tree->se_to_node.push_back(getptr());
    if (el_v[tree->atter] == ")") {
      tree->atter++;
      if (el_v[tree->atter] == ")") return;
    }
    name = "";
    blen = std::numeric_limits<double>::quiet_NaN(); // No branch length...
    if (el_v[tree->atter] != ",") {
      if (el_v[tree->atter] != ":") {
        name = el_v[tree->atter];
        tree->atter++;
      }
      if (el_v[tree->atter] == ":") {
        blen = std::atof(el_v[tree->atter + 1].c_str());
        tree->tblen += blen;
        tree->atter += 2;
      }
    }
    // if (!name.empty() && name[name.length() - 1] == '\n') {
    //   name.erase(name.length() - 1);
    // }
  } else {
    name = "";
    blen = std::numeric_limits<double>::quiet_NaN(); // No branch length...
    if (el_v[tree->atter] != ",") {
      if (el_v[tree->atter] != ":") {
        name = el_v[tree->atter];
        tree->atter++;
      }
      if (el_v[tree->atter] == ":") {
        blen = std::atof(el_v[tree->atter + 1].c_str());
        tree->tblen += blen;
        tree->atter += 2;
      }
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
  }
}

void Node::generate_tree(vec_str_iter name_first, vec_str_iter name_last)
{
  std::size_t diff_size = (name_last - name_first);
  if (diff_size == 1) {
    name = *(name_first);
    blen = 1.0;
    tree->tblen += blen;
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
      if (pix) {
        (children.back())->generate_tree(name_first, name_half);
      } else {
        (children.back())->generate_tree(name_half, name_last);
      }
      (children.back())->set_parent(getptr());
    }
    blen = 1.0;
    is_leaf = false;
    tree->nnodes++;
    se = tree->nnodes;
    // name = "N" + std::to_string(se - 1);
    name = "";
    tree->se_to_node.push_back(getptr());
    tree->tblen += blen;
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
  std::cout << "The tree rooted at : " << root->get_name() << std::endl;
  std::cout << "\tNumber of nodes : " << nnodes << std::endl;
  std::cout << "\tTotal branch length : " << tblen << std::endl;
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

void Tree::parse_lineages(std::ifstream& lineage_stream)
{
  root = std::make_shared<Node>(getptr(), "root", nullptr);
  root->set_rank("root");
  atter = 0, nnodes = 0, tblen = 0;
  subtree_root = root;
  parallel_flat_phmap<std::string, node_sptr_t> taxon_to_node = {};

  std::string line;
  while (std::getline(lineage_stream, line)) {
    line = std::regex_replace(line, std::regex("; "), ";");
    std::istringstream iss(line);
    std::string lineage, name;
    if (!(std::getline(iss, name, '\t') && std::getline(iss, lineage, '\t'))) {
      error_exit("Failed to reference to lineage mapping!");
    }

    std::stringstream lss(lineage);
    std::string taxon, rank;
    node_sptr_t parent = nullptr;
    while (std::getline(lss, taxon, ';')) {
      rank = std::regex_replace(taxon, std::regex("__.*"), "");
      taxon = std::regex_replace(taxon, std::regex(".__"), "");
      if (taxon.empty()) continue;
      if (!taxon_to_node.contains(taxon)) {
        taxon_to_node[taxon] = std::make_shared<Node>(getptr(), taxon, parent);
        if (parent) taxon_to_node[taxon]->set_parent(parent);
        taxon_to_node[taxon]->set_rank(rank);
      }
      parent = taxon_to_node[taxon];
    }

    if (!taxon_to_node.contains(name)) {
      taxon_to_node[name] = std::make_shared<Node>(getptr(), name, parent, true);
      taxon_to_node[name]->set_parent(parent);
    } else {
      error_exit("The same reference appears more than once in the lineage file.");
    }
  }

  for (auto& [taxon, nd] : taxon_to_node) {
    if (!nd->get_parent()) nd->set_parent(root);
  }

  reset_traversal();
  while (curr = next_post_order()) {
    nnodes++;
    se_to_node.push_back(curr);
    curr->set_se(nnodes);
  }
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
  nwk_str = std::string((std::istreambuf_iterator<char>(tree_stream)), std::istreambuf_iterator<char>());
  vec<std::string> el_v;
  split_nwk(el_v);
  root = std::make_shared<Node>(getptr());
  atter = 0, nnodes = 0, tblen = 0;
  root->parse(el_v);
  subtree_root = root;
  // compute_bdepth();
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
    if (curr->check_leaf() && curr->is_labeled()) {
      if (name_to_se.contains(curr->get_name())) {
        se_to_node[name_to_se[curr->get_name()]] = curr;
      } else {
        strstream parsed_tree;
        stream_nwk_basic(parsed_tree, root);
        std::cerr << "Parsed tree: " << parsed_tree.rdbuf() << "\n";
        std::cerr << "Unexpected leaf: " << curr->get_name() << "\n";
        error_exit("Given placement tree contains a reference that does not appear in the index.");
      }
    }
  }
  reset_traversal();
}
