#include "record.hpp"

Record::Record(node_sptr_t nd)
{
  subtree_root = nd;
  tree = subtree_root->get_tree();
  shash_to_node[nd->get_shash()] = nd;
  nd->get_tree()->set_subtree(nd);
  node_sptr_t nd_curr;
  sh_t ch;
  while (nd_curr = nd->get_tree()->next_post_order()) {
    shash_to_node[nd_curr->get_shash()] = nd_curr;
    if (nd_curr->check_leaf()) {
      ch = 0;
    } else {
      ch = (*nd_curr->get_children())->get_shash();
    }
    shash_to_subset[nd_curr->get_shash()] =
      std::make_shared<Subset>(nd_curr->get_shash(), ch, nd_curr->get_card());
  }
  nd->get_tree()->reset_traversal();
  if (subtree_root == subtree_root->get_tree()->get_root()) {
    while (check_tree_collision()) {
      rehash_tree();
    }
  }
  shash_to_subset[0] = std::make_shared<Subset>(0, 0, 0);
}

Record::Record(record_sptr_t source1, record_sptr_t source2)
{
  union_record(source1);
  union_record(source2);
  subtree_root = Tree::compute_lca(source1->subtree_root, source2->subtree_root);
  tree = subtree_root->get_tree();
}

bool Record::check_tree_collision()
{
  bool collision_free = true;
  subtree_root->get_tree()->set_subtree(subtree_root);
  node_sptr_t nd_curr;
  sh_t sh_curr;
  while (collision_free && (nd_curr = subtree_root->get_tree()->next_post_order())) {
    sh_curr = nd_curr->get_shash();
    if (!shash_to_node.contains(sh_curr) || (shash_to_node[sh_curr] != nd_curr) || (!sh_curr)) {
      collision_free = false;
    }
  }
  subtree_root->get_tree()->reset_traversal();
  return !collision_free;
}

void Record::rehash_tree()
{
  shash_to_node.clear();
  shash_to_subset.clear();
  subtree_root->get_tree()->set_subtree(subtree_root);
  node_sptr_t nd_curr;
  sh_t ch;
  while (nd_curr = subtree_root->get_tree()->next_post_order()) {
    if (nd_curr->check_leaf()) {
      nd_curr->set_shash(Subset::rehash(nd_curr->get_shash()));
    } else {
      nd_curr->set_shash(nd_curr->sum_children_shash());
    }
    shash_to_node[nd_curr->get_shash()] = nd_curr;
    if (nd_curr->check_leaf()) {
      ch = 0;
    } else {
      ch = (*nd_curr->get_children())->get_shash();
    }
    shash_to_subset[nd_curr->get_shash()] =
      std::make_shared<Subset>(nd_curr->get_shash(), ch, nd_curr->get_card());
  }
  subtree_root->get_tree()->reset_traversal();
}

void Record::union_record(record_sptr_t source)
{
  // TODO: check conflicts and resolve.
  // TODO: check if either of the records is empty.
  shash_to_node.insert(source->shash_to_node.begin(),
                       source->shash_to_node.end()); // TODO: Maybe remove.
  shash_to_subset.insert(source->shash_to_subset.begin(), source->shash_to_subset.end());
  subtree_root = Tree::compute_lca(subtree_root, source->subtree_root);
  tree = subtree_root->get_tree();
}

bool Record::check_subset_collision(sh_t sh, subset_sptr_t subset1, subset_sptr_t subset2)
{
  if (shash_to_subset.contains(sh)) {
    subset_sptr_t subset = shash_to_subset[sh];
    if (subset->ch == 0 || subset->ch == subset1->sh || subset->ch == subset2->sh) {
      return false;
    } else {
      // TODO: Missing cases where, i.e., subsets other than the kept partition.
      // Perphaps this is enough?
      return true;
    }
  } else {
    return false;
  }
}

sh_t Record::add_subset(sh_t sh1, sh_t sh2)
{
  if (!(shash_to_subset.contains(sh1) && shash_to_subset.contains(sh2))) {
    std::cerr << "Cannot make the subset for the partition: (" << sh1 << ", " << sh2 << ")\n";
    std::quick_exit(EXIT_FAILURE);
  }
  subset_sptr_t subset1 = shash_to_subset[sh1];
  subset_sptr_t subset2 = shash_to_subset[sh2];
  sh_t sh = sh1 + sh2;
  sh_t nonce = 0;
  while (((sh + nonce) == 0) || check_subset_collision(sh + nonce, subset1, subset2)) {
    nonce = Subset::rehash(nonce + sh1 * sh2);
  }
  sh += nonce;
  if (!shash_to_subset.contains(sh)) {
    shash_to_subset[sh] =
      std::make_shared<Subset>(sh,
                               subset1->card > subset2->card ? subset1->sh : subset2->sh,
                               subset1->card + subset2->card,
                               nonce);
  }
  return sh;
}

void Record::make_compact()
{ // TODO: Implement pruning, probably at this point.
  se_t limit_senum = std::numeric_limits<se_t>::max();
  se_t curr_senum = 1;
  tree->reset_traversal();
  node_sptr_t nd_curr;
  while (nd_curr = tree->next_post_order()) {
    if (curr_senum < limit_senum) {
      shash_to_senc[nd_curr->get_shash()] = curr_senum;
      curr_senum++;
    } else {
      std::cerr << "The current se_t size is too small to fit all nodes!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  tree->reset_traversal();
  for (auto& [sh, subset] : shash_to_subset) {
    if (curr_senum < limit_senum) {
      curr_senum += static_cast<se_t>(shash_to_senc.try_emplace(sh, curr_senum).second);
    } else {
      std::cerr << "The current se_t size is too small to fit all subsets observed!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  shash_to_senc[0] = 0;
}

CRecord::CRecord(record_sptr_t record)
{
  record->make_compact();
  tree = record->get_tree();
  tree->reset_traversal();
  se_t curr_senum = 1;
  node_sptr_t nd_curr;
  while (nd_curr = tree->next_post_order()) {
    senc_to_node[curr_senum] = nd_curr;
    curr_senum++;
  }
  tree->reset_traversal();
  for (auto& [sh, subset] : record->shash_to_subset) {
    senc_to_psenc.try_emplace(
      record->shash_to_senc[sh],
      std::make_pair(record->shash_to_senc[subset->ch],
                     record->shash_to_senc[subset->sh - subset->ch - subset->nonce]));
  }
  senc_to_psenc[0] = std::make_pair(0, 0);
  nsubsets = senc_to_psenc.size();
}

void CRecord::print_info()
{
  std::cout << "Total number of subsets excluding nodes: " << nsubsets << std::endl;
  std::cout << "Number of nodes: " << senc_to_node.size() << std::endl;
  for (auto [se, nd] : senc_to_node) {
    std::cout << se << ": " << nd->get_name() << "(" << nd->get_card() << ")" << std::endl;
  }
  for (auto [se, pse] : senc_to_psenc) {
    std::cout << se << ": " << pse.first << "+" << pse.second << std::endl;
  }
}

CRecord::CRecord(tree_sptr_t tree)
  : tree(tree)
{
  tree->reset_traversal();
  node_sptr_t nd_curr;
  se_t curr_senum = 1;
  while (nd_curr = tree->next_post_order()) {
    senc_to_node[curr_senum] = nd_curr;
    curr_senum++;
  }
  tree->reset_traversal();
  nsubsets = senc_to_node.size();
}

void CRecord::load(std::filesystem::path library_dir, std::string suffix)
{
  std::filesystem::path crecord_path = library_dir / ("crecord" + suffix);
  std::ifstream crecord_stream(crecord_path, std::ifstream::binary);
  if (!crecord_stream.is_open()) {
    std::cerr << "Failed to open " << crecord_path << std::endl;
    exit(EXIT_FAILURE);
  } else {
    crecord_stream.read(reinterpret_cast<char*>(&nsubsets), sizeof(se_t));
    std::vector<std::pair<se_t, std::pair<se_t, se_t>>> pse_pair_v(nsubsets);
    crecord_stream.read(reinterpret_cast<char*>(pse_pair_v.data()),
                        sizeof(std::pair<se_t, std::pair<se_t, se_t>>) * nsubsets);
    senc_to_psenc =
      parallel_flat_phmap<se_t, std::pair<se_t, se_t>>(pse_pair_v.begin(), pse_pair_v.end());
  }
  if (!crecord_stream.good()) {
    std::cerr << "Reading subset enumerations has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  crecord_stream.close();
}

void CRecord::save(std::filesystem::path library_dir, std::string suffix)
{
  std::ofstream crecord_stream(library_dir / ("crecord" + suffix), std::ofstream::binary);
  crecord_stream.write(reinterpret_cast<char*>(&nsubsets), sizeof(se_t));
  std::vector<std::pair<se_t, std::pair<se_t, se_t>>> pse_pair_v(senc_to_psenc.begin(),
                                                                 senc_to_psenc.end());
  crecord_stream.write(reinterpret_cast<char*>(pse_pair_v.data()),
                       sizeof(std::pair<se_t, std::pair<se_t, se_t>>) * nsubsets);
  if (!crecord_stream.good()) {
    std::cerr << "Writing subset enumerations has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  crecord_stream.close();
}

vec<node_sptr_t> Record::decode_shash(sh_t sh)
{
  vec<node_sptr_t> subset_v;
  std::queue<sh_t> qsubset;
  node_sptr_t nd;
  subset_sptr_t subset;
  qsubset.push(sh);
  while (!qsubset.empty()) {
    sh = qsubset.front();
    qsubset.pop();
    if (shash_to_node.contains(sh)) {
      nd = shash_to_node[sh];
      subset_v.push_back(nd);
    } else {
      subset = shash_to_subset[sh];
      qsubset.push(subset->ch);
      qsubset.push(subset->sh - subset->ch - subset->nonce);
    }
  }
  return subset_v;
}

vec<node_sptr_t> CRecord::decode_senc(se_t se)
{
  vec<node_sptr_t> subset_v;
  std::queue<se_t> qsubset;
  qsubset.push(se);
  while (!qsubset.empty()) {
    se = qsubset.front();
    qsubset.pop();
    if (se <= senc_to_node.size()) {
      subset_v.push_back(senc_to_node[se]);
    } else {
      qsubset.push(senc_to_psenc[se].first);
      qsubset.push(senc_to_psenc[se].second);
    }
  }
  return subset_v;
}

bool CRecord::check_compatible(crecord_sptr_t crecord) { return true; } // TODO: Implement this.

void CRecord::merge(crecord_sptr_t crecord)
{ // TODO: Make sure that nodes are mapped correctly.
  for (auto const& [se, pse] : crecord->senc_to_psenc) {
    senc_to_psenc[se] = pse; // TODO: Handle collisions across partials.
  }
  nsubsets = senc_to_psenc.size();
}
