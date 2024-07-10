#include "record.hpp"

Record::Record(node_sptr_t nd)
{
  subtree_root = nd;
  tree = subtree_root->get_tree();
  sh_to_node[nd->get_shash()] = nd;
  nd->get_tree()->set_subtree(nd);
  node_sptr_t nd_curr;
  while (nd_curr = nd->get_tree()->next_post_order()) {
    sh_to_node[nd_curr->get_shash()] = nd_curr;
    /* sh_to_subset[nd_curr->get_shash()] = // TODO: Consider an alternative. */
    /*   std::make_shared<Subset>(nd_curr->get_shash(), 0, nd_curr->get_card()); */
  }
  nd->get_tree()->reset_traversal();
  if (subtree_root == subtree_root->get_tree()->get_root()) {
    while (tree_with_collision()) {
      rehash_tree();
    }
  }
}

Record::Record(record_sptr_t source1, record_sptr_t source2)
{
  union_record(source1);
  union_record(source2);
  subtree_root = Tree::compute_lca(source1->subtree_root, source2->subtree_root);
  tree = subtree_root->get_tree();
}

bool Record::tree_with_collision()
{
  bool collision_free = true;
#pragma omp critical(RecordLock)
  {
    subtree_root->get_tree()->set_subtree(subtree_root);
    node_sptr_t nd_curr;
    while (collision_free && (nd_curr = subtree_root->get_tree()->next_post_order())) {
      if (sh_to_node.find(nd_curr->get_shash()) == sh_to_node.end() ||
          (sh_to_node.at(nd_curr->get_shash()) != nd_curr) || (!nd_curr->get_shash())) {
        collision_free = false;
      }
    }
    subtree_root->get_tree()->reset_traversal();
  }
  return !collision_free;
}

void Record::rehash_tree()
{
#pragma omp critical(RecordLock)
  {
    sh_to_node.clear();
    sh_to_subset.clear();
    subtree_root->get_tree()->set_subtree(subtree_root);
    node_sptr_t nd_curr;
    while (nd_curr = subtree_root->get_tree()->next_post_order()) {
      if (nd_curr->check_leaf()) {
        nd_curr->set_shash(Subset::rehash(nd_curr->get_shash()));
      } else {
        nd_curr->set_shash(nd_curr->sum_children_shash());
      }
      sh_to_node[nd_curr->get_shash()] = nd_curr;
      /* sh_to_subset[nd_curr->get_shash()] = // TODO: Consider an alternative. */
      /*   std::make_shared<Subset>(nd_curr->get_shash(), 0, nd_curr->get_card()); */
    }
    subtree_root->get_tree()->reset_traversal();
  }
}

void Record::union_record(record_sptr_t source)
{
#pragma omp critical(RecordLock)
  {
    // TODO: check conflicts and resolve.
    // TODO: check if either of the records is empty.
    sh_to_node.insert(source->sh_to_node.begin(), source->sh_to_node.end()); // TODO: Maybe remove.
    sh_to_subset.insert(source->sh_to_subset.begin(), source->sh_to_subset.end());
    subtree_root = Tree::compute_lca(subtree_root, source->subtree_root);
    tree = subtree_root->get_tree();
  }
}

void Record::add_subset(subset_sptr_t new_subset)
{
#pragma omp critical(RecordLock)
  {
    // TODO: check collisions and resolve.
    if (sh_to_node.find(new_subset->shash) ==
        sh_to_node.end()) { // TODO:Consider a faster alternative.
      sh_to_subset[new_subset->shash] = new_subset;
    } else {
      // TODO: Collision with a node or a real node.
    }
  }
}

Subset::Subset(sh_t shash1, sh_t shash2, record_sptr_t record)
{
#pragma omp critical(RecordLock)
  { // TODO: Consider a faster alternative.
    auto& subset_map = record->sh_to_subset;
    auto& node_map = record->sh_to_node;
    if ((subset_map.find(shash1) != subset_map.end()) &&
        (subset_map.find(shash2) != subset_map.end())) {
      card = subset_map[shash1]->card + subset_map[shash2]->card;
      shash = subset_map[shash1]->shash + subset_map[shash2]->shash;
      chash = subset_map[shash1]->card > subset_map[shash2]->card ? subset_map[shash2]->shash
                                                                  : subset_map[shash1]->shash;
    } else if ((node_map.find(shash1) != node_map.end()) &&
               (node_map.find(shash2) != node_map.end())) {
      card = node_map[shash1]->card + node_map[shash2]->card;
      shash = node_map[shash1]->shash + node_map[shash2]->shash;
      chash = node_map[shash1]->card > node_map[shash2]->card ? node_map[shash2]->shash
                                                              : node_map[shash1]->shash;
    } else if ((subset_map.find(shash1) != subset_map.end()) &&
               (node_map.find(shash2) != node_map.end())) {
      card = subset_map[shash1]->card + node_map[shash2]->card;
      shash = subset_map[shash1]->shash + node_map[shash2]->shash;
      chash = subset_map[shash1]->card > node_map[shash2]->card ? node_map[shash2]->shash
                                                                : subset_map[shash1]->shash;
    } else if ((node_map.find(shash1) != node_map.end()) &&
               (subset_map.find(shash2) != subset_map.end())) {
      card = node_map[shash1]->card + subset_map[shash2]->card;
      shash = node_map[shash1]->shash + subset_map[shash2]->shash;
      chash = node_map[shash1]->card > subset_map[shash2]->card ? subset_map[shash2]->shash
                                                                : node_map[shash1]->shash;
    } else {
      std::cerr << "Cannot constuct subset, given record lacks the partition!";
      std::quick_exit(EXIT_FAILURE);
    }
    // TODO: Consider adding another subsets to the record for further pruning.
  }
}

void Record::make_compact()
{ // TODO: Implement pruning, probably at this point.
#pragma omp critical(RecordLock)
  {
    se_t limit_senum = std::numeric_limits<se_t>::max();
    se_t curr_senum = 1;
    tree->reset_traversal();
    node_sptr_t nd_curr;
    while (nd_curr = tree->next_post_order()) {
      if (curr_senum < limit_senum) {
        sh_to_se[nd_curr->get_shash()] = curr_senum;
        curr_senum++;
      } else {
        std::cerr << "The current se_t size is too small to fit all nodes!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    tree->reset_traversal();
    for (auto& [shash, subset_sptr] : sh_to_subset) {
      if (curr_senum < limit_senum) {
        sh_to_se[shash] = curr_senum;
        curr_senum++;
      }
    }
  }
}

CRecord::CRecord(record_sptr_t record)
{
  record->make_compact();
  tree = record->get_tree();
  se_t curr_senum = 1;
  tree->reset_traversal();
  node_sptr_t nd_curr;
  while (nd_curr = tree->next_post_order()) {
    se_to_node[curr_senum] = nd_curr;
    curr_senum++;
  }
  tree->reset_traversal();
  for (auto& [shash, subset_sptr] : record->sh_to_subset) {
    se_to_pse[record->sh_to_se[shash]] =
      std::make_pair(record->sh_to_se[subset_sptr->chash],
                     record->sh_to_se[subset_sptr->shash - subset_sptr->chash]);
  }
  se_to_pse[0] = std::make_pair(0, 0);
  nsubsets = se_to_pse.size();
}

void CRecord::print_info()
{
  std::cout << "Total number of subsets excluding nodes: " << nsubsets << std::endl;
  std::cout << "Number of nodes: " << se_to_node.size() << std::endl;
  for (auto [se, nd] : se_to_node) {
    std::cout << se << ": " << nd->get_name() << "(" << nd->get_card() << ")" << std::endl;
  }
  for (auto [se, pse] : se_to_pse) {
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
    se_to_node[curr_senum] = nd_curr;
    curr_senum++;
  }
  tree->reset_traversal();
  nsubsets = 0;
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
    se_to_pse =
      std::unordered_map<se_t, std::pair<se_t, se_t>>(pse_pair_v.begin(), pse_pair_v.end());
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
  std::vector<std::pair<se_t, std::pair<se_t, se_t>>> pse_pair_v(se_to_pse.begin(),
                                                                 se_to_pse.end());
  crecord_stream.write(reinterpret_cast<char*>(pse_pair_v.data()),
                       sizeof(std::pair<se_t, std::pair<se_t, se_t>>) * nsubsets);
  if (!crecord_stream.good()) {
    std::cerr << "Writing subset enumerations has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  crecord_stream.close();
}

vec<se_t> CRecord::decode_se(se_t se)
{
  vec<se_t> subset;
  std::queue<se_t> q;
  q.push(se);
  while (!q.empty()) {
    se = q.front();
    q.pop();
    if (se_to_node.find(se) != se_to_node.end()) {
      subset.push_back(se);
    } else {
      q.push(se_to_pse[se].first);
      q.push(se_to_pse[se].second);
    }
  }
  return subset;
}

bool CRecord::check_compatible(crecord_sptr_t crecord) { return true; } // TODO: Implement this.

void CRecord::merge(crecord_sptr_t crecord)
{ // TODO: Make sure that nodes are mapped correctly.
  for (auto const& [se, pse] : crecord->se_to_pse) {
    se_to_pse[se] = pse; // TODO: Handle collisions across partials.
  }
  nsubsets = se_to_pse.size();
}
