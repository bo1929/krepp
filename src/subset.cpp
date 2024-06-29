#include "subset.hpp"

Record::Record(node_sptr_t nd) {
  subtree_root = nd;
  tree = subtree_root->get_tree();
  sh_to_node[nd->get_shash()] = nd;
  nd->get_tree()->set_subtree(nd);
  node_sptr_t nd_curr;
  while (nd_curr = nd->get_tree()->next_post_order()) {
    sh_to_node[nd_curr->get_shash()] = nd_curr;
  }
  nd->get_tree()->reset_traversal();
  if (subtree_root == subtree_root->get_tree()->get_root()) {
    while (tree_with_collision()) {
      rehash_tree();
    }
  }
}

Record::Record(record_sptr_t source1, record_sptr_t source2) {
  union_record(source1);
  union_record(source2);
  subtree_root =
      Tree::compute_lca(source1->subtree_root, source2->subtree_root);
  tree = subtree_root->get_tree();
}

bool Record::tree_with_collision() {
  bool collision_free = true;
  subtree_root->get_tree()->set_subtree(subtree_root);
  node_sptr_t nd_curr;
  while (collision_free &&
         (nd_curr = subtree_root->get_tree()->next_post_order())) {
    if (sh_to_node.find(nd_curr->get_shash()) == sh_to_node.end() ||
        (sh_to_node.at(nd_curr->get_shash()) != nd_curr) ||
        (!nd_curr->get_shash())) {
      collision_free = false;
    }
  }
  subtree_root->get_tree()->reset_traversal();
  return !collision_free;
}

void Record::rehash_tree() {
  sh_to_node.clear();
  sh_to_subset.clear();
  subtree_root->get_tree()->set_subtree(subtree_root);
  node_sptr_t nd_curr;
  while (nd_curr = subtree_root->get_tree()->next_post_order()) {
    if (nd_curr->check_leaf())
      nd_curr->set_shash(Subset::rehash(nd_curr->get_shash()));
    else
      nd_curr->set_shash(nd_curr->sum_children_shash());
    sh_to_node[nd_curr->get_shash()] = nd_curr;
  }
  subtree_root->get_tree()->reset_traversal();
}

void Record::union_record(record_sptr_t source) {
#pragma omp critical
  {
    // TODO: check conflicts and resolve.
    // TODO: check if either of the records is empty.
    sh_to_node.insert(source->sh_to_node.begin(), source->sh_to_node.end());
    sh_to_subset.insert(source->sh_to_subset.begin(),
                        source->sh_to_subset.end());
    subtree_root = Tree::compute_lca(subtree_root, source->subtree_root);
    tree = subtree_root->get_tree();
  }
}

void Record::add_subset(subset_sptr_t new_subset) {
#pragma omp critical
  {
    // TODO: check collisions and resolve.
    if (sh_to_node.find(new_subset->shash) == sh_to_node.end()) {
      sh_to_subset[new_subset->shash] = new_subset;
    } else {
      // TODO: Collision with a node or a real node.
    }
  }
}

Subset::Subset(sh_t shash1, sh_t shash2, record_sptr_t record) {
#pragma omp critical
  {
    auto &node_map = record->sh_to_node;
    auto &subset_map = record->sh_to_subset;
    if ((subset_map.find(shash1) != subset_map.end()) &&
        (subset_map.find(shash2) != subset_map.end())) {
      card = subset_map[shash1]->card + subset_map[shash2]->card;
      shash = subset_map[shash1]->shash + subset_map[shash2]->shash;
      chash = subset_map[shash1]->card > subset_map[shash2]->card
                  ? subset_map[shash2]->shash
                  : subset_map[shash1]->shash;
    } else if ((node_map.find(shash1) != node_map.end()) &&
               (node_map.find(shash2) != node_map.end())) {
      card = node_map[shash1]->card + node_map[shash2]->card;
      shash = node_map[shash1]->shash + node_map[shash2]->shash;
      chash = node_map[shash1]->card > node_map[shash2]->card
                  ? node_map[shash2]->shash
                  : node_map[shash1]->shash;
    } else if ((subset_map.find(shash1) != subset_map.end()) &&
               (node_map.find(shash2) != node_map.end())) {
      card = subset_map[shash1]->card + node_map[shash2]->card;
      shash = subset_map[shash1]->shash + node_map[shash2]->shash;
      chash = subset_map[shash1]->card > node_map[shash2]->card
                  ? node_map[shash2]->shash
                  : subset_map[shash1]->shash;
    } else if ((node_map.find(shash1) != node_map.end()) &&
               (subset_map.find(shash2) != subset_map.end())) {
      card = node_map[shash1]->card + subset_map[shash2]->card;
      shash = node_map[shash1]->shash + subset_map[shash2]->shash;
      chash = node_map[shash1]->card > subset_map[shash2]->card
                  ? subset_map[shash2]->shash
                  : node_map[shash1]->shash;
    } else {
      std::cerr << "Cannot constuct subset, given record lacks the partition!";
      std::quick_exit(EXIT_FAILURE);
    }
    // TODO: Consider adding another subsets to the record for further pruning.
  }
}
