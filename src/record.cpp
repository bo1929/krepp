#include "record.hpp"
#include <cstdint>

Record::Record(node_sptr_t nd)
{
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
      if (nd_curr->check_leaf())
        nd_curr->set_shash(Subset::rehash(nd_curr->get_shash()));
      else
        nd_curr->set_shash(nd_curr->sum_children_shash());
      sh_to_node[nd_curr->get_shash()] = nd_curr;
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
    sh_to_node.insert(source->sh_to_node.begin(), source->sh_to_node.end());
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
    if (sh_to_node.find(new_subset->shash) == sh_to_node.end()) {
      sh_to_subset[new_subset->shash] = new_subset;
    } else {
      // TODO: Collision with a node or a real node.
    }
  }
}

Subset::Subset(sh_t shash1, sh_t shash2, record_sptr_t record)
{
#pragma omp critical(RecordLock)
  {
    auto& node_map = record->sh_to_node;
    auto& subset_map = record->sh_to_subset;
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

void Record::save(std::filesystem::path library_dir, std::string suffix)
{
#pragma omp critical(RecordLock)
  {
    std::ofstream subset_stream(library_dir / ("subset_records" + suffix), std::ofstream::binary);
    uint64_t nsubsets = static_cast<uint64_t>(sh_to_subset.size());
    subset_stream.write(reinterpret_cast<char*>(&nsubsets), sizeof(uint64_t));
    for (auto& [shash, subset_sptr] : sh_to_subset) {
      subset_stream.write(reinterpret_cast<char*>(&subset_sptr->shash), sizeof(sh_t));
      subset_stream.write(reinterpret_cast<char*>(&subset_sptr->chash), sizeof(sh_t));
      subset_stream.write(reinterpret_cast<char*>(&subset_sptr->card), sizeof(tuint_t));
    }
    if (!subset_stream.good()) {
      std::cerr << "Writing the subset records of the tree has failed!" << std::endl;
      exit(EXIT_FAILURE);
    }
    subset_stream.close();

    std::ofstream node_stream(library_dir / ("node_records" + suffix), std::ofstream::binary);
    uint32_t name_length;
    uint64_t nnodes = static_cast<uint64_t>(sh_to_node.size());
    node_stream.write(reinterpret_cast<char*>(&nnodes), sizeof(uint64_t));
    for (auto& [shash, node_sptr] : sh_to_node) {
      name_length = static_cast<uint32_t>(node_sptr->name.length());
      node_stream.write(reinterpret_cast<char*>(&node_sptr->shash), sizeof(sh_t));
      node_stream.write(reinterpret_cast<char*>(&name_length), sizeof(uint32_t));
      node_stream.write(node_sptr->name.c_str(), name_length);
    }
    if (!node_stream.good()) {
      std::cerr << "Writing the node records of the tree has failed!" << std::endl;
      exit(EXIT_FAILURE);
    }
    node_stream.close();
  }
}

void Record::load(std::filesystem::path library_dir, std::string suffix)
{
#pragma omp critical(RecordLock)
  {
    std::filesystem::path subset_path = library_dir / ("subset_records" + suffix);
    std::ifstream subset_stream(subset_path, std::ifstream::binary);
    if (!subset_stream.is_open()) {
      std::cerr << "Failed to open " << subset_path << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sh_t shash, chash, card;
      uint64_t nsubsets;
      subset_stream.read(reinterpret_cast<char*>(&nsubsets), sizeof(uint64_t));
      for (uint64_t six = 0; six < nsubsets; ++six) {
        subset_stream.read(reinterpret_cast<char*>(&shash), sizeof(sh_t));
        subset_stream.read(reinterpret_cast<char*>(&chash), sizeof(sh_t));
        subset_stream.read(reinterpret_cast<char*>(&card), sizeof(tuint_t));
        sh_to_subset.emplace(shash, std::make_shared<Subset>(shash, chash, card));
      }
    }
    if (!subset_stream.good()) {
      std::cerr << "Reading subset records of the tree has failed!" << std::endl;
      exit(EXIT_FAILURE);
    }
    subset_stream.close();

    std::filesystem::path node_path = library_dir / ("node_records" + suffix);
    std::ifstream node_stream(node_path, std::ifstream::binary);
    if (!node_stream.is_open()) {
      std::cerr << "Failed to open " << node_path << std::endl;
      exit(EXIT_FAILURE);
    } else {
      sh_t shash;
      uint32_t name_length;
      std::string name;
      uint64_t nnodes;
      node_stream.read(reinterpret_cast<char*>(&nnodes), sizeof(uint64_t));
      for (uint64_t nix = 0; nix < nnodes; ++nix) {
        node_stream.read(reinterpret_cast<char*>(&shash), sizeof(sh_t));
        node_stream.read(reinterpret_cast<char*>(&name_length), sizeof(uint32_t));
        name.resize(name_length);
        node_stream.read(name.data(), name_length);
        if (!(sh_to_node.find(shash) != sh_to_node.end() && sh_to_node[shash]->name == name)) {
          std::cerr << "Inconsistency between loaded and assigned hash values." << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    if (!node_stream.good()) {
      std::cerr << "Reading node records of the tree has failed!" << std::endl;
      exit(EXIT_FAILURE);
    }
    subset_stream.close();
  }
}
