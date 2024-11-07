#include "record.hpp"

Record::Record(tree_sptr_t tree)
  : tree(tree)
{
  node_sptr_t nd_curr;
  sh_t ch;
  tree->reset_traversal();
  while (nd_curr = tree->next_post_order()) {
    sh_to_node[nd_curr->get_sh()] = nd_curr;
    if (nd_curr->check_leaf()) {
      ch = 0;
    } else {
      ch = (*nd_curr->get_children())->get_sh();
    }
    sh_to_subset[nd_curr->get_sh()] =
      std::make_shared<Subset>(nd_curr->get_sh(), ch, nd_curr->get_card());
  }
  tree->reset_traversal();
  while (check_tree_collision()) {
    rehash_tree();
  }
  sh_to_subset[0] = std::make_shared<Subset>(0, 0, 0);
}

Record::Record(record_sptr_t source1, record_sptr_t source2)
{
  union_record(source1);
  union_record(source2);
  node_sptr_t root =
    Tree::compute_lca(source1->get_tree()->get_root(), source2->get_tree()->get_root());
  tree = root->get_tree();
}

bool Record::check_tree_collision()
{
  bool collision_free = true;
  node_sptr_t nd_curr;
  sh_t sh_curr;
  tree->reset_traversal();
  while (collision_free && (nd_curr = tree->next_post_order())) {
    sh_curr = nd_curr->get_sh();
    if (!sh_to_node.contains(sh_curr) || (sh_to_node[sh_curr] != nd_curr) || (!sh_curr)) {
      collision_free = false;
    }
  }
  tree->reset_traversal();
  return !collision_free;
}

void Record::rehash_tree()
{
  sh_to_node.clear();
  sh_to_subset.clear();
  node_sptr_t nd_curr;
  sh_t ch;
  sh_t ah;
  tree->reset_traversal();
  while (nd_curr = tree->next_post_order()) {
    if (nd_curr->check_leaf()) {
      nd_curr->set_sh(++ah + Subset::rehash(nd_curr->get_sh()));
    } else {
      nd_curr->set_sh(nd_curr->sum_children_sh());
    }
    sh_to_node[nd_curr->get_sh()] = nd_curr;
    if (nd_curr->check_leaf()) {
      ch = 0;
    } else {
      ch = (*nd_curr->get_children())->get_sh();
    }
    sh_to_subset[nd_curr->get_sh()] =
      std::make_shared<Subset>(nd_curr->get_sh(), ch, nd_curr->get_card());
  }
  tree->reset_traversal();
}

sh_t Record::add_subset(sh_t sh1, sh_t sh2)
{
  if (!(sh_to_subset.contains(sh1) && sh_to_subset.contains(sh2))) {
    std::cerr << "Cannot make the subset for the partition: (" << sh1 << ", " << sh2 << ")\n";
    std::quick_exit(EXIT_FAILURE);
  }
  subset_sptr_t subset1 = sh_to_subset[sh1];
  subset_sptr_t subset2 = sh_to_subset[sh2];
  sh_t sh = sh1 + sh2;
  sh_t nonce = 0;
  while (check_subset_collision(sh + nonce, subset1, subset2)) {
    nonce = Subset::rehash(nonce++ * sh1 * sh2);
  }
  sh += nonce;
  if (!sh_to_subset.contains(sh)) {
    sh_to_subset[sh] =
      std::make_shared<Subset>(sh,
                               subset1->card > subset2->card ? subset1->sh : subset2->sh,
                               subset1->card + subset2->card,
                               nonce);
  }
  return sh;
}

void Record::union_record(record_sptr_t source)
{
  // TODO: check conflicts and resolve.
  // TODO: check if either of the records is empty.
  sh_to_node.insert(source->sh_to_node.begin(), source->sh_to_node.end());
  sh_to_subset.insert(source->sh_to_subset.begin(), source->sh_to_subset.end());
  node_sptr_t root = Tree::compute_lca(tree->get_root(), source->get_tree()->get_root());
  tree = root->get_tree();
}

bool Record::check_subset_collision(sh_t sh, subset_sptr_t subset1, subset_sptr_t subset2)
{
  if (sh == 0) {
    return true;
  } else if (sh_to_subset.contains(sh)) {
    subset_sptr_t subset = sh_to_subset[sh];
    if (subset->ch == 0) {
      return true;
    } else if ((subset->ch == subset1->sh || subset->ch == subset2->sh)) {
      return false;
    } else {
      return true;
    }
  } else {
    return false;
  }
}

void Record::make_compact()
{ // TODO: Implement pruning, probably at this point.
  se_t limit_senum = std::numeric_limits<se_t>::max();
  se_t curr_senum = 1;
  node_sptr_t nd_curr;
  tree->reset_traversal();
  while (nd_curr = tree->next_post_order()) {
    if (curr_senum < limit_senum) {
      sh_to_se[nd_curr->get_sh()] = curr_senum++;
    } else {
      std::cerr << "The current se_t size is too small to fit all nodes!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  tree->reset_traversal();
  for (auto& [sh, subset] : sh_to_subset) {
    if (curr_senum < limit_senum) {
      curr_senum += static_cast<se_t>(sh_to_se.try_emplace(sh, curr_senum).second);
    } else {
      std::cerr << "The current se_t size is too small to fit all subsets observed!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  sh_to_se[0] = 0;
}

CRecord::CRecord(record_sptr_t record)
{
  record->make_compact();
  tree = record->get_tree();
  se_t curr_senum = 1;
  node_sptr_t nd_curr;
  tree->reset_traversal();
  nsubsets = record->sh_to_se.size() + 1;
  nnodes = record->sh_to_node.size() + 1;
  se_to_pse.resize(nsubsets);
  se_to_rho.resize(nnodes);
  while (nd_curr = tree->next_post_order()) {
    se_to_rho[nd_curr->get_se()] = record->sh_to_rho[nd_curr->get_sh()];
  }
  tree->reset_traversal();
  for (auto& [sh, subset] : record->sh_to_subset) {
    se_to_pse[record->sh_to_se[sh]] = std::make_pair(
      record->sh_to_se[subset->ch], record->sh_to_se[sh - subset->ch - subset->nonce]);
  }
  se_to_pse[0] = std::make_pair(0, 0);
}

void CRecord::print_info()
{
  std::cout << "Total number of subsets excluding nodes: " << nsubsets << std::endl;
  std::cout << "Number of nodes: " << tree->get_nnodes() << std::endl;
  for (uint32_t se = 1; se < nnodes; ++se) {
    node_sptr_t nd = tree->get_node(se);
    std::cout << se << ": " << nd->get_name() << "(" << nd->get_card() << ")" << std::endl;
  }
  for (uint32_t se = 1; se < nsubsets; ++se) {
    pse_t& pse = se_to_pse[se];
    std::cout << se << ": " << pse.first << "+" << pse.second << std::endl;
  }
}

CRecord::CRecord(tree_sptr_t tree)
  : tree(tree)
{
  node_sptr_t nd_curr;
  se_t curr_senum = 1;
  tree->reset_traversal();
  nnodes = tree->get_nnodes() + 1;
  nsubsets = nnodes;
  se_to_rho.resize(nnodes);
  while (nd_curr = tree->next_post_order()) {
    se_to_rho[curr_senum] = 0;
    curr_senum++;
  }
  tree->reset_traversal();
}

void CRecord::load(std::filesystem::path library_dir, std::string suffix)
{
  std::filesystem::path crecord_path = library_dir / ("crecord" + suffix);
  std::ifstream crecord_stream(crecord_path, std::ifstream::binary);
  if (!crecord_stream.is_open()) {
    std::cerr << "Failed to open " << crecord_path << std::endl;
    exit(EXIT_FAILURE);
  } else {
    crecord_stream.read(reinterpret_cast<char*>(&nnodes), sizeof(se_t));
    crecord_stream.read(reinterpret_cast<char*>(&nsubsets), sizeof(se_t));
    se_to_pse.resize(nsubsets);
    crecord_stream.read(reinterpret_cast<char*>(se_to_pse.data()), sizeof(pse_t) * nsubsets);
    se_to_rho.resize(nnodes);
    crecord_stream.read(reinterpret_cast<char*>(se_to_rho.data()), sizeof(double) * nnodes);
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
  crecord_stream.write(reinterpret_cast<char*>(&nnodes), sizeof(se_t));
  crecord_stream.write(reinterpret_cast<char*>(&nsubsets), sizeof(se_t));
  crecord_stream.write(reinterpret_cast<char*>(se_to_pse.data()), sizeof(pse_t) * nsubsets);
  crecord_stream.write(reinterpret_cast<char*>(se_to_rho.data()), sizeof(double) * nnodes);
  if (!crecord_stream.good()) {
    std::cerr << "Writing subset enumerations has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  crecord_stream.close();
}

void Record::decode_sh(sh_t sh, vec<node_sptr_t>& subset_v)
{
  std::queue<sh_t> qsubset;
  qsubset.push(sh);
  subset_sptr_t subset;
  while (!qsubset.empty()) {
    sh = qsubset.front();
    qsubset.pop();
    if (sh_to_node.contains(sh)) {
      subset_v.push_back(sh_to_node[sh]);
    } else {
      subset = sh_to_subset[sh];
      qsubset.push(subset->ch);
      qsubset.push(subset->sh - subset->ch - subset->nonce);
    }
  }
}

void CRecord::decode_se(se_t se, vec<node_sptr_t>& subset_v)
{
  std::queue<se_t> qsubset;
  qsubset.push(se);
  pse_t pse;
  while (!qsubset.empty()) {
    se = qsubset.front();
    qsubset.pop();
    if (tree->get_node(se)) {
      subset_v.push_back(tree->get_node(se));
    } else {
      pse = se_to_pse[se];
      qsubset.push(pse.first);
      qsubset.push(pse.second);
    }
  }
}

void CRecord::display_info(uint32_t r, vec<uint64_t>& se_to_count)
{
  std::cout << r << "\t#NODES\t" << nnodes << std::endl;
  std::cout << r << "\t#COLORS\t" << nsubsets << std::endl;
  vec<uint64_t> se_to_outdegree(se_to_pse.size());
  for (int ix = 1; ix < se_to_pse.size(); ++ix) {
    se_to_outdegree[se_to_pse[ix].first]++;
    se_to_outdegree[se_to_pse[ix].second]++;
  }
  flat_phmap<uint64_t, uint32_t> outdegree_hist;
  for (uint64_t ix = 1; ix < se_to_outdegree.size(); ++ix) {
    outdegree_hist[se_to_outdegree[ix]]++;
  }
  flat_phmap<uint64_t, uint32_t> count_hist;
  for (uint64_t ix = 1; ix < se_to_pse.size(); ++ix) {
    count_hist[se_to_count[ix]]++;
  }
  for (auto const& [key, val] : count_hist) {
    std::cout << r << "\t#COLOR_COUNT" << key << "\t" << val << std::endl;
  }
  for (auto const& [key, val] : outdegree_hist) {
    std::cout << r << "\t#COLOR_OUTDEGREE\t" << key << "\t" << val << std::endl;
  }
}
