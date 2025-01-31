#include "index.hpp"

void Index::generate_partial_tree(std::string suffix)
{
  wtree = false;
  std::ifstream reflist_file(index_dir / ("reflist" + suffix));
  std::string name;
  std::vector<std::string> names_v;
  if (reflist_file.is_open()) {
    while (std::getline(reflist_file, name)) {
      names_v.push_back(name);
    }
    reflist_file.close();
  } else {
    std::cerr << "Unable to open reference list file for an index without a tree." << std::endl;
    exit(EXIT_FAILURE);
  }

  tree_sptr_t curr_tree = std::make_shared<Tree>();
  curr_tree->generate_tree(names_v);

#pragma omp critical
  {
    if (curr_tree->check_compatible(tree)) {
      tree = !tree ? curr_tree : tree;
    } else {
      std::cerr << "Partial libraries are based on different references." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

void Index::add_partial_tree(std::string suffix)
{
  wtree = true;
  tree_sptr_t curr_tree = std::make_shared<Tree>();
  curr_tree->load(index_dir, suffix);

#pragma omp critical
  {
    if (curr_tree->check_compatible(tree)) {
      tree = !tree ? curr_tree : tree;
    } else {
      std::cerr << "Partial libraries are based on different trees." << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

void Index::add_partial_index(std::string suffix)
{
  std::filesystem::path metadata_path = index_dir / ("metadata" + suffix);
  std::ifstream metadata_stream(metadata_path, std::ifstream::binary);
  uint8_t k_curr, w, h_curr;
  uint32_t m_curr, r, nrows_partial;
  bool frac;
  metadata_stream.read(reinterpret_cast<char*>(&k_curr), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&h_curr), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&m_curr), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
  metadata_stream.read(reinterpret_cast<char*>(&nrows_partial), sizeof(uint32_t));
  vec<uint8_t> ppos_v(h_curr), npos_v(k_curr - h_curr);
  metadata_stream.read(reinterpret_cast<char*>(ppos_v.data()), ppos_v.size() * sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(npos_v.data()), npos_v.size() * sizeof(uint8_t));
  if (!metadata_stream.good()) {
    std::cerr << "Reading the metadata for the partial index has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  metadata_stream.close();

  lshf_sptr_t curr_lshf = std::make_shared<LSHF>(m_curr, ppos_v, npos_v);

#pragma omp critical
  {
    if (curr_lshf->check_compatible(lshf)) {
      lshf = !lshf ? curr_lshf : lshf;
      k = k_curr;
      h = h_curr;
      m = m_curr;
      nrows = pow(2, 2 * h);
    } else {
      std::cerr << "Partial libraries have incompatible hash functions." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  crecord_sptr_t curr_crecord;
  flatht_sptr_t curr_flatht;

#pragma omp critical
  {
    curr_crecord = std::make_shared<CRecord>(tree);
    curr_flatht = std::make_shared<FlatHT>(tree, curr_crecord);
  }

  curr_crecord->load(index_dir, suffix);
  curr_flatht->load(index_dir, suffix);

#pragma omp critical
  {
    if (frac) {
      for (uint32_t ix = 0; ix <= r; ++ix) {
        r_to_flatht[ix] = curr_flatht;
        r_to_numerator[ix] = r + 1;
      }
    } else {
      r_to_flatht[r] = curr_flatht;
      r_to_numerator[r] = 1;
    }
  }
}

std::vector<cmer_t>::const_iterator Index::bucket_start()
{
  return r_to_flatht[rix_res]->bucket_start(offset);
}

std::vector<cmer_t>::const_iterator Index::bucket_next()
{
  return r_to_flatht[rix_res]->bucket_next(offset);
}

crecord_sptr_t Index::get_crecord() { return r_to_flatht[rix_res]->get_crecord(); }

void Index::make_rho_partial()
{
  double ratio_m = static_cast<double>(r_to_flatht.size()) / static_cast<double>(m);
  flat_phmap<flatht_sptr_t, bool> flatht_to_applied;
  for (auto [r, flatht] : r_to_flatht) {
    flatht_to_applied[flatht] = false;
  }
  for (auto [r, flatht] : r_to_flatht) {
    if (!flatht_to_applied[flatht]) {
      flatht->get_crecord()->apply_rho_coef(ratio_m);
      flatht_to_applied[flatht] = true;
    }
  }
}
