#include "index.hpp"

void Index::generate_partial_tree(std::string suffix)
{
  wbackbone = false;
  std::ifstream reflist_file(index_dir / ("reflist" + suffix));
  std::string name;
  std::vector<std::string> names_v;
  if (reflist_file.is_open()) {
    while (std::getline(reflist_file, name)) {
      names_v.push_back(name);
    }
    reflist_file.close();
  } else {
    error_exit("Unable to open reference list file for an index without a tree.");
  }
  tree_sptr_t curr_tree = std::make_shared<Tree>();
  curr_tree->generate_tree(names_v);
#pragma omp critical
  {
    if (curr_tree->check_compatible(tree)) {
      tree = !tree ? curr_tree : tree;
    } else {
      error_exit("Partial libraries are based on different references.");
    }
  }
}

void Index::load_partial_tree(std::string suffix)
{
  wbackbone = true;
  tree_sptr_t curr_tree = std::make_shared<Tree>();
  std::filesystem::path nwk_path = index_dir / ("tree" + suffix);
  std::ifstream tree_stream(nwk_path);
  if (!tree_stream.is_open()) {
    error_exit(std::string("Failed to open ") + nwk_path.string());
  }
  curr_tree->load(tree_stream);
  CHECK_STREAM_OR_EXIT(tree_stream, "Failed to read the backbone tree of a partial index!");
  tree_stream.close();
#pragma omp critical
  {
    if (curr_tree->check_compatible(tree)) {
      tree = !tree ? curr_tree : tree;
    } else {
      error_exit("Partial libraries are based on different trees!");
    }
  }
}

void Index::load_partial_index(std::string suffix)
{ // TODO: Split for each file (e.g, metadata, crecord etc.)
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
  CHECK_STREAM_OR_EXIT(metadata_stream, "Failed to read the metadata of a partial skecth!");
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
      error_exit("Partial libraries have incompatible hash functions!");
    }
  }

  crecord_sptr_t curr_crecord;
  flatht_sptr_t curr_flatht;
#pragma omp critical
  {
    curr_crecord = std::make_shared<CRecord>(tree);
    curr_flatht = std::make_shared<FlatHT>(tree, curr_crecord);
  }

  std::filesystem::path mer_path = index_dir / ("cmer" + suffix);
  std::ifstream mer_stream(mer_path, std::ifstream::binary);
  if (!mer_stream.is_open()) {
    error_exit(std::string("Failed to open ") + mer_path.string());
  }
  std::filesystem::path inc_path = index_dir / ("inc" + suffix);
  std::ifstream inc_stream(inc_path, std::ifstream::binary);
  if (!inc_stream.is_open()) {
    error_exit(std::string("Failed to open ") + inc_path.string());
  }
  curr_flatht->load(mer_stream, inc_stream);
  CHECK_STREAM_OR_EXIT(mer_stream, "Failed to read the k-mer vector of a partial index!");
  mer_stream.close();
  CHECK_STREAM_OR_EXIT(inc_stream, "Failed to read the offset array of a partial index!");
  inc_stream.close();

  std::filesystem::path crecord_path = index_dir / ("crecord" + suffix);
  std::ifstream crecord_stream(crecord_path, std::ifstream::binary);
  if (!crecord_stream.is_open()) {
    error_exit(std::string("Failed to open ") + crecord_path.string());
  }
  curr_crecord->load(crecord_stream);
  CHECK_STREAM_OR_EXIT(crecord_stream, "Failed to read the color array of a partial index!");
  crecord_stream.close();

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

std::pair<vec_cmer_it, vec_cmer_it> Index::bucket_indices(uint32_t rix)
{
  uint32_t rix_res = rix % m;
  uint32_t offset = rix / m;
  if (r_to_numerator[rix_res] > 1) {
    offset = offset * r_to_numerator[rix_res] + rix_res;
  }
  return std::make_pair(r_to_flatht[rix_res]->bucket_start(offset), r_to_flatht[rix_res]->bucket_next(offset));
}

crecord_sptr_t Index::get_crecord(uint32_t rix) { return r_to_flatht[rix % m]->get_crecord(); }

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
