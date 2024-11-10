#include "index.hpp"

void Index::add_partial_index(std::string suffix)
{
  std::filesystem::path metadata_path = index_dir / ("metadata" + suffix);
  std::ifstream metadata_stream(metadata_path, std::ifstream::binary);
  uint8_t k_lshf, w, h_lshf;
  uint32_t m_lshf, r, nrows_partial;
  bool frac;
  metadata_stream.read(reinterpret_cast<char*>(&k_lshf), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&h_lshf), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&m_lshf), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
  metadata_stream.read(reinterpret_cast<char*>(&nrows_partial), sizeof(uint32_t));
  vec<uint8_t> ppos_v(h_lshf), npos_v(k_lshf - h_lshf);
  metadata_stream.read(reinterpret_cast<char*>(ppos_v.data()), ppos_v.size() * sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(npos_v.data()), npos_v.size() * sizeof(uint8_t));
  if (!metadata_stream.good()) {
    std::cerr << "Reading the metadata for the partial index has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  metadata_stream.close();

  lshf_sptr_t curr_lshf = std::make_shared<LSHF>(m_lshf, ppos_v, npos_v);

#pragma omp critical
  {
    if (curr_lshf->check_compatible(lshf)) {
      lshf = !lshf ? curr_lshf : lshf;
      k = k_lshf;
      h = h_lshf;
      m = m_lshf;
      nrows = pow(2, 2 * h);
    } else {
      std::cerr << "Partial libraries have incompatible hash functions." << std::endl;
      exit(EXIT_FAILURE);
    }
  }

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
      }
    } else {
      r_to_flatht[r] = curr_flatht;
    }
  }
}

std::vector<cmer_t>::const_iterator Index::get_first(uint32_t rix)
{
  return rix < m ? (r_to_flatht[rix % m])->begin() : (r_to_flatht[rix % m])->at(rix / m - 1);
}

std::vector<cmer_t>::const_iterator Index::get_next(uint32_t rix)
{
  return rix < nrows ? (r_to_flatht[rix % m])->at(rix / m) : (r_to_flatht[rix % m])->end();
}

crecord_sptr_t Index::get_crecord(uint32_t rix)
{
  return rix < nrows ? (r_to_flatht[rix % m])->get_crecord()
                     : (r_to_flatht[rix % m])->get_crecord();
}
