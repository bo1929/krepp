#include "table.hpp"

FlatHT::FlatHT(dynht_sptr_t source)
{
  nkmers = source->nkmers;
  nrows = source->nrows;
  inc_v.resize(nrows);
  cmer_v.reserve(nkmers);
  tree = source->tree;
  crecord = std::make_shared<CRecord>(source->get_record());
  inc_t limit_inc = std::numeric_limits<inc_t>::max();
  inc_t copy_inc;
  inc_t lix = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    copy_inc = std::min(limit_inc, source->mer_vvec[rix].size());
    for (inc_t i = 0; i < copy_inc; ++i) {
      cmer_v.emplace_back(source->conv_mer_cmer(source->mer_vvec[rix][i]));
    }
    lix += copy_inc;
    inc_v[rix] = lix;
    source->mer_vvec[rix].clear();
  }
}

void FlatHT::load(std::filesystem::path library_dir, std::string suffix)
{
  std::filesystem::path mer_path = library_dir / ("cmer" + suffix);
  std::ifstream mer_stream(mer_path, std::ifstream::binary);
  if (!mer_stream.is_open()) {
    std::cerr << "Failed to open " << mer_path << std::endl;
    exit(EXIT_FAILURE);
  } else {
    mer_stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
    cmer_v.resize(nkmers);
    mer_stream.read(reinterpret_cast<char*>(cmer_v.data()), nkmers * sizeof(cmer_t));
    assert(nkmers == cmer_v.size());
  }
  if (!mer_stream.good()) {
    std::cerr << "Reading k-mer vector to has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  mer_stream.close();

  std::filesystem::path inc_path = library_dir / ("inc" + suffix);
  std::ifstream inc_stream(inc_path, std::ifstream::binary);
  if (!inc_stream.is_open()) {
    std::cerr << "Failed to open " << inc_path << std::endl;
    exit(EXIT_FAILURE);
  } else {
    inc_stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
    inc_v.resize(nrows);
    inc_stream.read(reinterpret_cast<char*>(inc_v.data()), nrows * sizeof(inc_t));
    assert(nrows == inc_v.size());
  }
  if (!inc_stream.good()) {
    std::cerr << "Reading index-increment vector has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  inc_stream.close();
}

void FlatHT::save(std::filesystem::path library_dir, std::string suffix)
{
  std::ofstream mer_stream(library_dir / ("cmer" + suffix), std::ofstream::binary);
  mer_stream.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  mer_stream.write(reinterpret_cast<const char*>(cmer_v.data()), sizeof(cmer_t) * nkmers);
  if (!mer_stream.good()) {
    std::cerr << "Writing k-mer vector has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  mer_stream.close();

  std::ofstream inc_stream(library_dir / ("inc" + suffix), std::ofstream::binary);
  inc_stream.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  inc_stream.write(reinterpret_cast<const char*>(inc_v.data()), sizeof(inc_t) * nrows);
  if (!inc_stream.good()) {
    std::cerr << "Writing index-increment vector has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  inc_stream.close();
}

void DynHT::print_info()
{
  update_size_hist();
  std::cout << "size: " << nkmers << "\t";
  for (auto kv : size_hist) {
    std::cout << "H(" << kv.first << ")=" << kv.second << "/";
  }
}

void DynHT::clear_rows()
{
  mer_vvec.clear();
  size_hist.clear();
  nkmers = 0;
}

void DynHT::sort_columns()
{
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if (!mer_vvec[i].empty()) {
      std::sort(mer_vvec[i].begin(), mer_vvec[i].end(), comp_encoding);
    }
  }
}

void DynHT::ensure_sorted_columns()
{
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if ((!mer_vvec[i].empty()) &&
        !std::is_sorted(mer_vvec[i].begin(), mer_vvec[i].end(), comp_encoding)) {
      std::sort(mer_vvec[i].begin(), mer_vvec[i].end(), comp_encoding);
    }
  }
}

void DynHT::update_nkmers()
{
  nkmers = 0;
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    nkmers += mer_vvec[i].size();
  }
}

void DynHT::update_size_hist()
{
  size_hist.clear();
  nkmers = 0;
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    size_hist[mer_vvec[i].size()]++;
    nkmers += mer_vvec[i].size();
  }
}

void DynHT::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if (!mer_vvec[i].empty()) {
      mer_vvec[i].erase(std::unique(mer_vvec[i].begin(), mer_vvec[i].end(), eq_encoding),
                        mer_vvec[i].end());
    }
    nkmers += mer_vvec[i].size();
  }
}

void DynHT::prune_columns(size_t max_size)
{
  nkmers = 0;
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if (mer_vvec[i].size() > max_size) {
      vec<mer_t> tmp_v;
      tmp_v.reserve(max_size);
      std::sample(mer_vvec[i].begin(), mer_vvec[i].end(), std::back_inserter(tmp_v), max_size, gen);
      mer_vvec[i] = std::move(tmp_v);
    }
    nkmers += max_size;
  }
}

void DynHT::union_table(dynht_sptr_t source)
{
  assertm(nrows == source->nrows, "Two tables differ in size.");
  if (source->mer_vvec.empty()) {
    return;
  } else if (mer_vvec.empty()) {
    mer_vvec = std::move(source->mer_vvec);
    size_hist = std::move(source->size_hist);
    nkmers = source->nkmers;
    return;
  } else {
    nkmers = 0;
    for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
      if (!source->mer_vvec[i].empty() && !mer_vvec[i].empty()) {
        DynHT::union_row(mer_vvec[i], source->mer_vvec[i]);
      } else if (!source->mer_vvec[i].empty()) {
        mer_vvec[i] = std::move(source->mer_vvec[i]);
      } else {
      }
      nkmers += mer_vvec[i].size();
    }
  }
}

void DynHT::union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v)
{
  vec<mer_t> temp_v;
  temp_v.reserve(source_v.size() + dest_v.size());
  auto it_dest = dest_v.begin(), it_source = source_v.begin();
  for (; it_source != source_v.end(); ++it_source) {
    while (it_dest != dest_v.end() && it_source->encoding > it_dest->encoding) {
      temp_v.push_back(*it_dest);
      it_dest++;
    }
    while (it_dest != dest_v.end() && it_source->encoding == it_dest->encoding) {
      it_source->sh = record->add_subset(it_dest->sh, it_source->sh);
      it_dest++;
    }
    temp_v.push_back(*it_source);
  }
  for (; it_dest != dest_v.end(); ++it_dest) {
    temp_v.push_back(*it_dest);
  }
  dest_v = std::move(temp_v);
}

/*
void DynHT::union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v)
{
  dest_v.insert(dest_v.end(), source_v.begin(), source_v.end());
  std::inplace_merge(dest_v.begin(),
                     std::next(dest_v.begin(), dest_v.size() - source_v.size()),
                     dest_v.end(),
                     comp_encoding);
  auto it_dest = dest_v.begin(), result = dest_v.begin();
  while (++it_dest != dest_v.end()) {
    if (result->encoding == it_dest->encoding) {
      result->sh = record->add_subset(it_dest->sh, result->sh);
    } else if (++result != it_dest) {
      *result = *it_dest;
    } else {
      result = it_dest;
    }
  }
  if (result != dest_v.end()) {
    dest_v.erase(++result, dest_v.end());
  }
}
*/

void DynHT::fill_table(rseq_sptr_t rs)
{
  mer_vvec.resize(nrows);
  while (rs->read_next_seq() && rs->set_curr_seq()) {
    rs->extract_mers(mer_vvec);
  }
  // update_nkmers();
  sort_columns();
  make_unique();
}
