#include "table.hpp"

namespace {
  /* Reads exactly count bytes; a short read means the file was truncated. */
  void read_exact(std::ifstream& stream, void* dst, size_t count, const std::string& what)
  {
    stream.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(count));
    if (static_cast<size_t>(stream.gcount()) != count) {
      error_exit("Truncated " + what);
    }
  }
} // namespace

SFlatHT::SFlatHT(sdynht_sptr_t source)
{
  nkmers = source->nkmers;
  nrows = source->enc_vvec.size();
  inc_ov.resize(nrows);
  enc_ov.reserve(nkmers);
  inc_t limit_inc = std::numeric_limits<inc_t>::max();
  inc_t copy_inc;
  inc_t lix = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    copy_inc = std::min(limit_inc, static_cast<inc_t>(source->enc_vvec[rix].size()));
    for (inc_t i = 0; i < copy_inc; ++i) {
      enc_ov.push_back(source->enc_vvec[rix][i]);
    }
    lix += copy_inc;
    inc_ov[rix] = lix;
    source->enc_vvec[rix].clear();
  }
  enc_vv = enc_ov.data();
  inc_vv = reinterpret_cast<const char*>(inc_ov.data());
}

size_t SFlatHT::load(std::ifstream& sketch_stream)
{
  read_exact(sketch_stream, &nkmers, sizeof(uint64_t), "sketch file");
  enc_ov.resize(nkmers);
  read_exact(sketch_stream, enc_ov.data(), nkmers * sizeof(enc_t), "sketch file");
  assert(nkmers == enc_ov.size());
  read_exact(sketch_stream, &nrows, sizeof(uint32_t), "sketch file");
  inc_ov.resize(nrows);
  read_exact(sketch_stream, inc_ov.data(), nrows * sizeof(inc_t), "sketch file");
  assert(nrows == inc_ov.size());
  enc_vv = enc_ov.data();
  inc_vv = reinterpret_cast<const char*>(inc_ov.data());
  return static_cast<size_t>(sketch_stream.tellg());
}

size_t SFlatHT::load(std::ifstream& sketch_stream, const std::filesystem::path& path)
{
  if (use_mmap) {
    map = krepp::FileMap(path);
    if (map.is_open()) {
      size_t offset = 0;
      uint64_t nk = 0;
      std::memcpy(&nk, map.data() + offset, sizeof(nk));
      offset += sizeof(nk);
      const size_t enc_bytes = static_cast<size_t>(nk) * sizeof(enc_t);
      if (map.size() < offset + enc_bytes + sizeof(uint32_t)) {
        error_exit("Truncated sketch file: " + path.string());
      }
      enc_vv = reinterpret_cast<const enc_t*>(map.data() + offset);
      offset += enc_bytes;
      uint32_t nr = 0;
      std::memcpy(&nr, map.data() + offset, sizeof(nr));
      offset += sizeof(nr);
      if (map.size() < offset + static_cast<size_t>(nr) * sizeof(inc_t) + sizeof(uint32_t) + 3 * sizeof(uint8_t)) {
        error_exit("Truncated sketch file: " + path.string());
      }
      nkmers = nk;
      nrows = nr;
      inc_vv = map.data() + offset;
      offset += static_cast<size_t>(nr) * sizeof(inc_t);
      return offset;
    }
  }
  return load(sketch_stream);
}

void SFlatHT::save(std::ofstream& sketch_stream)
{
  sketch_stream.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  sketch_stream.write(reinterpret_cast<const char*>(enc_vv), sizeof(enc_t) * nkmers);
  sketch_stream.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  sketch_stream.write(inc_vv, sizeof(inc_t) * nrows);
}

FlatHT::FlatHT(dynht_sptr_t source)
{
  nkmers = source->nkmers;
  nrows = source->nrows;
  inc_ov.resize(nrows);
  cmer_ov.reserve(nkmers);
  tree = source->tree;
  crecord = std::make_shared<CRecord>(source->get_record());
  inc_t limit_inc = std::numeric_limits<inc_t>::max();
  inc_t copy_inc;
  inc_t lix = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    copy_inc = std::min(limit_inc, static_cast<inc_t>(source->mer_vvec[rix].size()));
    for (inc_t i = 0; i < copy_inc; ++i) {
      cmer_ov.emplace_back(source->conv_mer_cmer(source->mer_vvec[rix][i]));
    }
    lix += copy_inc;
    inc_ov[rix] = lix;
    source->mer_vvec[rix].clear();
  }
  bind();
}

void FlatHT::bind()
{
  cmer_vv = cmer_ov.data();
  inc_vv = reinterpret_cast<const char*>(inc_ov.data());
}

void FlatHT::load(std::ifstream& mer_stream, std::ifstream& inc_stream)
{
  read_exact(mer_stream, &nkmers, sizeof(uint64_t), "k-mer array");
  cmer_ov.resize(nkmers);
  read_exact(mer_stream, cmer_ov.data(), nkmers * sizeof(cmer_t), "k-mer array");
  assert(nkmers == cmer_ov.size());
  read_exact(inc_stream, &nrows, sizeof(uint32_t), "offset array");
  inc_ov.resize(nrows);
  read_exact(inc_stream, inc_ov.data(), nrows * sizeof(inc_t), "offset array");
  assert(nrows == inc_ov.size());
  bind();
}

void FlatHT::load(const std::filesystem::path& mer_path, const std::filesystem::path& inc_path)
{
  if (!use_mmap) {
    std::ifstream mer_stream(mer_path, std::ifstream::binary);
    std::ifstream inc_stream(inc_path, std::ifstream::binary);
    if (!mer_stream.is_open() || !inc_stream.is_open()) {
      error_exit("Failed to open " + mer_path.string() + " or " + inc_path.string());
    }
    load(mer_stream, inc_stream);
    return;
  }
  mer_map = krepp::FileMap(mer_path);
  inc_map = krepp::FileMap(inc_path);
  if (!mer_map.is_open() || !inc_map.is_open()) {
    // Some filesystems cannot map: fall back to reading.
    std::ifstream mer_stream(mer_path, std::ifstream::binary);
    std::ifstream inc_stream(inc_path, std::ifstream::binary);
    if (!mer_stream.is_open() || !inc_stream.is_open()) {
      error_exit("Failed to open " + mer_path.string() + " or " + inc_path.string());
    }
    load(mer_stream, inc_stream);
    return;
  }
  uint64_t nk = 0;
  std::memcpy(&nk, mer_map.data(), sizeof(nk));
  uint32_t nr = 0;
  std::memcpy(&nr, inc_map.data(), sizeof(nr));
  // Validate before touching the arrays.
  if (mer_map.size() < sizeof(nk) + static_cast<size_t>(nk) * sizeof(cmer_t)) {
    error_exit("Truncated k-mer array in " + mer_path.string());
  }
  if (inc_map.size() < sizeof(nr) + static_cast<size_t>(nr) * sizeof(inc_t)) {
    error_exit("Truncated offset array in " + inc_path.string());
  }
  nkmers = nk;
  nrows = nr;
  cmer_vv = reinterpret_cast<const cmer_t*>(mer_map.data() + sizeof(nk));
  inc_vv = inc_map.data() + sizeof(nr);
}

void FlatHT::save(std::ofstream& mer_stream, std::ofstream& inc_stream)
{
  mer_stream.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  mer_stream.write(reinterpret_cast<const char*>(cmer_vv), sizeof(cmer_t) * nkmers);
  inc_stream.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  inc_stream.write(inc_vv, sizeof(inc_t) * nrows);
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

void SDynHT::sort_columns()
{
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      std::sort(enc_vvec[i].begin(), enc_vvec[i].end());
    }
  }
}

void DynHT::sort_columns()
{
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if (!mer_vvec[i].empty()) {
      std::sort(mer_vvec[i].begin(), mer_vvec[i].end(), comp_encoding);
    }
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

void SDynHT::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      enc_vvec[i].erase(std::unique(enc_vvec[i].begin(), enc_vvec[i].end()), enc_vvec[i].end());
    }
    nkmers += enc_vvec[i].size();
  }
}

void DynHT::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < mer_vvec.size(); ++i) {
    if (!mer_vvec[i].empty()) {
      mer_vvec[i].erase(std::unique(mer_vvec[i].begin(), mer_vvec[i].end(), eq_encoding), mer_vvec[i].end());
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
    nkmers += mer_vvec[i].size();
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
    source->nkmers = 0;
    source->size_hist.clear();
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
  temp_v.reserve(dest_v.size() + source_v.size());

  auto id = dest_v.begin(), is = source_v.begin();
  const auto ed = dest_v.end(), es = source_v.end();

  while (id != ed && is != es) {
    if (id->encoding < is->encoding) {
      temp_v.push_back(*id++);
    } else if (is->encoding < id->encoding) {
      temp_v.push_back(*is++);
    } else {
      mer_t x = *is++;
      while (id != ed && id->encoding == x.encoding) {
        x.sh = record->add_subset(id->sh, x.sh);
        ++id;
      }
      temp_v.push_back(x);
    }
  }
  temp_v.insert(temp_v.end(), id, ed);
  temp_v.insert(temp_v.end(), is, es);

  dest_v.swap(temp_v);
}

void SDynHT::fill_table(uint32_t nrows, rseq_sptr_t rs)
{
  enc_vvec.resize(nrows);
  while (rs->read_next_seq()) {
    if (rs->set_curr_seq()) {
      rs->extract_mers(enc_vvec);
    }
  }
  sort_columns();
  make_unique();
  /* update_nkmers(); */
  rs->compute_rho();
}

void DynHT::fill_table(sh_t sh, rseq_sptr_t rs, bool curr)
{
  mer_vvec.resize(nrows);
  if (curr) {
    rs->reset_estimates();
    rs->extract_mers(mer_vvec, sh);
  } else {
    while (rs->read_next_seq()) {
      if (rs->set_curr_seq()) {
        rs->extract_mers(mer_vvec, sh);
      }
    }
  }
  sort_columns();
  make_unique();
  /* update_nkmers(); */
  rs->compute_rho();
}

void FlatHT::display_info(std::ostream* output_stream, uint32_t r)
{
  vec<uint64_t> se_to_count;
  se_to_count.resize(crecord->get_nsubsets());
  for (uint64_t ix = 0; ix < nkmers; ++ix) {
    se_to_count[cmer_vv[ix].second]++;
  }
  crecord->display_info(output_stream, r, se_to_count);
}
