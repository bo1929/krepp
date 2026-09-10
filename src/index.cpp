#include "index.hpp"
#include <atomic>

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
  bool compatible = false;
#pragma omp critical
  {
    compatible = curr_tree->check_compatible(tree);
    if (compatible) tree = !tree ? curr_tree : tree;
  }
  if (!compatible) error_exit("Partial libraries are based on different references.");
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
  bool compatible = false;
#pragma omp critical
  {
    compatible = curr_tree->check_compatible(tree);
    if (compatible) tree = !tree ? curr_tree : tree;
  }
  if (!compatible) error_exit("Partial libraries are based on different trees!");
}

void Index::load_partial_index(std::string suffix)
{ // TODO: Split for each file (e.g, metadata, crecord etc.)
  std::filesystem::path metadata_path = index_dir / ("metadata" + suffix);
  std::ifstream metadata_stream(metadata_path, std::ifstream::binary);
  if (!metadata_stream.is_open()) {
    error_exit(std::string("Failed to open ") + metadata_path.string());
  }
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

  lshf_sptr_t curr_lshf = std::make_shared<LSHF>(m_curr, ppos_v, npos_v, r, frac);
  bool compatible = false;
#pragma omp critical
  {
    compatible = curr_lshf->check_compatible(lshf);
    if (compatible) {
      lshf = !lshf ? curr_lshf : lshf;
      k = k_curr;
      h = h_curr;
      m = m_curr;
      nrows = pow(2, 2 * h);
    }
  }
  if (!compatible) error_exit("Partial libraries have incompatible hash configurations!");

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

  std::string info_str;
  std::filesystem::path info_path = index_dir / ("metadata" + suffix + ".txt");
  if (std::filesystem::exists(info_path)) {
    std::ifstream info_stream(info_path);
    info_str.assign((std::istreambuf_iterator<char>(info_stream)), std::istreambuf_iterator<char>());
    CHECK_STREAM_OR_EXIT(info_stream, "Failed to metadata info_str text of a partial index!");
  } else {
    info_str += "krepp version: ?\n";
    info_str += "date: ?\n";
    info_str += "seed: ?\n";
    info_str += "k: " + std::to_string(static_cast<uint32_t>(k)) + "\n";
    info_str += "w: " + std::to_string(static_cast<uint32_t>(w)) + "\n";
    info_str += "h: " + std::to_string(static_cast<uint32_t>(h)) + "\n";
    info_str += "m: " + std::to_string(m) + "\n";
    info_str += frac ? "frac: true\n" : "frac: false\n";
    info_str += "ppos_v: " + vec_to_str(curr_lshf->get_ppos()) + "\n";
    info_str += "npos_v: " + vec_to_str(curr_lshf->get_npos()) + "\n";
    info_str += "nrows: " + std::to_string(nrows) + "\n";
    info_str += "total_num_kmers: " + std::to_string(curr_flatht->get_nkmers()) + "\n";
    info_str += "sdust-t: ?\n";
    info_str += "sdust-w: ?\n";
  }

#pragma omp critical
  {
    if (frac) {
      for (uint32_t ix = 0; ix <= r; ++ix) {
        r_to_flatht[ix] = curr_flatht;
        r_to_numerator[ix] = r + 1;
        r_to_info[ix] = info_str;
      }
    } else {
      r_to_flatht[r] = curr_flatht;
      r_to_numerator[r] = 1;
      r_to_info[r] = info_str;
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

void Index::display_info(std::ostream* output_stream)
{
  if (wbackbone) {
    strstream newick_stream;
    tree->stream_nwk_basic(newick_stream, tree->get_root());
    (*output_stream) << "Backbone tree: " << newick_stream.rdbuf() << "\n";
  } else {
    (*output_stream) << "Backbone tree: NA\n";
  }
  for (auto const& [key, val] : r_to_info) {
    (*output_stream) << "======= Partial index: " << key << " =======\n";
    (*output_stream) << r_to_info[key];
    r_to_flatht[key]->display_info(output_stream, key);
  }
}

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

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h, m, r, frac); }

void BaseLSH::set_nrows()
{
  uint32_t hash_size = pow(2, 2 * h);
  uint32_t full_residue = hash_size % m;
  if (frac) {
    nrows = (hash_size / m) * (r + 1);
    nrows = full_residue > r ? nrows + (r + 1) : nrows + full_residue;
  } else {
    nrows = (hash_size / m);
    nrows = full_residue > r ? nrows + 1 : nrows;
  }
}

void BaseLSH::save_configuration(std::ofstream& cfg_stream)
{
  cfg_stream.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(&frac), sizeof(bool));
  cfg_stream.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(lshf->ppos_data()), (h) * sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
}

IndexMultiple::IndexMultiple(const IndexConfig& config)
{
  input = config.input;
  index_dir = config.index_dir;
  nwk_path = config.nwk_path;
  k = config.k;
  w = config.w.value_or(k + 6);
  h = config.h.value_or(k - 16);
  m = config.m;
  r = config.r;
  frac = config.frac;
  sdust_t = config.sdust_t;
  sdust_w = config.sdust_w;
  if (!validate_configuration()) {
    error_exit("Invalid configuration!");
  }
  nrows = pow(2, 2 * h - 1);
  std::filesystem::create_directory(index_dir);
  suffix = "-";
  suffix += "m" + std::to_string(m) + "r" + std::to_string(r);
  suffix += frac ? "-frac" : "-no_frac";
}

void IndexMultiple::obtain_build_tree()
{
  tree = std::make_shared<Tree>();
  if (per_sequence) {
    if (!nwk_path.empty()) {
      error_exit("A guide tree (-t) is incompatible with a per sequence indexing.");
    }
    std::cerr << "No guide tree for per sequence indexing." << std::endl;
    tree->generate_tree(names_v);
  } else if (nwk_path.empty()) {
    std::cerr << "No tree has given as a guide, the color index could be suboptimal." << std::endl;
    tree->generate_tree(names_v);
  } else {
    std::ifstream tree_stream(nwk_path);
    CHECK_STREAM_OR_EXIT(tree_stream, (std::string("Error opening ") + nwk_path.string()));
    tree->load(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to read the backbone tree of the index!");
    tree_stream.close();
  }
  tree->reset_traversal();
}

void IndexMultiple::read_input_file()
{
  gzFile gfile = gzopen(input.c_str(), "rb");
  if (gfile != nullptr) {
    int c;
    while ((c = gzgetc(gfile)) != -1 && (c == '\n' || c == '\r' || c == ' ' || c == '\t')) {
    }
    gzrewind(gfile);
    kseq_t* kseq = nullptr;
    bool is_fastx = false;
    if (c == '>' || c == '@') {
      kseq = kseq_init(gfile);
      int osk = kseq_read(kseq);
      if (osk < -1) {
        error_exit("Error reading the input (truncated FASTQ record?).");
      }
      is_fastx = osk >= 0 && kseq->name.l > 0 && (c == '>' || (kseq->seq.l > 0 && kseq->qual.l == kseq->seq.l));
    }
    if (is_fastx) {
      gzrewind(gfile);
      kseq_rewind(kseq);
      per_sequence = true;
      flat_phmap<std::string, bool> seen_names;
      uint64_t r_offset = 0, n_offset = 0;
      int kret;
      while ((kret = kseq_read(kseq)) >= 0) {
        std::string name(kseq->name.s);
        if (name.empty()) {
          error_exit("Empty reference ID in the input!");
        }
        if (seen_names.contains(name)) {
          error_exit("Duplicate reference ID \"" + name + "\" in the input!");
        }
        seen_names[name] = true;
        fastx_names.push_back(name);
        fastx_offsets.push_back(r_offset);
        if (kseq->seq.l < w) {
          std::cerr << "[WARNING] Skipping \"" << name << "\" it is shorter than the minimizer window." << std::endl;
        } else {
          names_v.push_back(name);
        }
        n_offset = gztell(kseq->f->f) - (kseq->f->end - kseq->f->begin) - (kseq->last_char ? 1 : 0);
        r_offset = n_offset;
      }
      if (kret < -1) {
        error_exit("Error reading the input (truncated record?).");
      }
      fastx_offsets.push_back(n_offset);
      kseq_destroy(kseq);
      gzclose(gfile);
      if (names_v.empty()) {
        error_exit("No sequences of longer the minimizer length found in the input!");
      }
      return;
    }
    kseq_destroy(kseq);
    gzclose(gfile);
  }
  std::ifstream input_stream(input);
  CHECK_STREAM_OR_EXIT(input_stream, (std::string("Error opening ") + input.string()));
  std::string line;
  while (std::getline(input_stream, line)) {
    std::istringstream iss(line);
    std::string input, name;
    if (!(std::getline(iss, name, '\t') && std::getline(iss, input, '\t'))) {
      error_exit("Failed to read the reference name to path/URL mapping!");
    }
    if (name_to_path.contains(name)) {
      error_exit("Duplicate reference ID \"" + name + "\" in the input map file!");
    }
    name_to_path[name] = input;
    names_v.push_back(name);
  }
  input_stream.close();
}

void IndexMultiple::build_index()
{
  build_count = 0;
  if (per_sequence) {
    index_sequences();
  } else {
    index_files();
  }
}

void IndexMultiple::save_info(std::ofstream& info_stream)
{
  info_stream << "krepp version: " << VERSION << "\n";
  std::time_t t = std::time(nullptr);
  info_stream << "date: " << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "\n";
  info_stream << "seed: " << seed << "\n";
  info_stream << "k: " << static_cast<uint32_t>(k) << "\n";
  info_stream << "w: " << static_cast<uint32_t>(w) << "\n";
  info_stream << "h: " << static_cast<uint32_t>(h) << "\n";
  info_stream << "m: " << m << "\n";
  info_stream << "frac: " << (frac ? "true" : "false") << "\n";
  info_stream << "ppos_v: " << vec_to_str(lshf->get_ppos()) << "\n";
  info_stream << "npos_v: " << vec_to_str(lshf->get_npos()) << "\n";
  info_stream << "nrows: " << nrows << "\n";
  info_stream << "total_num_kmers: " << root_flatht->get_nkmers() << "\n";
  info_stream << "sdust-t: " << sdust_t << "\n";
  info_stream << "sdust-w: " << sdust_w << "\n";
}

void IndexMultiple::save_index()
{
  std::ofstream mer_stream(index_dir / ("cmer" + suffix), std::ofstream::binary);
  std::ofstream inc_stream(index_dir / ("inc" + suffix), std::ofstream::binary);
  root_flatht->save(mer_stream, inc_stream);
  CHECK_STREAM_OR_EXIT(mer_stream, "Failed to write the k-mer array of the index!");
  CHECK_STREAM_OR_EXIT(inc_stream, "Failed to read the offset array of a partial index!");
  inc_stream.close();
  mer_stream.close();

  std::ofstream crecord_stream(index_dir / ("crecord" + suffix), std::ofstream::binary);
  root_flatht->get_crecord()->save(crecord_stream);
  CHECK_STREAM_OR_EXIT(crecord_stream, "Failed to write the color array of the index!");
  crecord_stream.close();

  std::ofstream reflist_stream(index_dir / ("reflist" + suffix));
  std::ostream_iterator<std::string> reflist_iterator(reflist_stream, "\n");
  std::copy(std::begin(names_v), std::end(names_v), reflist_iterator);
  CHECK_STREAM_OR_EXIT(reflist_stream, "Failed to write the reference list of the index!");
  reflist_stream.close();
  if (!nwk_path.empty()) {
    std::ofstream tree_stream(index_dir / ("tree" + suffix));
    root_flatht->get_tree()->save(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to write the backbone tree of the index!");
    tree_stream.close();
  } else {
    std::cerr << "Skipped saving a backbone for the index!" << std::endl;
  }

  std::filesystem::path metadata_path = index_dir / ("metadata" + suffix);
  std::ofstream metadata_stream(metadata_path, std::ofstream::binary);
  save_configuration(metadata_stream);
  CHECK_STREAM_OR_EXIT(metadata_stream, "Failed to write the metadata of the index!");
  metadata_stream.close();

  // Save human-readable info
  std::ofstream info_stream(index_dir / ("metadata" + suffix + ".txt"));
  save_info(info_stream);
  CHECK_STREAM_OR_EXIT(info_stream, "Failed to write the text metadata of the index!");
  info_stream.close();
}

void IndexMultiple::build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht, ErrorRelay& relay)
{
  if (nd->check_leaf()) {
    sh_t sh = nd->get_sh();
    if (name_to_path.find(nd->get_name()) != name_to_path.end()) {
      rseq_sptr_t rs = std::make_shared<RSeq>(name_to_path[nd->get_name()], lshf, w, r, frac, sdust_t, sdust_w);
      dynht->fill_table(sh, rs);
      dynht->get_record()->insert_rho(nd->get_sh(), rs->get_rho());
#pragma omp critical
      {
        relay.guard([&] {
          std::cerr << "\33[2K\r" << std::flush;
          std::cerr << "Leaf node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                    << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
        });
      }
    } else {
#pragma omp critical
      {
        relay.guard([&] {
          std::cerr << "\33[2K\r" << std::flush;
          std::cerr << "Genome skipped: " << nd->get_name() << "\r" << std::flush;
          build_count++;
        });
      }
    }
  } else {
    assert(nd->get_nchildren() > 0);
    vec<dynht_sptr_t> children_dynht_v;
    children_dynht_v.reserve(nd->get_nchildren());
#if defined(_OPENMP) && _WOPENMP == 1
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
#endif
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
      bool prepared = false;
      relay.guard([&] {
        children_dynht_v.emplace_back(std::make_shared<DynHT>(nrows, tree, dynht->get_record()));
        prepared = true;
      });
      if (!prepared) break;
      node_sptr_t child = *std::next(nd->get_children(), i);
      dynht_sptr_t child_dynht = children_dynht_v[i];
#if defined(_OPENMP) && _WOPENMP == 1
  #pragma omp task shared(dynht, parent_lock, relay) firstprivate(child, child_dynht)
#endif
      {
        relay.guard([&] { build_for_subtree(child, child_dynht, relay); });
#if defined(_OPENMP) && _WOPENMP == 1
        omp_set_lock(&parent_lock);
#endif
        /* guard() is noexcept, so the lock is released on every path. */
        relay.guard([&] { dynht->union_table(child_dynht); });
#if defined(_OPENMP) && _WOPENMP == 1
        omp_unset_lock(&parent_lock);
#endif
      }
    }
#pragma omp taskwait
#if defined(_OPENMP) && _WOPENMP == 1
    omp_destroy_lock(&parent_lock);
#endif
#pragma omp critical
    {
      relay.guard([&] {
        std::cerr << "\33[2K\r" << std::flush;
        std::cerr << "Internal node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                  << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
      });
    }
  }
}

void IndexMultiple::index_sequences()
{
  ErrorRelay relay;
  record_sptr_t record = std::make_shared<Record>(tree);
  const uint32_t nrec = static_cast<uint32_t>(fastx_names.size());
  vec<std::atomic<int32_t>> pending_children(tree->get_nnodes() + 1);
  flat_phmap<std::string, node_sptr_t> name_to_leaf;
  tree->reset_traversal();
  node_sptr_t nd;
  while ((nd = tree->next_post_order())) {
    if (nd->check_leaf()) {
      name_to_leaf[nd->get_name()] = nd;
    }
    pending_children[nd->get_se()] = nd->get_nchildren();
  }
  vec<dynht_sptr_t> se_to_table(tree->get_nnodes() + 1);
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  const uint32_t nsl = std::min<uint32_t>(std::max<uint32_t>(num_threads, 1), nrec);
  vec<std::pair<uint32_t, uint32_t>> slices;
  {
    uint32_t b = 0;
    for (uint32_t s = 1; s <= nsl; ++s) {
      uint64_t target = (s == nsl) ? std::numeric_limits<uint64_t>::max()
                                   : (fastx_offsets[nrec] - fastx_offsets[0]) * s / nsl + fastx_offsets[0];
      uint32_t e = b;
      while (e < nrec && fastx_offsets[e] < target)
        e++;
      if (e > b) slices.emplace_back(b, e);
      b = e;
    }
  }
  auto fill_slice = [&](uint32_t b, uint32_t e) {
    rseq_sptr_t rs = std::make_shared<RSeq>(input.string(), lshf, w, r, frac, sdust_t, sdust_w, fastx_offsets[b]);
    for (uint32_t i = b; i < e; ++i) {
      if (!rs->read_next_seq()) {
        error_exit("FASTX record missing during the build (offset desync at record " + fastx_names[i] + ").");
      }
      bool usable = rs->set_curr_seq();
      if (fastx_names[i] != rs->get_name()) {
        error_exit("FASTX record mismatch at offset " + std::to_string(fastx_offsets[i]) + "; expected \"" + fastx_names[i] +
                   "\" but read \"" + rs->get_name() + "\".");
      }
      auto lit = name_to_leaf.find(fastx_names[i]);
      if (lit == name_to_leaf.end() || !usable) {
        continue; // Short record skipped from the index; it has been consumed.
      }
      node_sptr_t lf = lit->second;
      dynht_sptr_t ltab = std::make_shared<DynHT>(nrows, tree, record);
      ltab->fill_table(lf->get_sh(), rs, true);
      record->insert_rho(lf->get_sh(), rs->get_rho());
      se_to_table[lf->get_se()] = ltab;
#pragma omp critical
      {
        relay.guard([&] {
          std::cerr << "\33[2K\r" << std::flush;
          std::cerr << "Leaf node: " << lf->get_name() << "\tsize: " << ltab->get_nkmers()
                    << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
        });
      }
      node_sptr_t done = lf;
      while (true) {
        node_sptr_t parent = done->get_parent();
        if (!parent) break;
        if (pending_children[parent->get_se()].fetch_sub(1, std::memory_order_acq_rel) != 1) break;
        dynht_sptr_t acc;
        const tuint_t nch = parent->get_nchildren();
        for (tuint_t cix = 0; cix < nch; ++cix) {
          node_sptr_t child = *std::next(parent->get_children(), cix);
          dynht_sptr_t ctab = se_to_table[child->get_se()];
          se_to_table[child->get_se()] = nullptr;
          if (!ctab) continue;
          if (!acc) {
            acc = ctab;
            continue;
          }
          acc->union_table(ctab);
        }
        se_to_table[parent->get_se()] = acc;
#pragma omp critical
        {
          relay.guard([&] {
            std::cerr << "\33[2K\r" << std::flush;
            std::cerr << "Internal node: " << parent->get_name() << "\tsize: " << (acc ? acc->get_nkmers() : 0)
                      << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
          });
        }
        done = parent;
      }
    }
  };
#pragma omp parallel for num_threads(nsl) schedule(static)
  for (uint32_t six = 0; six < static_cast<uint32_t>(slices.size()); ++six) {
    const std::pair<uint32_t, uint32_t> slice = slices[six];
    relay.guard([&] { fill_slice(slice.first, slice.second); });
  }
  relay.rethrow_error();
  dynht_sptr_t root_dynht = se_to_table[tree->get_root()->get_se()];
  assertm(root_dynht && root_dynht->get_nkmers() > 0, "No k-mers to index!");
  root_flatht = std::make_shared<FlatHT>(root_dynht);
}

void IndexMultiple::index_files()
{
  ErrorRelay relay;
  record_sptr_t record = std::make_shared<Record>(tree);
  dynht_sptr_t root_dynht = std::make_shared<DynHT>(nrows, tree, record);
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
  #if _OPENMP >= 202011
  omp_set_max_active_levels(2);
  #else
  omp_set_nested(1);
  #endif
#endif
#pragma omp parallel
  {
#pragma omp single
    {
      relay.guard([&] { build_for_subtree(tree->get_root(), root_dynht, relay); });
    }
  }
  relay.rethrow_error();
  assertm(root_dynht->get_nkmers() > 0, "No k-mers to index!");
  root_flatht = std::make_shared<FlatHT>(root_dynht);
}
