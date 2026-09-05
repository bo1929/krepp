#include "krepp.hpp"

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h, m); }

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

void TargetSketch::load_sketch()
{
  sketch->load_full_sketch();
  sketch->make_rho_partial();
}

void TargetIndex::read_lineages()
{
  qtree = std::make_shared<Tree>();
  std::ifstream lineage_stream(lineage_path);
  CHECK_STREAM_OR_EXIT(lineage_stream, (std::string("Error opening ") + lineage_path.string()));
  qtree->parse_lineages(lineage_stream);
  qtree->reset_traversal();
  index->get_tree()->map_to_qtree(qtree);
  lineage_stream.close();
}

void TargetIndex::ensure_backbone()
{
  if (!nwk_path.empty()) {
    qtree = std::make_shared<Tree>();
    std::ifstream tree_stream(nwk_path);
    CHECK_STREAM_OR_EXIT(tree_stream, (std::string("Error opening ") + nwk_path.string()));
    qtree->load(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to read the given backbone tree!");
    qtree->reset_traversal();
    index->get_tree()->map_to_qtree(qtree);
    tree_stream.close();
  } else if (index->check_wbackbone()) {
    qtree = index->get_tree();
  } else {
    error_exit("Given index lacks a tree and no backbone tree is provided...");
  }
}

void TargetIndex::load_index()
{
  const std::set<std::string> lall{"cmer", "crecord", "inc", "metadata", "tree", "reflist"};
  const std::set<std::string> lall_wbackbone{"cmer", "crecord", "inc", "metadata", "tree"};
  const std::set<std::string> lall_wobackbone{"cmer", "crecord", "inc", "metadata", "reflist"};
  node_phmap<std::string, std::set<std::string>> suffix_to_ltype;
  for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    std::string filename, ltype, mrcfg, fracv;
    size_t pos1, pos2;
    filename = entry.path().filename();
    pos1 = filename.find("-", 0);
    pos2 = filename.find("-", pos1 + 1);
    ltype = filename.substr(0, pos1);
    if (lall.find(ltype) != lall.end()) {
      if (entry.path().extension().empty()) {
        mrcfg = filename.substr(pos1, pos2 - pos1);
        fracv = filename.substr(pos2, filename.size() - pos2);
        suffix_to_ltype[mrcfg + fracv].insert(ltype);
      }
    } else {
    }
  }
  std::vector<std::string> suffixes;
  for (auto const& [suffix, ltypes] : suffix_to_ltype) {
    suffixes.push_back(suffix);
  }
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (uint32_t lix = 0; lix < suffixes.size(); ++lix) {
    const std::set<std::string>& ltype = suffix_to_ltype[suffixes[lix]];
    bool wobackbone = std::includes(ltype.begin(), ltype.end(), lall_wobackbone.begin(), lall_wobackbone.end());
    bool wbackbone = std::includes(ltype.begin(), ltype.end(), lall_wbackbone.begin(), lall_wbackbone.end());
    if (wbackbone) {
      index->load_partial_tree(suffixes[lix]);
      index->load_partial_index(suffixes[lix]);
    } else if (wobackbone) {
      index->generate_partial_tree(suffixes[lix]);
      index->load_partial_index(suffixes[lix]);
    } else {
      error_exit("There is a partial index with a missing file!");
    }
  }
  index->make_rho_partial();
}

void SketchSingle::create_sketch()
{
  rseq_sptr_t rs = std::make_shared<RSeq>(input, lshf, w, r, frac, sdust_t, sdust_w);
  sdynht_sptr_t sdynht = std::make_shared<SDynHT>();
  sdynht->fill_table(nrows, rs);
  sketch_sflatht = std::make_shared<SFlatHT>(sdynht);
  rho = rs->get_rho();
  std::cerr << "Total number of k-mers included in the sketch: " << sdynht->get_nkmers() << std::endl;
  std::cerr << "Subsampling rate (rho) is: " << rho << std::endl;
}

void SketchSingle::save_sketch()
{
  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  sketch_sflatht->save(sketch_stream);
  save_configuration(sketch_stream);
  sketch_stream.write(reinterpret_cast<char*>(&rho), sizeof(double));
  CHECK_STREAM_OR_EXIT(sketch_stream, "Failed to write the sketch!");
  sketch_stream.close();
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

void IndexMultiple::build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht)
{
  if (nd->check_leaf()) {
    sh_t sh = nd->get_sh();
    if (name_to_path.find(nd->get_name()) != name_to_path.end()) {
      rseq_sptr_t rs = std::make_shared<RSeq>(name_to_path[nd->get_name()], lshf, w, r, frac, sdust_t, sdust_w);
      dynht->fill_table(sh, rs);
      dynht->get_record()->insert_rho(nd->get_sh(), rs->get_rho());
#pragma omp critical
      {
        std::cerr << "\33[2K\r" << std::flush;
        std::cerr << "Leaf node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                  << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
      }
    } else {
#pragma omp critical
      {
        std::cerr << "\33[2K\r" << std::flush;
        std::cerr << "Genome skipped: " << nd->get_name() << "\r" << std::flush;
        build_count++;
      }
    }
  } else {
    assert(nd->get_nchildren() > 0);
    vec<dynht_sptr_t> children_dynht_v;
#if defined(_OPENMP) && _WOPENMP == 1
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
#endif
    children_dynht_v.reserve(nd->get_nchildren());
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
      children_dynht_v.emplace_back(std::make_shared<DynHT>(nrows, tree, dynht->get_record()));
#pragma omp task shared(dynht)
      {
        build_for_subtree(*std::next(nd->get_children(), i), children_dynht_v[i]);
#if defined(_OPENMP) && _WOPENMP == 1
        omp_set_lock(&parent_lock);
#endif
        dynht->union_table(children_dynht_v[i]);
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
      std::cerr << "\33[2K\r" << std::flush;
      std::cerr << "Internal node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
    }
  }
}

void IndexMultiple::index_sequences()
{
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
        std::cerr << "\33[2K\r" << std::flush;
        std::cerr << "Leaf node: " << lf->get_name() << "\tsize: " << ltab->get_nkmers() << "\tprogress: " << (++build_count)
                  << "/" << tree->get_nnodes() << "\r" << std::flush;
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
          std::cerr << "\33[2K\r" << std::flush;
          std::cerr << "Internal node: " << parent->get_name() << "\tsize: " << (acc ? acc->get_nkmers() : 0)
                    << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r" << std::flush;
        }
        done = parent;
      }
    }
  };
#pragma omp parallel for num_threads(nsl) schedule(static)
  for (uint32_t six = 0; six < static_cast<uint32_t>(slices.size()); ++six) {
    fill_slice(slices[six].first, slices[six].second);
  }
  dynht_sptr_t root_dynht = se_to_table[tree->get_root()->get_se()];
  assertm(root_dynht && root_dynht->get_nkmers() > 0, "No k-mers to index!");
  root_flatht = std::make_shared<FlatHT>(root_dynht);
}

void IndexMultiple::index_files()
{
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
      build_for_subtree(tree->get_root(), root_dynht);
    }
  }
  assertm(root_dynht->get_nkmers() > 0, "No k-mers to index!");
  root_flatht = std::make_shared<FlatHT>(root_dynht);
}

void QuerySketch::header_dreport(strstream& dreport_stream)
{
  dreport_stream << "# software: krepp\tversion: " VERSION "\tinvocation :" + invocation;
  dreport_stream << "\nSEQ_ID\tDIST\n";
}

void QueryIndex::header_dreport(strstream& dreport_stream)
{
  dreport_stream << "# software: krepp\tversion: " VERSION "\tinvocation :" + invocation;
  if (summarize) {
    dreport_stream << "\nREFERENCE_NAME\tWEIGHTED_COUNT\tSEQUENCE_ABUNDANCE\n";
  } else {
    dreport_stream << "\nSEQ_ID\tREFERENCE_NAME\tDIST\n";
  }
}

void QuerySketch::seek_sequences()
{
  strstream dreport_stream;
  header_dreport(dreport_stream);
  (*output_stream) << dreport_stream.rdbuf();
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      bool cont_reading = false;
      while ((cont_reading = qs->read_next_batch()) || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        auto sb = std::make_shared<SBatch>(sketch, qs, hdist_th);
#pragma omp task firstprivate(sb)
        {
          sb->seek_sequences(*output_stream);
        }
      }
#pragma omp taskwait
    }
  }
}

void QueryIndex::estimate_distances()
{
  parallel_flat_phmap<node_sptr_t, double> node_to_wcount = {};
  double twcount = 0;
  strstream dreport_stream;
  dreport_stream.precision(STRSTREAM_PRECISION);
  dreport_stream << std::fixed;
  header_dreport(dreport_stream);
  (*output_stream) << dreport_stream.rdbuf();
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      bool cont_reading = false;
      while ((cont_reading = qs->read_next_batch()) || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        auto ib = std::make_shared<IBatch>(index, qs, hdist_th, chisq_value, dist_max, tau, no_filter, multi, summarize);
#pragma omp task firstprivate(ib)
        {
          strstream batch_stream;
          ib->estimate_distances(batch_stream);
#pragma omp critical
          {
            if (summarize) {
              for (auto& [nd, wcount] : ib->get_summary()) {
                twcount += wcount;
                node_to_wcount[nd] += wcount;
              }
            } else {
              (*output_stream) << batch_stream.rdbuf();
            }
          }
        }
      }
#pragma omp taskwait
    }
  }
  if (summarize) {
    for (auto& [nd, wcount] : node_to_wcount) {
      dreport_stream << nd->get_name() << "\t" << wcount << "\t" << wcount / twcount << "\n";
    }
    (*output_stream) << dreport_stream.rdbuf();
  }
}

void QueryIndex::header_preport(strstream& dreport_stream)
{
  dreport_stream << "# software: krepp\tversion: " VERSION "\tinvocation :" + invocation;
  dreport_stream << "\n# ";
  qtree->stream_nwk_jplace(dreport_stream, qtree->get_root());
  if (summarize) {
    dreport_stream << "\nDISTAL_NODE\tEDGE_NUM\tWEIGHTED_COUNT\tSEQUENCE_ABUNDANCE\n";
  } else if (tabular) {
    dreport_stream << "\nSEQ_ID\tDISTAL_NODE\tEDGE_NUM\tLWR\tDIST\n";
  } else {
    (void)0; // TODO: Perhaps introduce an alternative placement report format.
  }
}

void QueryIndex::end_jplace(strstream& jplace_stream)
{
  jplace_stream << "],\n";
  jplace_stream << "\t\"metadata\" : {\n"
                << "\t\t\"software\" : \"krepp\",\n"
                << "\t\t\"version\" : \"" VERSION "\",\n"
                << "\t\t\"repository\" : \"https://github.com/bo1929/krepp\",\n"
                << "\t\t\"num_queries\" : \"" << total_qseq << "\",\n"
                << "\t\t\"invocation\" : \"";
  jplace_stream << invocation;
  jplace_stream << "\"\n\t},\n";
  jplace_stream << "\t\"tree\" : \"";
  qtree->stream_nwk_jplace(jplace_stream, qtree->get_root());
  jplace_stream << "\"\n}";
}

void QueryIndex::begin_jplace(strstream& jplace_stream)
{
  // Keep it compatible with jplace standard.
  jplace_stream
    << "{\n\t\"version\" : 3,\n\t"
       "\"fields\" : [\"edge_num\", \"pendant_length\", \"distal_length\", \"likelihood\", \"like_weight_ratio\", \"distance\"],\n\t\"placements\" : [\n";
}

void QueryIndex::place_sequences()
{
  parallel_flat_phmap<node_sptr_t, double> node_to_wcount = {};
  double twcount = 0;
  strstream preport_stream;
  preport_stream.precision(STRSTREAM_PRECISION);
  preport_stream << std::fixed;
  if (summarize || tabular) {
    header_preport(preport_stream);
    (*output_stream) << preport_stream.rdbuf();
  } else {
    begin_jplace(preport_stream);
    (*output_stream) << preport_stream.rdbuf();
  }
  bool has_previous = false;
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      bool cont_reading = false;
      while ((cont_reading = qs->read_next_batch()) || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        auto ib = std::make_shared<IBatch>(index, qs, hdist_th, chisq_value, dist_max, tau, no_filter, multi, summarize);
#pragma omp task firstprivate(ib)
        {
          strstream batch_stream;
          ib->place_sequences(batch_stream, tabular);
#pragma omp critical
          {
            if (summarize) {
              for (auto& [nd, wcount] : ib->get_summary()) {
                twcount += wcount;
                node_to_wcount[nd] += wcount;
              }
            } else {
              if (tabular) {
                (*output_stream) << batch_stream.rdbuf();
              } else {
                if (batch_stream.tellp() != std::streampos(0)) {
                  if (has_previous) {
                    (*output_stream) << ",\n";
                  }
                  (*output_stream) << batch_stream.rdbuf();
                  has_previous = true;
                }
              }
            }
          }
        }
      }
#pragma omp taskwait
    }
  }
  preport_stream.str("");
  preport_stream.clear();
  if (summarize) {
    for (auto& [nd, wcount] : node_to_wcount) {
      preport_stream << nd->get_name(true) << "\t" << nd->get_en() << "\t" << wcount << "\t" << wcount / twcount << "\n";
    }
    (*output_stream) << preport_stream.rdbuf();
  } else if (tabular) {
    (void)0;
  } else {
    end_jplace(preport_stream);
    (*output_stream) << preport_stream.rdbuf();
  }
}

void InfoIndex::display_info() { index->display_info(output_stream); }

InfoIndex::InfoIndex(CLI::App& sc)
{
  sc.add_option("-i,--index-dir", index_dir, "Directory <path> containing reference index.")
    ->required()
    ->check(CLI::ExistingDirectory);
  sc.callback([&]() { index = std::make_shared<Index>(index_dir); });
}

SketchSingle::SketchSingle(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-file", input, "Input FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file.")->required();
  sc.add_option("-k,--kmer-len", k, "Length of k-mers. [26]")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of minimizer window (w>=k). [k+6]");
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH. [k-16]");
  sc.add_option("-m,--modulo-lsh", m, "Modulo value to partition LSH space. [4]")->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m. [1]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m. [true]");
  sc.add_option("--sdust-t", sdust_t, "SDUST threshold (NCBI dustmasker: 20). [0]")->check(CLI::NonNegativeNumber);
  sc.add_option("--sdust-w", sdust_w, "SDUST window (NCBI dustmasker: 64). [0]")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
    }
    if (!(sc.count("-h") + sc.count("--num-positions"))) {
      h = std::max(k - 16, 9);
    }
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
  });
}

QuerySketch::QuerySketch(CLI::App& sc)
{
  sc.add_option("-q,--query", query, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query.")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path>. [stdout]");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match. [4]")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
    sketch = std::make_shared<Sketch>(sketch_path);
  });
}

IndexMultiple::IndexMultiple(CLI::App& sc)
{
  set_index_defaults();
  sc.add_option(
      "-i,--input-file",
      input,
      "TSV file <path> mapping reference IDs to (gzip compatible) paths/URLs, or a single (gzip compatible) FASTA/FASTQ.")
    ->required()
    ->check(CLI::ExistingFile);
  sc.add_option("-o,--index-dir", index_dir, "Directory <path> in which the index will be stored.")->required();
  sc.add_option("-t,--nwk-file", nwk_path, "Path to the Newick file for the guide tree (must be rooted).")
    ->check(CLI::ExistingFile);
  sc.add_option("-k,--kmer-len", k, "Length of k-mers. [29]")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of minimizer window (w>k). [k+6]");
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH. [k-16]");
  sc.add_option("-m,--modulo-lsh", m, "Mudulo value to partition LSH space. [4]")->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m. [1]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m. [true]");
  sc.add_option("--sdust-t", sdust_t, "SDUST threshold (NCBI dustmasker: 20). [0]")->check(CLI::NonNegativeNumber);
  sc.add_option("--sdust-w", sdust_w, "SDUST window (NCBI dustmasker: 64). [0]")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
    }
    if (!(sc.count("-h") + sc.count("--num-positions"))) {
      h = k - 16;
    }
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    std::filesystem::create_directory(index_dir);
    suffix = "-";
    suffix += "m" + std::to_string(m) + "r" + std::to_string(r);
    suffix += frac ? "-frac" : "-no_frac";
  });
}

void QueryIndex::init_sc_place(CLI::App& sc)
{
  sc.add_option(
      "-t,--nwk-file",
      nwk_path,
      "Path to the Newick file for the (rooted) placement tree (overrides if the index has a backbone tree). "
      "May be a jplace-style decorated tree, whose edge numbers ({N}) are used directly in the output (all nodes must be decorated).")
    ->check(CLI::ExistingFile);
  sc.add_option(
      "-l,--lineage-file",
      lineage_path,
      "Path to the Greengenes/GTDB style taxonomic lineage file, the first column has to match reference IDs present in the index (tolerates missing IDs).")
    ->excludes("--nwk-file")
    ->excludes("-t")
    ->check(CLI::ExistingFile);
  sc.add_option("--tau", tau, "Highest Hamming distance for placement threshold (increase to relax). [2]")
    ->check(CLI::NonNegativeNumber);
  multi = true;
  sc.add_flag("--multi,!--no-multi",
              multi,
              "Output all candidate placements satisfying the filters (not just the largest clade). [true]");
  filter = true;
  sc.add_flag("--filter,!--no-filter",
              filter,
              "Filter a placement when there is not enough k-mer matches below threshold tau. [true]");
  sc.add_flag("--tabular,!--no-tabular",
              tabular,
              "Output the per query sequence placements (taxonomic/phylogenetic) in a tab-separated format. [false]");
  sc.callback([&]() {
    no_filter = !filter;
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
    index = std::make_shared<Index>(index_dir);
    if (!validate_configuration_place()) {
      error_exit("Invalid configuration!");
    }
  });
}

void QueryIndex::init_sc_dist(CLI::App& sc)
{
  sc.add_option("--dist-max", dist_max, "Maximum distance to report for matching references.")->check(CLI::Range(1e-8, 0.33));
  multi = true;
  sc.add_flag(
    "--multi,!--no-multi", multi, "Output all distances satisfying the filters (not just the closest reference). [true]");
  filter = false;
  sc.add_flag(
    "--filter,!--no-filter",
    filter,
    "Filter a hit if its distance is too high compared to the best hit (based on the statistical significance). [false]");
  sc.callback([&]() {
    no_filter = !filter;
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
    index = std::make_shared<Index>(index_dir);
    if (!validate_configuration_dist()) {
      error_exit("Invalid configuration!");
    }
  });
}

QueryIndex::QueryIndex(CLI::App& sc)
{
  sc.add_option("-q,--query", query, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--index-dir", index_dir, "Directory <path> containing the reference index.")
    ->required()
    ->check(CLI::ExistingDirectory);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path>. [stdout]");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match. [4]")->check(CLI::NonNegativeNumber);
  sc.add_option("--chisq",
                chisq_value,
                "Chi-square value for statistical distinguishability test, default correspons to alpha=90%. [2.706]")
    ->check(CLI::PositiveNumber);
  sc.add_flag(
    "--summarize,!--no-summarize",
    summarize,
    "Summarize results into a table of read counts. If a read is mapped/placed to n references/edges, each gets 1/n. Overrides --no-multi and --no-filter. [false]");
  /* sc.add_option("--leave-out-ref", leave_out_ref, "Reference ID to exclude, useful for testing."); */
}

int main(int argc, char** argv)
{
  PRINT_VERSION
  std::ios::sync_with_stdio(false);
  CLI::App app{"krepp: a tool for k-mer-based search, distance estimation & phylogenetic placement."};
  app.set_help_flag("--help");
  app.fallthrough();

  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report.");
  app.require_subcommand();
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness. [0]");
  app.callback([&]() {
    if (app.count("--seed")) {
      gen.seed(seed);
    }
  });
  app.add_option("--num-threads", num_threads, "Number of threads to use in OpenMP-based parallelism. [1]");

  auto& sc_index = *app.add_subcommand("index", "Build an index from k-mers of reference genomes.");
  auto& sc_place = *app.add_subcommand("place", "Place queries on a tree with respect to an index.");
  auto& sc_dist = *app.add_subcommand("dist", "Estimate distances of queries to genomes in an index.");
  auto& sc_inspect = *app.add_subcommand("inspect", "Display statistics and information for a given index.");
  auto& sc_sketch = *app.add_subcommand("sketch", "Create a sketch from k-mers in a single FASTA/FASTQ file.");
  auto& sc_seek = *app.add_subcommand("seek", "Seek query sequences in a sketch and estimate distances.");

  IndexMultiple krepp_index(sc_index);
  QueryIndex krepp_place(sc_place);
  krepp_place.init_sc_place(sc_place);
  QueryIndex krepp_dist(sc_dist);
  krepp_dist.init_sc_dist(sc_dist);
  InfoIndex krepp_inspect(sc_inspect);
  SketchSingle krepp_sketch(sc_sketch);
  QuerySketch krepp_seek(sc_seek);

  CLI11_PARSE(app, argc, argv);
  for (int i = 0; i < argc; ++i) {
    invocation += std::string(argv[i]) + " ";
  }
  invocation.pop_back();

  auto tstart = std::chrono::system_clock::now();
  std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  std::cerr << "Invocation: " << invocation << "\n";
  std::cerr << std::ctime(&tstart_f);

  if (sc_index.parsed()) {
    std::cerr << "Reading the tree and initializing the index..." << std::endl;
    krepp_index.set_nrows();
    krepp_index.set_lshf();
    krepp_index.read_input_file();
    krepp_index.obtain_build_tree();
    std::cerr << "Building the index..." << std::endl;
    krepp_index.build_index();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    std::cerr << std::endl << "Finished indexing, elapsed: " << es_b.count() << " sec" << std::endl;
    krepp_index.save_index();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done converting & saving, elapsed: " << es_s.count() << " sec" << std::endl;
  }

  if (sc_place.parsed()) {
    std::cerr << "Loading the index and the backbone tree..." << std::endl;
    krepp_place.load_index();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    if (sc_place.count("--lineage-file") + sc_place.count("-l")) {
      krepp_place.read_lineages();
      std::cerr << "Placing given sequences on the taxonomic lineage..." << std::endl;
    } else {
      krepp_place.ensure_backbone();
      std::cerr << "Placing given sequences on the backbone tree..." << std::endl;
    }
    krepp_place.place_sequences();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done placing queries, elapsed: " << es_s.count() << " sec" << std::endl;
    std::cerr << "Total number of sequences queried: " << krepp_place.get_total_qseq() << std::endl;
  }

  if (sc_dist.parsed()) {
    std::cerr << "Loading the index and initializing..." << std::endl;
    krepp_dist.load_index();
    std::cerr << "Estimating distances between given sequences and references..." << std::endl;
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    krepp_dist.estimate_distances();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done estimating distances, elapsed: " << es_s.count() << " sec" << std::endl;
    std::cerr << "Total number of sequences queried: " << krepp_dist.get_total_qseq() << std::endl;
  }

  if (sc_inspect.parsed()) {
    std::cerr << "Inspecting the index..." << std::endl;
    krepp_inspect.load_index();
    krepp_inspect.display_info();
    std::cerr << "Done reporting the index information..." << std::endl;
  }

  if (sc_sketch.parsed()) {
    std::cerr << "Initializing the sketch..." << std::endl;
    krepp_sketch.set_nrows();
    krepp_sketch.set_lshf();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    krepp_sketch.create_sketch();
    krepp_sketch.save_sketch();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done sketching & saving, elapsed: " << es_s.count() << " sec" << std::endl;
  }

  if (sc_seek.parsed()) {
    std::cerr << "Loading the sketch..." << std::endl;
    krepp_seek.load_sketch();
    std::cerr << "Seeking query sequences in the sktech..." << std::endl;
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    krepp_seek.seek_sequences();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done seeking sequences, elapsed: " << es_s.count() << " sec" << std::endl;
    std::cerr << "Total number of sequences queried: " << krepp_seek.get_total_qseq() << std::endl;
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
