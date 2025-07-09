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
void TargetIndex::ensure_wbackbone()
{
  if (!index->check_wbackbone()) {
    error_exit("Given index lacks a backbone tree required for this subcommand...");
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
    if (lall.find(ltype) != lall.end() && entry.path().extension().empty()) {
      mrcfg = filename.substr(pos1, pos2 - pos1);
      fracv = filename.substr(pos2, filename.size() - pos2);
      suffix_to_ltype[mrcfg + fracv].insert(ltype);
    }
  }
  std::vector<std::string> suffixes;
  for (auto const& [suffix, ltypes] : suffix_to_ltype) {
    suffixes.push_back(suffix);
  }
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (uint32_t lix = 0; lix < suffixes.size(); ++lix) {
    if (suffix_to_ltype[suffixes[lix]] == lall_wbackbone) {
      index->load_partial_tree(suffixes[lix]);
      index->load_partial_index(suffixes[lix]);
    } else if (suffix_to_ltype[suffixes[lix]] == lall_wobackbone) {
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
  if (nwk_path.empty()) {
    std::cerr << "No tree has given as a guide, the color index could be suboptimal." << std::endl;
    tree->generate_tree(names_v);
  } else {
    std::ifstream tree_stream(nwk_path);
    CHECK_STREAM_OR_EXIT(tree_stream, (std::string("Error opening ") + nwk_path.string()));
    tree->load(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to read the backbone tree of the index!");
  }
  tree->reset_traversal();
}

void IndexMultiple::read_input_file()
{
  std::ifstream input_stream(input);
  CHECK_STREAM_OR_EXIT(input_stream, (std::string("Error opening ") + input.string()));
  std::string line;
  while (std::getline(input_stream, line)) {
    std::istringstream iss(line);
    std::string input, name;
    if (!(std::getline(iss, name, '\t') && std::getline(iss, input, '\t'))) {
      error_exit("Failed to read the reference name to path/URL mapping!");
    }
    name_to_path[name] = input;
    names_v.push_back(name);
  }
  input_stream.close();
}

void IndexMultiple::build_index()
{
  record_sptr_t record = std::make_shared<Record>(tree);
  dynht_sptr_t root_dynht = std::make_shared<DynHT>(nrows, tree, record);
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
  omp_set_nested(1);
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

static std::string vec_to_str(const std::vector<uint8_t>& v)
{
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i > 0) oss << ", ";
    oss << static_cast<int>(v[i]);
  }
  oss << "]";
  return oss.str();
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
  info_stream << "SDUST-T: " << sdust_t << "\n";
  info_stream << "SDUST-W: " << sdust_w << "\n";
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

  if (nwk_path.empty()) {
    std::cerr << "Skipped saving a backbone for the index!" << std::endl;
    std::ofstream reflist_stream(index_dir / ("reflist" + suffix));
    std::ostream_iterator<std::string> reflist_iterator(reflist_stream, "\n");
    std::copy(std::begin(names_v), std::end(names_v), reflist_iterator);
    CHECK_STREAM_OR_EXIT(reflist_stream, "Failed to write the reference list of the index!");
    reflist_stream.close();
  } else {
    std::ofstream tree_stream(index_dir / ("tree" + suffix));
    root_flatht->get_tree()->save(tree_stream);
    CHECK_STREAM_OR_EXIT(tree_stream, "Failed to write the backbone tree of the index!");
    tree_stream.close();
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
      rseq_sptr_t rs =
        std::make_shared<RSeq>(name_to_path[nd->get_name()], lshf, w, r, frac, sdust_t, sdust_w);
      dynht->fill_table(sh, rs);
      dynht->get_record()->insert_rho(nd->get_sh(), rs->get_rho());
#pragma omp critical
      {
        std::cerr << "\33[2K\r" << std::flush;
        std::cerr << "Leaf node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                  << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r"
                  << std::flush;
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
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
      children_dynht_v.emplace_back(std::make_shared<DynHT>(nrows, tree, dynht->get_record()));
#pragma omp task untied shared(dynht)
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
                << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r"
                << std::flush;
    }
  }
}

void QuerySketch::header_dreport(strstream& dreport_stream)
{
  dreport_stream << "#software: krepp\t#version: " VERSION "\t#invocation :" + invocation;
  dreport_stream << "\nSEQ_ID\tDIST\n";
}

void QueryIndex::header_dreport(strstream& dreport_stream)
{
  dreport_stream << "#software: krepp\t#version: " VERSION "\t#invocation :" + invocation;
  dreport_stream << "\nSEQ_ID\tREFERENCE_NAME\tDIST\n";
}

void QuerySketch::seek_sequences()
{
  strstream dreport_stream;
  header_dreport(dreport_stream);
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      while (qs->read_next_batch() || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        SBatch sb(sketch, qs, hdist_th);
#pragma omp task untied
        {
          sb.seek_sequences(*output_stream);
        }
      }
#pragma omp taskwait
    }
  }
}

void QueryIndex::estimate_distances()
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
      while (qs->read_next_batch() || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        IBatch ib(index, qs, hdist_th, chisq_value, dist_max, tau, no_filter, multi);
#pragma omp task untied
        {
          ib.estimate_distances(*output_stream);
        }
      }
#pragma omp taskwait
    }
  }
}

void QueryIndex::end_jplace(strstream& jplace_stream)
{
  jplace_stream << "\t\t\t{\"n\" : [\"NaN\"], \"p\" : [["
                << index->get_tree()->get_root()->get_se() - 1 << ", 0, 0, 0, 0, 0]]}\n";
  jplace_stream << "\t],\n";
  jplace_stream << "\t\"metadata\" : {\n"
                << "\t\t\"software\" : \"krepp\",\n"
                << "\t\t\"version\" : \"" VERSION "\",\n"
                << "\t\t\"repository\" : \"https://github.com/bo1929/krepp\",\n"
                << "\t\t\"num_queries\" : \"" << total_qseq << "\",\n"
                << "\t\t\"invocation\" : \"";
  jplace_stream << invocation;
  jplace_stream << "\"\n\t},\n";
  jplace_stream << "\t\"tree\" : \"";
  index->get_tree()->stream_newick_str(jplace_stream, index->get_tree()->get_root());
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
  strstream jplace_stream;
  begin_jplace(jplace_stream);
  (*output_stream) << jplace_stream.rdbuf();
#if defined(_OPENMP) && _WOPENMP == 1
  omp_set_num_threads(num_threads);
#endif
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      while (qs->read_next_batch() || !qs->is_batch_finished()) {
        total_qseq += qs->get_cbatch_size();
        IBatch ib(index, qs, hdist_th, chisq_value, dist_max, tau, no_filter, multi);
#pragma omp task untied
        {
          ib.place_sequences(*output_stream);
        }
      }
#pragma omp taskwait
    }
  }
  // output_stream->seekp(-1, output_stream->cur);
  // output_stream->seekp(-1, std::ios_base::end);
  jplace_stream.str("");
  end_jplace(jplace_stream);
  (*output_stream) << jplace_stream.rdbuf();
}

void InfoIndex::display_info() { index->display_info(); }

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
  sc.add_option(
      "-i,--input-file", input, "Input FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file.")
    ->required();
  sc.add_option("-k,--kmer-len", k, "Length of k-mers [26].")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of minimizer window (w>=k) [k+6].");
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16].");
  sc.add_option("-m,--modulo-lsh", m, "Mudulo value to partition LSH space [5].")
    ->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m [1].")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m [true].");
  sc.add_option("--sdust-t", sdust_t, "SDUST threshold [20].")->check(CLI::NonNegativeNumber);
  sc.add_option("--sdust-w", sdust_w, "SDUST window [64].")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
      h = k - 16;
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
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query.")
    ->required()
    ->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout].");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4].")
    ->check(CLI::NonNegativeNumber);
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
  sc.add_option("-i,--input-file",
                input,
                "TSV file <path> mapping reference IDs to (gzip compatible) paths/URLs.")
    ->required()
    ->check(CLI::ExistingFile);
  sc.add_option("-o,--index-dir", index_dir, "Directory <path> in which the index will be stored.")
    ->required();
  sc.add_option(
      "-t,--nwk-file", nwk_path, "Path to the Newick file for the reference tree (must be rooted).")
    ->check(CLI::ExistingFile);
  sc.add_option("-k,--kmer-len", k, "Length of k-mers [29].")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of minimizer window (w>k) [k+6].");
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16].");
  sc.add_option("-m,--modulo-lsh", m, "Mudulo value to partition LSH space [5].")
    ->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m [1].")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m [true].");
  sc.add_option("--sdust-t", sdust_t, "SDUST threshold [20].")->check(CLI::NonNegativeNumber);
  sc.add_option("--sdust-w", sdust_w, "SDUST window [64].")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
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
      "--tau", tau, "Highest Hamming distance for placement threshold (increase to relax) [2].")
    ->check(CLI::NonNegativeNumber);
  multi = false;
  sc.add_flag(
    "--multi",
    multi,
    "Output all candidate placements satisfying the filters (not just the largest clade). [false]");
}

void QueryIndex::init_sc_dist(CLI::App& sc)
{
  sc.add_option(
      "--dist-max",
      dist_max,
      "Maximum distance to report for matching references, the output may become too large if high [0.2].")
    ->check(CLI::Range(1e-8, 0.33));
  multi = true;
  sc.add_flag(
    "--multi", multi, "Output all distances satisfying the filter (not just the best one). [true]");
}

QueryIndex::QueryIndex(CLI::App& sc)
{
  sc.add_option("-q,--query", query, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--index-dir", index_dir, "Directory <path> containing the reference index.")
    ->required()
    ->check(CLI::ExistingDirectory);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout].");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4].")
    ->check(CLI::NonNegativeNumber);
  sc.add_option(
      "--chisq",
      chisq_value,
      "Chi-square value for statistical distinguishability test, default correspons to alpha=90% [2.706].")
    ->check(CLI::PositiveNumber);
  sc.add_flag(
    "--no-filter,!--filter",
    no_filter,
    "Report everything; regardless of the statistical significance, match count, or maximum distance.");
  /* sc.add_option("--leave-out-ref", leave_out_ref, "Reference ID to exclude, useful for testing."); */
  sc.callback([&]() {
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
    index = std::make_shared<Index>(index_dir);
  });
}

int main(int argc, char** argv)
{
  PRINT_VERSION
  std::ios::sync_with_stdio(false);
  CLI::App app{
    "krepp: a tool for k-mer-based search, distance estimation & phylogenetic placement."};
  app.set_help_flag("--help");
  app.fallthrough();

  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report.");
  app.require_subcommand();
  app.add_option(
    "--seed", seed, "Random seed for the LSH and other parts that require randomness [0].");
  app.callback([&]() {
    if (app.count("--seed")) {
      gen.seed(seed);
    }
  });
  app.add_option(
    "--num-threads", num_threads, "Number of threads to use in OpenMP-based parallelism [1].");

  auto& sc_index = *app.add_subcommand("index", "Build an index from k-mers of reference genomes.");
  auto& sc_place =
    *app.add_subcommand("place", "Place queries on a tree with respect to an index.");
  auto& sc_dist =
    *app.add_subcommand("dist", "Estimate distances of queries to genomes in an index.");
  auto& sc_inspect =
    *app.add_subcommand("inspect", "Display statistics and information for a given index.");
  auto& sc_sketch =
    *app.add_subcommand("sketch", "Create a sketch from k-mers in a single FASTA/FASTQ file.");
  auto& sc_seek =
    *app.add_subcommand("seek", "Seek query sequences in a sketch and estimate distances.");

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
    krepp_place.ensure_wbackbone();
    std::cerr << "Placing given sequences on the backbone tree..." << std::endl;
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
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
    std::cerr << "Done skething & saving, elapsed: " << es_s.count() << " sec" << std::endl;
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
