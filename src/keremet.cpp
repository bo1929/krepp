#include "keremet.hpp"

void Bkrmt::build_library()
{
  std::cout << "Building the library..." << std::endl;
  auto start_b = std::chrono::system_clock::now();

  record = std::make_shared<Record>(ref_tree->get_root());
  DynHT root_dt(nrows, ref_tree, record);
  omp_set_num_threads(num_threads);
  omp_set_nested(1);
#pragma omp parallel
  {
#pragma omp single
    {
      build_for_subtree(ref_tree->get_root(), root_dt);
    }
  }

  auto end_b = std::chrono::system_clock::now();
  std::chrono::duration<double> es_b = end_b - start_b;
  std::cout << "Finished building, elapsed: " << es_b.count() << " seconds" << std::endl;

  FlatHT root_flatht(root_dt);
  root_flatht.save(library_dir, suffix);
  root_flatht.get_tree()->save(library_dir, suffix);
  root_flatht.get_crecord()->save(library_dir, suffix);

  auto end_c = std::chrono::system_clock::now();
  std::chrono::duration<double> es_s = end_c - end_b;
  std::cout << "Done converting & saving, elapsed: " << es_s.count() << " seconds" << std::endl;

  std::time_t end_time = std::chrono::system_clock::to_time_t(end_c);
  std::cout << std::ctime(&end_time);
}

void Bkrmt::build_for_subtree(node_sptr_t nd, DynHT& dynht)
{
  if (nd->check_leaf()) {
    sh_t shash = nd->get_shash();
    if (name_to_input.find(nd->get_name()) != name_to_input.end()) {
      rseq_sptr_t rs =
        std::make_shared<RSeq>(w, r, frac, shash, lshashf, name_to_input[nd->get_name()]);
      dynht.fill_table(rs);
#pragma omp critical
      {
        std::cout << "Genome processed: " << nd->get_name() << "\t";
        dynht.print_info();
      }
    } else {
#pragma omp critical
      {
        std::cout << "Genome skipped: " << nd->get_name() << std::endl;
      }
    }
  } else {
    assert(nd->get_nchildren() > 0);
    vec<node_sptr_t> children_nd_v = nd->get_children();
    vec<DynHT> children_dt_v;
    children_dt_v.assign(nd->get_nchildren(), DynHT(nrows, ref_tree, record));
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
#pragma omp task untied shared(dynht)
      {
        build_for_subtree(children_nd_v[i], children_dt_v[i]);
        omp_set_lock(&parent_lock);
        dynht.union_table(children_dt_v[i]);
        omp_unset_lock(&parent_lock);
      }
    }
#pragma omp taskwait
    omp_destroy_lock(&parent_lock);
#pragma omp critical
    {
      std::cout << "Internal node processed: " << nd->get_name() << "\t";
      dynht.print_info();
    }
  }
}

void Bkrmt::save_metadata()
{
  std::filesystem::path metadata_path = library_dir / ("metadata" + suffix);
  std::ofstream metadata_stream(metadata_path, std::ofstream::binary);
  metadata_stream.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(&frac), sizeof(frac));
  metadata_stream.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(lshashf->ppos_data()), (h) * sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(lshashf->npos_data()), (k - h) * sizeof(uint8_t));
  if (!metadata_stream.good()) {
    std::cerr << "Writing the metadata for the library has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  metadata_stream.close();
}

void Bkrmt::parse_newick_tree()
{
  ref_tree = std::make_shared<Tree>();
  ref_tree->parse(nwk_path);
  ref_tree->reset_traversal();
}

void Bkrmt::read_input_file()
{
  std::ifstream input_file(input_path);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_path << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string input, name;
    if (!(std::getline(iss, name, '\t') && std::getline(iss, input, '\t'))) {
      std::cerr << "Failed to read file for mapping of reference names to paths" << std::endl;
      exit(EXIT_FAILURE);
    }
    name_to_input[name] = input;
  }
  input_file.close();
}

void Bkrmt::set_lshashf() { lshashf = std::make_shared<LSHF>(k, h, m); }

Bkrmt::Bkrmt(CLI::App& sub_build)
{
  sub_build
    .add_option(
      "-l,--library-dir", library_dir, "Path to the directory in which the library is stored.")
    ->required();
  sub_build
    .add_option("-i,--input-file",
                input_path,
                "Path to the tsv-file containing paths/urls and names of references.")
    ->required()
    ->check(CLI::ExistingFile);
  sub_build
    .add_option("-t,--nwk-file", nwk_path, "Path to the Newick file for the reference tree.")
    ->required()
    ->check(CLI::ExistingFile);
  sub_build.add_option("-k,--kmer-len", k, "Length of k-mers. Default: 29.");
  sub_build.add_option("-w,--win-len", w, "Length of minimizer windows (w>k). Default: k+3.");
  sub_build.add_option("-h,--num-positions", h, "Number of positions for the LSH. Default: 13.");
  sub_build.add_option(
    "-m,--modulo-lsh", m, "Mudulo value to partition LSH space, must be smaller. Default: 2.");
  sub_build.add_option(
    "-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m. Default: 1.");
  sub_build.add_flag(
    "--frac,!--no-frac", frac, "If --frac, all k-mers with r < LSH(x) mod m will be excluded.");
  sub_build.callback([&]() {
    uint32_t hash_size = pow(2, 2 * h);
    uint32_t full_residue = hash_size % m;
    if (frac) {
      nrows = (hash_size / m) * r;
      nrows = full_residue >= r ? nrows + r : nrows + full_residue;
    } else {
      nrows = (hash_size / m);
      nrows = full_residue > r ? nrows + 1 : nrows;
    }
    if (!(sub_build.count("-w") + sub_build.count("--win-len"))) {
      w = k + 3;
    }
    std::filesystem::create_directory(library_dir);
    suffix = "-";
    /* suffix += "k" + std::to_string(k) + "w" + std::to_string(w) + "h" + std::to_string(h); */
    suffix += "m" + std::to_string(m) + "r" + std::to_string(r);
    suffix += frac ? "-frac" : "-no_frac";
  });
}

void Library::add_partial_flatht(std::string suffix)
{
  std::filesystem::path metadata_path = library_dir / ("metadata" + suffix);
  std::ifstream metadata_stream(metadata_path, std::ifstream::binary);
  uint8_t k_lshashf, w, h_lshashf;
  uint32_t m_lshashf, r, nrows_partial;
  bool frac;
  metadata_stream.read(reinterpret_cast<char*>(&k_lshashf), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&h_lshashf), sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(&m_lshashf), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  metadata_stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
  metadata_stream.read(reinterpret_cast<char*>(&nrows_partial), sizeof(uint32_t));
  vec<uint8_t> ppos_v(h_lshashf), npos_v(k_lshashf - h_lshashf);
  metadata_stream.read(reinterpret_cast<char*>(ppos_v.data()), ppos_v.size() * sizeof(uint8_t));
  metadata_stream.read(reinterpret_cast<char*>(npos_v.data()), npos_v.size() * sizeof(uint8_t));
  if (!metadata_stream.good()) {
    std::cerr << "Reading the metadata for the partial library has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  metadata_stream.close();

  lshf_sptr_t curr_lshashf = std::make_shared<LSHF>(m_lshashf, ppos_v, npos_v);
#pragma omp critical
  {
    if (lshashf != nullptr) {
      if (!lshashf->check_compatible(curr_lshashf)) {
        std::cerr << "Partial libraries have incompatible hash functions." << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      lshashf = curr_lshashf;
      k = k_lshashf;
      h = h_lshashf;
      m = m_lshashf;
      nrows = pow(2, 2 * h);
    }
  }

  flatht_sptr_t curr_flatht = std::make_shared<FlatHT>(tree, crecord);
  curr_flatht->load(library_dir, suffix);
#pragma omp critical
  {
    if (frac) {
      for (uint32_t ix = 0; ix < r; ++ix) {
        r_to_flatht[ix] = curr_flatht;
      }
    } else {
      r_to_flatht[r] = curr_flatht;
    }
  }
}

void Pkrmt::load_library()
{
  std::unordered_map<std::string, std::set<std::string>> suffix_to_ltype;
  for (const auto& entry : std::filesystem::directory_iterator(library_dir)) {
    std::string filename, ltype, mrcfg, fracv;
    size_t pos1, pos2;
    filename = entry.path().filename();
    pos1 = filename.find("-", 0);
    pos2 = filename.find("-", pos1 + 1);
    ltype = filename.substr(0, pos1);
    mrcfg = filename.substr(pos1, pos2 - pos1);
    fracv = filename.substr(pos2, filename.size() - pos2);
    suffix_to_ltype[mrcfg + fracv].insert(ltype);
  }
  std::vector<std::string> suffixes;
  for (auto const& [suffix, ltypes] : suffix_to_ltype) {
    suffixes.push_back(suffix);
    library->add_partial_tree(suffix);
  }
  for (auto const& [suffix, ltypes] : suffix_to_ltype) {
    library->add_partial_crecord(suffix);
  }
  std::set<std::string> lall{"cmer", "crecord", "inc", "metadata", "tree"};
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (uint32_t lix = 0; lix < suffixes.size(); ++lix) {
    if (suffix_to_ltype[suffixes[lix]] == lall) {
      library->add_partial_flatht(suffixes[lix]);
    } else {
      std::cerr << "There is a partial library with a missing file!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

void Library::add_partial_tree(std::string suffix)
{
  tree_sptr_t curr_tree = std::make_shared<Tree>();
  curr_tree->load(library_dir, suffix);
  if (tree != nullptr) {
    if (!tree->check_compatible(curr_tree)) {
      std::cerr << "Partial libraries are based on different trees." << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    tree = curr_tree;
  }
}

void Library::add_partial_crecord(std::string suffix)
{
  crecord_sptr_t curr_crecord = std::make_shared<CRecord>(tree);
  curr_crecord->load(library_dir, suffix);
  if (crecord != nullptr) {
    if (!crecord->check_compatible(curr_crecord)) {
      std::cerr << "Partial libraries have incompatible subset records." << std::endl;
      exit(EXIT_FAILURE);
    } else {
      crecord->merge(curr_crecord);
    }
  } else {
    crecord = curr_crecord;
  }
}

std::vector<cmer_t>::const_iterator Library::begin(uint32_t rix)
{
  return rix < m ? (r_to_flatht[rix % m])->begin() : (r_to_flatht[rix % m])->at(rix / m - 1);
}

std::vector<cmer_t>::const_iterator Library::end(uint32_t rix)
{
  /* return rix < nrows ? (r_to_flatht[rix % m])->at(rix / m) : (r_to_flatht[rix % m])->end(); */
  return (r_to_flatht[rix % m])->at(rix / m);
}

void Pkrmt::place_sequences()
{
  qseq_sptr_t qs = std::make_shared<QSeq>(query_path);
  while (qs->read_next_batch() || !qs->is_batch_finished()) {
    QBatch qb(library, qs);
    qb.search_batch(5);
  }
}

Pkrmt::Pkrmt(CLI::App& sub_place)
{
  sub_place
    .add_option(
      "-l,--library-dir", library_dir, "Path to the directory containing reference library.")
    ->required()
    ->check(CLI::ExistingDirectory);
  sub_place.add_option("-o,--output-dir",
                       output_dir,
                       "Path to the directory to output place results. "
                       "Default:the current working directory.");
  sub_place
    .add_option(
      "-q,--query-file", query_path, "Path to FASTA/FASTQ query file to place on the tree.")
    ->required()
    ->check(CLI::ExistingFile);
  sub_place.callback([&]() {
    std::filesystem::create_directory(output_dir);
    library = std::make_shared<Library>(library_dir);
  });
}

int main(int argc, char** argv)
{
  CLI::App app{"Keremet: "
               "a tool for k-mer-based search in large genome collections & "
               "metagenomic analysis!"};
  app.set_help_flag("--help");
  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report.");
  app.require_subcommand();
  uint32_t seed = 0;
  app.add_option(
    "--seed", seed, "Random seed for the LSH and other parts that require randomness.");
  app.callback([&]() {
    if (app.count("--seed")) {
      gen.seed(seed);
    }
  });
  app.add_option(
    "--num-threads", num_threads, "Number of threads to use in OpenMP-based parallelism.");

  auto& sub_build =
    *app.add_subcommand("build", "Builds a library using k-mers of reference genomes.");
  Bkrmt b(sub_build);

  auto& sub_place =
    *app.add_subcommand("place", "Place given sequences with respect to reference libraries.");
  Pkrmt p(sub_place);

  CLI11_PARSE(app, argc, argv);

  if (sub_build.parsed()) {
    b.set_lshashf();
    b.read_input_file();
    b.parse_newick_tree();
    b.build_library();
    b.save_metadata();
  }
  if (sub_place.parsed()) {
    p.load_library();
    p.place_sequences();
  }

  return 0;
}
