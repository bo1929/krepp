#include "krepp.hpp"
#include "common.hpp"

void Bkrepp::set_lshf() { lshf = std::make_shared<LSHF>(k, h, m); }

void Bkrepp::parse_newick_tree()
{
  tree = std::make_shared<Tree>();
  tree->parse(nwk_path);
  tree->reset_traversal();
}

void Bkrepp::read_input_file()
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
      std::cerr << "Failed to read file for reference name to path mapping" << std::endl;
      exit(EXIT_FAILURE);
    }
    name_to_path[name] = input;
  }
  input_file.close();
}

void Bkrepp::build_index()
{
  record_sptr_t record = std::make_shared<Record>(tree);
  root_dynht = std::make_shared<DynHT>(nrows, tree, record);
  omp_set_num_threads(num_threads);
  omp_set_nested(1);
#pragma omp parallel
  {
#pragma omp single
    {
      build_for_subtree(tree->get_root(), root_dynht);
    }
  }
}

void Bkrepp::save_index()
{
  assertm(root_dynht->get_nkmers() > 0, "No k-mers to to save!");
  flatht_sptr_t root_flatht = std::make_shared<FlatHT>(root_dynht);
  root_flatht->save(index_dir, suffix);
  root_flatht->get_tree()->save(index_dir, suffix);
  root_flatht->get_crecord()->save(index_dir, suffix);
}

void Bkrepp::build_for_subtree(node_sptr_t nd, dynht_sptr_t dynht)
{
  if (nd->check_leaf()) {
    sh_t sh = nd->get_sh();
    if (name_to_path.find(nd->get_name()) != name_to_path.end()) {
      rseq_sptr_t rs = std::make_shared<RSeq>(w, r, frac, sh, lshf, name_to_path[nd->get_name()]);
      dynht->fill_table(rs);
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
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
      children_dynht_v.emplace_back(std::make_shared<DynHT>(nrows, tree, dynht->get_record()));
#pragma omp task untied shared(dynht)
      {
        build_for_subtree(*std::next(nd->get_children(), i), children_dynht_v[i]);
        omp_set_lock(&parent_lock);
        dynht->union_table(children_dynht_v[i]);
        omp_unset_lock(&parent_lock);
      }
    }
#pragma omp taskwait
    omp_destroy_lock(&parent_lock);
#pragma omp critical
    {
      std::cerr << "\33[2K\r" << std::flush;
      std::cerr << "Internal node: " << nd->get_name() << "\tsize: " << dynht->get_nkmers()
                << "\tprogress: " << (++build_count) << "/" << tree->get_nnodes() << "\r"
                << std::flush;
    }
  }
}

void Bkrepp::save_metadata()
{
  std::filesystem::path metadata_path = index_dir / ("metadata" + suffix);
  std::ofstream metadata_stream(metadata_path, std::ofstream::binary);
  metadata_stream.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(&frac), sizeof(frac));
  metadata_stream.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  metadata_stream.write(reinterpret_cast<char*>(lshf->ppos_data()), (h) * sizeof(uint8_t));
  metadata_stream.write(reinterpret_cast<char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
  if (!metadata_stream.good()) {
    std::cerr << "Writing the metadata for the index has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  metadata_stream.close();
}

Bkrepp::Bkrepp(CLI::App& sub_build)
{
  sub_build
    .add_option(
      "-l,--index-dir", index_dir, "Path to the directory in which the index will be stored.")
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
  sub_build.add_option("-k,--kmer-len", k, "Length of k-mers [30].");
  sub_build.add_option("-w,--win-len", w, "Length of minimizer window (w>k) [k+3].");
  sub_build.add_option("-h,--num-positions", h, "Number of positions for the LSH [14].");
  sub_build.add_option("-m,--modulo-lsh", m, "Mudulo value to partition LSH space [2].");
  sub_build.add_option(
    "-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m [1].");
  sub_build.add_flag(
    "--frac,!--no-frac", frac, "If given, k-mers with r > LSH(x) mod m will be excluded [false].");
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
    std::filesystem::create_directory(index_dir);
    suffix = "-";
    /* suffix += "k" + std::to_string(k) + "w" + std::to_string(w) + "h" + std::to_string(h); */
    suffix += "m" + std::to_string(m) + "r" + std::to_string(r);
    suffix += frac ? "-frac" : "-no_frac";
  });
}

void WLkrepp::load_index()
{
  node_phmap<std::string, std::set<std::string>> suffix_to_ltype;
  for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
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
  }
  std::set<std::string> lall{"cmer", "crecord", "inc", "metadata", "tree"};
#pragma omp parallel for num_threads(num_threads), schedule(static)
  for (uint32_t lix = 0; lix < suffixes.size(); ++lix) {
    if (suffix_to_ltype[suffixes[lix]] == lall) {
      index->add_partial_index(suffixes[lix]);
    } else {
      std::cerr << "There is a partial index with a missing file!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

void Qkrepp::header_dreport(strstream& dreport_stream)
{
  dreport_stream << "#software: krepp\t#version: " VERSION "\t#invocation :" + invocation;
  dreport_stream << "\nSEQ_ID\tREFERENCE_NAME\tDIST\n";
}

void Qkrepp::estimate_distances()
{
  strstream dreport_stream;
  header_dreport(dreport_stream);
  std::cout << dreport_stream.rdbuf();
  omp_set_num_threads(num_threads);
  qseq_sptr_t qs = std::make_shared<QSeq>(query_path);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      while (qs->read_next_batch() || !qs->is_batch_finished()) {
        QBatch qb(index, qs, hdist_th);
#pragma omp task untied
        {
          qb.estimate_distances();
        }
      }
#pragma omp taskwait
    }
  }
}

void Qkrepp::end_jplace(strstream& jplace_stream)
{
  jplace_stream << "\t\t\t{\"p\" : [ ], \"n\" : [\"NaN\"]}\n";
  jplace_stream << "\t],\n";
  jplace_stream << "\t\"tree\" : \"";
  index->get_tree()->stream_newick_str(jplace_stream, index->get_tree()->get_root());
  jplace_stream << "\"\n}";
}

void Qkrepp::begin_jplace(strstream& jplace_stream)
{
  jplace_stream
    << "{\n\t\"version\" : 3,\n\t"
       "\"fields\" : [\"edge_num\", \"like_weight_ratio\", \"likelihood\", \"pavg_distance\", \"pendant_length\", \"distal_length\"],\n"
       "\t\"metadata\" : {\n"
       "\t\t\"software\" : \"krepp\",\n"
       "\t\t\"version\" : \"" VERSION "\",\n"
       "\t\t\"repository\" : \"https://github.com/bo1929/krepp\",\n"
       "\t\t\"invocation\" : \"";
  jplace_stream << invocation;
  jplace_stream << "\"\n\t},\n\t\"placements\" :\n\t\t[\n";
}

void Qkrepp::place_sequences()
{
  strstream jplace_stream;
  begin_jplace(jplace_stream);
  std::cout << jplace_stream.rdbuf();
  omp_set_num_threads(num_threads);
  qseq_sptr_t qs = std::make_shared<QSeq>(query_path);
#pragma omp parallel shared(qs)
  {
#pragma omp single
    {
      while (qs->read_next_batch() || !qs->is_batch_finished()) {
        QBatch qb(index, qs, hdist_th);
#pragma omp task untied
        {
          qb.place_sequences();
        }
      }
#pragma omp taskwait
    }
  }
  jplace_stream.str("");
  end_jplace(jplace_stream);
  std::cout << jplace_stream.rdbuf();
}

Ikrepp::Ikrepp(CLI::App& sub_info)
{
  sub_info
    .add_option("-l,--index-dir", index_dir, "Path to the directory containing reference index.")
    ->required()
    ->check(CLI::ExistingDirectory);
  sub_info.callback([&]() { index = std::make_shared<Index>(index_dir); });
}

void Ikrepp::display_info() { index->display_info(); }

Qkrepp::Qkrepp(CLI::App& sub_query)
{
  sub_query
    .add_option(
      "-l,--index-dir", index_dir, "Path to the directory containing the reference index.")
    ->required()
    ->check(CLI::ExistingDirectory);
  /* sub_query.add_option( */
  /*   "-o,--output-dir", output_dir, "Path to the directory to output results [./]."); */
  sub_query.add_option("-q,--query-file", query_path, "Path to FASTA/FASTQ query file.")
    ->required()
    ->check(CLI::ExistingFile);
  sub_query.add_option(
    "--hdist-th", hdist_th, "The maximum Hamming distance for a k-mer to match [4].");
  sub_query.add_option(
    "--leave-out-ref",
    leave_out_ref,
    "The reference taxon to be excluded during query, useful for benchmarking and testing.");
  sub_query.callback([&]() {
    /* std::filesystem::create_directory(output_dir); */
    index = std::make_shared<Index>(index_dir);
  });
}

int main(int argc, char** argv)
{
  std::ios::sync_with_stdio(false);
  CLI::App app{"Keremet: "
               "a tool for k-mer-based search in large genome collections & metagenomic analysis."};
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

  auto& sub_build = *app.add_subcommand("build", "Build an index for k-mers of reference genomes.");
  Bkrepp b(sub_build);

  auto& sub_place = *app.add_subcommand(
    "place", "Place query sequences on the backbone tree using reference k-mers.");
  Qkrepp qplace(sub_place);

  auto& sub_dist = *app.add_subcommand(
    "dist", "Estimate distances between query sequences and matching references.");
  Qkrepp qdist(sub_dist);

  auto& sub_info =
    *app.add_subcommand("info", "Display statistics and information for a given index.");
  Ikrepp i(sub_info);

  CLI11_PARSE(app, argc, argv);
  for (int i = 0; i < argc; ++i) {
    invocation += std::string(argv[i]) + " ";
  }
  invocation.pop_back();

  auto tstart = std::chrono::system_clock::now();
  std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  std::cerr << std::ctime(&tstart_f) << "\n";

  if (sub_build.parsed()) {
    std::cerr << "Reading the tree and initializing the index..." << std::endl;
    b.set_lshf();
    b.read_input_file();
    b.parse_newick_tree();

    std::cerr << "Building the index..." << std::endl;
    b.build_index();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    std::cerr << "\nFinished building, elapsed: " << es_b.count() << " seconds" << std::endl;

    b.save_index();
    b.save_metadata();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "\nDone converting & saving, elapsed: " << es_s.count() << " seconds" << std::endl;
  }

  if (sub_place.parsed()) {
    std::cerr << "Loading the index and the tree..." << std::endl;
    qplace.load_index();
    std::cerr << "Placing given sequences on the backbone tree..." << std::endl;
    qplace.place_sequences();
  }

  if (sub_dist.parsed()) {
    std::cerr << "Loading the index and the tree..." << std::endl;
    qdist.load_index();
    std::cerr << "Estimating distances between given sequences and references..." << std::endl;
    qdist.estimate_distances();
  }

  if (sub_info.parsed()) {
    i.load_index();
    i.display_info();
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << "\n" << std::ctime(&tend_f);

  return 0;
}
