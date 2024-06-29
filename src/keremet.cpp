#include "keremet.hpp"

void Bkrmt::build_library() {
  DynTable root_dyntable(nrows, ref_tree->get_record());
  omp_set_num_threads(num_threads);
  omp_set_nested(1);
#pragma omp parallel
  {
#pragma omp single
    { build_for_subtree(ref_tree->get_root(), root_dyntable); }
  }
}

void Bkrmt::build_for_subtree(node_sptr_t nd, DynTable &dt) {
  if (nd->check_leaf()) {
    sh_t shash = nd->get_shash();
    if (name_to_gpath.find(nd->get_name()) != name_to_gpath.end()) {
      std::string gpath = name_to_gpath[nd->get_name()];
      refseq_sptr_t rs =
          std::make_shared<RefSeq>(k, w, shash, gpath, hash_func);
      dt.fill_table(rs);
#pragma omp critical
      {
        std::cout << "Genome processed: " << nd->get_name() << "\t";
        dt.print_info();
      }
    } else {
#pragma omp critical
      { std::cout << "Genome skipped: " << nd->get_name() << std::endl; }
    }
  } else {
    assert(nd->get_nchildren() > 0);
    vec<node_sptr_t> children_nd_v = nd->get_children();
    vec<DynTable> children_dt_v;
    children_dt_v.assign(nd->get_nchildren(),
                         DynTable(nrows, ref_tree->get_record()));
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
    for (tuint i = 0; i < nd->get_nchildren(); ++i) {
#pragma omp task untied shared(dt)
      {
        build_for_subtree(children_nd_v[i], children_dt_v[i]);
        omp_set_lock(&parent_lock);
        dt.union_table(children_dt_v[i]);
        omp_unset_lock(&parent_lock);
      }
    }
#pragma omp taskwait
#pragma omp critical
    {
      std::cout << "Internal node processed: " << nd->get_name() << "\t";
      dt.print_info();
    }
  }
}

void Bkrmt::parse_newick_tree() {
  ref_tree = std::make_shared<Tree>(nwk_filepath);
  ref_tree->parse();
  ref_tree->reset_traversal();
}

void Bkrmt::read_input_file() {
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string gpath, name;
    if (!(std::getline(iss, name, '\t') && std::getline(iss, gpath, '\t'))) {
      std::cerr << "Failed to read file for mapping of reference names to paths"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    name_to_gpath[name] = gpath;
  }
  input_file.close();
  nleaves = name_to_gpath.size();
}

void Bkrmt::get_random_positions(uint8_t k, uint8_t h,
                                 std::vector<uint8_t> &npos_v,
                                 std::vector<uint8_t> &ppos_v) {
  uint8_t n;
  assert(h <= 16);
  assert(h < k);
  std::uniform_int_distribution<uint8_t> distrib(0, k - 1);
  for (uint8_t m = 0; m < h; m++) {
    n = distrib(gen);
    if (std::count(ppos_v.begin(), ppos_v.end(), n)) {
      m -= 1;
    } else {
      ppos_v.push_back(n);
    }
  }
  std::sort(ppos_v.begin(), ppos_v.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < k; ++i) {
    if (i != ppos_v[ix_pos])
      npos_v.push_back(i);
    else
      ix_pos++;
  }
}

Bkrmt::Bkrmt(CLI::App &app) {
  CLI::App *sub_build = app.add_subcommand(
      "build", "Builds a library using k-mers of reference genomes.");
  sub_build
      ->add_option("-l,--library-dir", library_dir,
                   "Path to the directory in which the library is stored.")
      ->required();
  sub_build
      ->add_option("-i,--input-file", input_filepath,
                   "Path to the tsv-file containing paths/urls and names of "
                   "references.")
      ->required()
      ->check(CLI::ExistingFile);
  sub_build
      ->add_option("-t,--nwk-file", nwk_filepath,
                   "Path to the Newick file for the reference tree.")
      ->required()
      ->check(CLI::ExistingFile);
  sub_build->add_option("-k,--kmer-length", k,
                        "Length of k-mers. Default: 29.");
  sub_build->add_option("-w,--window-length", w,
                        "Length of minimizer window. Default: k+3.");
  sub_build->add_option("-h,--num-positions", h,
                        "Number of positions for the LSH. Default: 13.");
  nrows = pow(2, 2 * h);
}

void Bkrmt::set_hash_func() {
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  get_random_positions(k, h, npos_v, ppos_v);
  hash_func = std::make_shared<LSHF>(ppos_v, npos_v);
}

Qkrmt::Qkrmt(CLI::App &app) {
  CLI::App *sub_query = app.add_subcommand(
      "query", "Query given sequences with respect to reference libraries.");
  sub_query
      ->add_option(
          "-l,--library-dir", library_dir_v,
          "Path(s) to the directory containing reference library. "
          "Note that multiple libraries could be given to this option.")
      ->required();
  sub_query
      ->add_option("-o,--output-dir", output_dir,
                   "Path to the directory to output query results. "
                   "Default: the current working directory.")
      ->check(CLI::ExistingDirectory);
  sub_query
      ->add_option("-q,--query-file", query_file,
                   "Path to query FASTA/FASTQ files.")
      ->required()
      ->check(CLI::ExistingFile);
}

int main(int argc, char **argv) {
  CLI::App app{"Keremet: "
               "a tool for k-mer-based search in large genome collections & "
               "metagenomic analysis!"};
  app.set_help_flag("--help");
  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose,
               "Increased verbosity and progress report.");
  app.require_subcommand();
  uint32_t seed = 0;
  app.add_option(
      "--seed", seed,
      "Random seed for the LSH and other parts that require randomness.");
  app.callback([&]() {
    if (app.count("--seed"))
      gen.seed(seed);
  });
  app.add_option("--num-threads", num_threads,
                 "Number of threads to use in OpenMP-based parallelism.");
  Bkrmt b(app);
  Qkrmt q(app);
  // sub_build->callback([&]() {
  //   if (!(sub_build->count("-w") + sub_build->count("--window-length")))
  //     w = k + 3;
  // });

  CLI11_PARSE(app, argc, argv);

  b.set_hash_func();
  b.read_input_file();
  b.parse_newick_tree();
  std::cout << "Building the library..." << std::endl;
  b.build_library();

  return 0;
}
