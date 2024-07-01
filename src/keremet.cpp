#include "keremet.hpp"
#include "table.hpp"
#include <string>

void Bkrmt::build_library()
{
  std::cout << "Building the library..." << std::endl;

  auto start_b = std::chrono::system_clock::now();
  DynTable root_dyntable(nrows, ref_tree->get_record());
  omp_set_num_threads(num_threads);
  omp_set_nested(1);
#pragma omp parallel
  {
#pragma omp single
    {
      build_for_subtree(ref_tree->get_root(), root_dyntable);
    }
  }
  auto end_b = std::chrono::system_clock::now();
  std::chrono::duration<double> es_b = end_b - start_b;
  std::cout << "Finished building, elapsed: " << es_b.count() << " seconds" << std::endl;

  auto start_c = std::chrono::system_clock::now();
  FlatTable root_flattable(root_dyntable);
  root_flattable.save(library_dir, suffix);
  ref_tree->save(library_dir, suffix);
  auto end_c = std::chrono::system_clock::now();
  std::chrono::duration<double> es_s = end_c - start_c;
  std::cout << "Done converting & saving, elapsed: " << es_s.count() << " seconds" << std::endl;

  std::time_t end_time = std::chrono::system_clock::to_time_t(end_c);
  std::cout << std::ctime(&end_time);
}

void Bkrmt::build_for_subtree(node_sptr_t nd, DynTable& dt)
{
  if (nd->check_leaf()) {
    sh_t shash = nd->get_shash();
    if (name_to_gpath.find(nd->get_name()) != name_to_gpath.end()) {
      std::string gpath = name_to_gpath[nd->get_name()];
      refseq_sptr_t rs = std::make_shared<RefSeq>(k, w, shash, gpath, hash_func);
      dt.fill_table(rs);
#pragma omp critical
      {
        std::cout << "Genome processed: " << nd->get_name() << "\t";
        dt.print_info();
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
    vec<DynTable> children_dt_v;
    children_dt_v.assign(nd->get_nchildren(), DynTable(nrows, ref_tree->get_record()));
    omp_lock_t parent_lock;
    omp_init_lock(&parent_lock);
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
#pragma omp task untied shared(dt)
      {
        build_for_subtree(children_nd_v[i], children_dt_v[i]);
        omp_set_lock(&parent_lock);
        dt.union_table(children_dt_v[i]);
        omp_unset_lock(&parent_lock);
      }
    }
#pragma omp taskwait
    omp_destroy_lock(&parent_lock);
#pragma omp critical
    {
      std::cout << "Internal node processed: " << nd->get_name() << "\t";
      dt.print_info();
    }
  }
}

void Bkrmt::parse_newick_tree()
{
  ref_tree = std::make_shared<Tree>();
  ref_tree->parse(nwk_filepath);
  ref_tree->reset_traversal();
}

void Bkrmt::read_input_file()
{
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
      std::cerr << "Failed to read file for mapping of reference names to paths" << std::endl;
      exit(EXIT_FAILURE);
    }
    name_to_gpath[name] = gpath;
  }
  input_file.close();
}

Bkrmt::Bkrmt(CLI::App& sub_build)
{
  sub_build
    .add_option(
      "-l,--library-dir", library_dir, "Path to the directory in which the library is stored.")
    ->required();
  sub_build
    .add_option("-i,--input-file",
                input_filepath,
                "Path to the tsv-file containing paths/urls and names of "
                "references.")
    ->required()
    ->check(CLI::ExistingFile);
  sub_build
    .add_option("-t,--nwk-file", nwk_filepath, "Path to the Newick file for the reference tree.")
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

void Bkrmt::set_hash_func() { hash_func = std::make_shared<LSHF>(k, h, m, r, frac); }

Qkrmt::Qkrmt(CLI::App& sub_query)
{
  sub_query
    .add_option("-l,--library-dir",
                library_dir_v,
                "Path(s) to the directory containing reference library. "
                "Note that multiple libraries could be given to this option.")
    ->required();
  sub_query
    .add_option("-o,--output-dir",
                output_dir,
                "Path to the directory to output query results. "
                "Default: the current working directory.")
    ->check(CLI::ExistingDirectory);
  sub_query.add_option("-q,--query-file", query_file, "Path to query FASTA/FASTQ files.")
    ->required()
    ->check(CLI::ExistingFile);
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

  auto& sub_query =
    *app.add_subcommand("query", "Query given sequences with respect to reference libraries.");
  Qkrmt q(sub_query);

  CLI11_PARSE(app, argc, argv);

  b.set_hash_func();
  b.read_input_file();
  b.parse_newick_tree();
  b.build_library();

  return 0;
}
