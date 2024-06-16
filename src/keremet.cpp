#include "keremet.hpp"

int main(int argc, char **argv) {
  // {{{ Definitions and arguments for CLI
  CLI::App app{"Keremet: "
               "a tool for k-mer-based search in large genome collections & "
               "metagenomic analysis!"};
  app.set_help_flag("--help");
  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose,
               "Increased verbosity and progress report.");
  app.require_subcommand();
  int seed = 0;
  app.add_option(
      "--seed", seed,
      "Random seed for the LSH and other parts that require randomness.");
  app.callback([&]() {
    if (app.count("--seed"))
      gen.seed(seed);
  });
  app.add_option("--num-threads", num_threads,
                 "Number of threads to use for OpenMP-based parallelism.");

  CLI::App *sub_build = app.add_subcommand(
      "build", "Builds a library using k-mers of reference genomes.");
  std::string library_dir;
  sub_build
      ->add_option("-l,--library-dir", library_dir,
                   "Path to the directory in which the library is stored.")
      ->required();
  std::string input_filepath;
  sub_build
      ->add_option(
          "-i,--input-file", input_filepath,
          "Path to the tsv-file containing paths/urls and IDs of references.")
      ->required()
      ->check(CLI::ExistingFile);
  std::string nwk_filepath;
  sub_build
      ->add_option("-t,--nwk-file", nwk_filepath,
                   "Path to the Newick file for the reference tree.")
      ->required()
      ->check(CLI::ExistingFile);
  uint8_t k = 29;
  sub_build->add_option("-k,--kmer-length", k,
                        "Length of k-mers. Default: 29.");
  uint8_t w = k + 3;
  sub_build->add_option("-w,--window-length", w,
                        "Length of minimizer window. Default: k+3.");
  uint8_t h = 13;
  sub_build->add_option("-h,--num-positions", h,
                        "Number of positions for the LSH. Default: 13.");
  sub_build->callback([&]() {
    if (!(sub_build->count("-w") + sub_build->count("--window-length")))
      w = k + 3;
  });

  CLI::App *sub_query = app.add_subcommand(
      "query", "Query given sequences with respect to reference libraries.");
  std::vector<std::string> library_dir_v;
  sub_query
      ->add_option(
          "-l,--library-dir", library_dir_v,
          "Path(s) to the directory containing reference library. "
          "Note that multiple libraries could be given to this option.")
      ->required();
  std::string output_dir = "./";
  sub_query
      ->add_option("-o,--output-dir", output_dir,
                   "Path to the directory to output query results. "
                   "Default: the current working directory.")
      ->check(CLI::ExistingDirectory);
  std::string query_file;
  sub_query
      ->add_option("-q,--query-file", query_file,
                   "Path to query FASTA/FASTQ files.")
      ->required()
      ->check(CLI::ExistingFile);

  CLI11_PARSE(app, argc, argv);
  // }}}

  return 0;
}
