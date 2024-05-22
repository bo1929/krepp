#include "keremet.h"
#include <cstdint>
#include <unordered_set>

int main(int argc, char **argv) {
  // {{{ Initialization
  PRINT_VERSION
  std::random_device rd;
  std::mt19937 g(rd());
  // }}}

  // {{{ Definitions and arguments for CLI
  CLI::App app{
      "Keremet: "
      "a tool for k-mer-based search across a large genome collection!"};
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
                   "Path to the directory containing the library.")
      ->required();
  std::string input_filepath;
  sub_build
      ->add_option(
          "-i,--input-file", input_filepath,
          "Path to the tsv-file containing paths and IDs of references.")
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
  uint8_t b = 16;
  sub_build->add_option("-b,--num-columns", b,
                        "Number of columns of the table. Default: 16.");
  sub_build->callback([&]() {
    if (!(sub_build->count("-w") + sub_build->count("--window-length")))
      w = k + 3;
  });

  CLI::App *sub_query = app.add_subcommand(
      "query",
      "Query given sequences with respect to given reference libraries.");
  std::vector<std::string> library_dir_v;
  sub_query
      ->add_option("-l,--library-dir", library_dir_v,
                   "Path(s) to the directory containing the library. "
                   "Note that multiple libraries can be given to this option.")
      ->required();
  std::string output_dir = "./";
  sub_query
      ->add_option("-o,--output-dir", output_dir,
                   "Path to the directory to output result files. "
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

  // {{{ Getting LSH positions
  vec_uint8 ppos;
  vec_uint8 npos;
  maskLSH lsh_vg;
  get_randpos(k, h, ppos, npos);
  lsh_vg = generateMaskLSH(ppos);
  // }}}

  // {{{ Parsing the Newick tree
  std::ifstream tree_file(nwk_filepath.c_str());
  std::string newick_tree((std::istreambuf_iterator<char>(tree_file)),
                          std::istreambuf_iterator<char>());
  // fprintf(stderr, "newick: %s", newick_tree.c_str());
  std::vector<std::string> n_vec;
  mysplit2(strdup(newick_tree.c_str()), n_vec);
  // for (int i = 0; 0 && i < n_vec.size(); i++)
  //   fprintf(stderr, "%d) %s\n", i, n_vec[i].c_str());
  node_t *root = parse(n_vec);
  // node_t **lst = new node_t *[serial];
  // serialize(root, lst);
  // fprintf(stderr, "serial: %d\n", serial);
  // for (int i = 0; i < serial; i++)
  //   print_node(stderr, lst[i]);
  // }}}

  // {{{ Read reference genomes and correspoding paths
  std::unordered_map<std::string, std::string> le_fpaths;
  read_refmap(le_fpaths, input_filepath);
  // }}}

  // {{{ Initialize main data structures
  libtab root_lt;
  root_lt.k = k;
  root_lt.w = w;
  root_lt.h = h;
  root_lt.lsh_vg = lsh_vg;
  root_lt.npos = npos;
  um_h_esset shash2desc;
  // }}}

  uint64_t MAXGN = 600;
  std::vector<node_t *> po_vec;
  post_order_traversal(root, po_vec, root_lt, shash2desc, le_fpaths, MAXGN);

  po_vec.resize(MAXGN);
  // std::shuffle(po_vec.begin(), po_vec.end(), g);

  // int num_gread = 0;
  // for (auto &nd_ : po_vec) {
  //   std::cout << nd_->name << ", hash: " << nd_->sh << std::endl;
  //   if (le_fpaths.contains(nd_->name)) {
  //     inputHandler<encT> pI({le_fpaths[nd_->name]}, k, w, h, &lsh_vg, &npos);
  //     uint64_t total_genome_len = pI.extractInput(1);
  //     num_gread++;
  //     process_genome(nd_, pI, root_lt.table, shash2desc);
  //   }
  //   std::cout << "genome read: " << num_gread
  //             << ", record size:" << shash2desc.size() << std::endl;
  // }

  std::unordered_map<h_t, std::string> shash2names;
  for (auto nd_ : po_vec) {
    shash2names[nd_->sh] = nd_->name;
  }

  std::cout << "Number of subsets seen : " << shash2desc.size() << std::endl;

  for (auto &nd_ : po_vec) {
    shash2desc.erase(nd_->sh);
  }

  std::cout << "Number of records excluding tree : " << shash2desc.size()
            << std::endl;

  // std::unordered_map<h_t, h_t> shash2rr, shash2ss;
  // for (auto nd_1 : po_vec) {
  //   for (auto nd_2 : po_vec) {
  //     if ((nd_1->sh > nd_2->sh) && (nd_1->nbranches > nd_2->nbranches) &&
  //         shash2desc.contains(nd_1->sh - nd_2->sh)) {
  //       shash2desc.erase(nd_1->sh - nd_2->sh);
  //       shash2rr[nd_1->sh - nd_2->sh] = nd_2->sh;
  //     }
  //     if ((nd_1->sh > nd_2->sh) && shash2desc.contains(nd_1->sh + nd_2->sh))
  //     {
  //       shash2desc.erase(nd_1->sh + nd_2->sh);
  //       shash2ss[nd_1->sh + nd_2->sh] = nd_1->sh;
  //     }
  //   }
  // }

  // std::cout << "Number of records removing trivials : " << shash2desc.size()
  //           << std::endl;

  uint32_t nonzero = 0;
  std::map<uint32_t, std::vector<h_t>, std::greater<uint32_t>> m_card_sh;
  for (auto it = shash2desc.begin(); it != shash2desc.end(); ++it) {
    m_card_sh[it->second.a].push_back(it->first);
    nonzero += (it->second.c > 0);
  }

  std::cout << "Number of nonzero : " << nonzero << std::endl;

  for (auto kv : m_card_sh) {
    for (auto v : kv.second) {
      if (shash2desc[v].c > 0 || shash2desc[v].ke > 0) {
        if (shash2desc.contains(shash2desc[v].ss))
          shash2desc[shash2desc[v].ss].ke++;
        if (shash2desc.contains(v - shash2desc[v].ss))
          shash2desc[v - shash2desc[v].ss].ke++;
        shash2desc[v].ke++;
      }
    }
  }
  for (auto kv : m_card_sh) {
    for (auto v : kv.second) {
      if (shash2desc[v].ke == 0) {
        shash2desc.erase(v);
      }
    }
  }

  std::cout << "Number of subsets kept : " << shash2desc.size() << std::endl;

  for (auto kv : shash2desc) {
    std::cout << kv.first << ": ";
    std::unordered_set<std::string> sh_name;
    get_set(shash2desc, shash2names, sh_name, kv.first);
    for (auto hv : sh_name)
      std::cout << " " << hv;
    std::cout << std::endl;
  }

  return 0;

  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //

  // std::vector<std::pair<std::string, std::string>> kvs{gID2gpath.begin(),
  //                                                      gID2gpath.end()};
  // #pragma omp parallel for schedule(dynamic) num_threads(16)
  // for (unsigned int fix = 0; fix < gID2gpath.size(); ++fix) {
  //   std::string gID;
  //   gID = kvs[fix].first;
  //   h_t ghash = gID2ghash[gID];
  //   vec_str filepath_v;
  //   filepath_v.push_back(kvs[fix].second);
  //   inputHandler<encT> pI(filepath_v, k, w, h, &lsh_vg, &npos);
  //   uint64_t total_genome_len = pI.extractInput(1);
  // #pragma omp critical
  //   {
  //   }
  // }

  // /* report_histogram(table); */
  // /* report_ghashes(gID2ghash); */
  // /* report_shashes(shash2desc, ghash2gID); */

  // // for (auto it1 = shash2desc.begin(); it1 != shash2desc.end(); ++it1) {
  // //   for (auto it2 = shash2desc.begin();
  // //        (it2 != shash2desc.end()) && (it1->first > it2->first); ++it2)
  // {
  // //     if (shash2desc.contains(it1->first + it2->first)) {
  // //       auto el = shash2desc[it1->first + it2->first];
  // //       auto d = (el.ss > el.tt) ? el.ss - el.tt : el.tt - el.ss;
  // //       auto dn = ((it1->second).a > (it2->second).a)
  // //                     ? (it1->second).a - (it2->second).a
  // //                     : (it2->second).a - (it1->second).a;
  // //       /* if (d > dn) { */
  // //       shash2desc[shash2desc[it1->first + it2->first].ss].ie--;
  // //       shash2desc[shash2desc[it1->first + it2->first].tt].ie--;
  // //       it1->second.ie++;
  // //       it2->second.ie++;
  // //       shash2desc[it1->first + it2->first].ss = it1->first;
  // //       shash2desc[it1->first + it2->first].tt = it2->first;
  // //       assert((it2->second.a + it1->second.a) ==
  // //              shash2desc[it1->first + it2->first].a);
  // //       std::cout << shash2desc[it1->first + it2->first].a << "="
  // //                 << it2->second.a << "+" << it1->second.a << std::endl;
  // //       /* } */
  // //     }
  // //   }
  // // }

  // //  for (auto it1 = shash2desc.begin(); it1 != shash2desc.end(); ++it1) {
  // //    if ((it1->second.c > 0) && (it1->second.a > 1)) {
  // //      for (auto it2 = shash2desc.begin();
  // //           (it2 != shash2desc.end()) && (it1->first > it2->first);
  // ++it2)
  // {
  // //        if ((it2->second.c > 0) && (it2->second.a > 1)) {
  // //          std::unordered_set<h_t> sg_s1;
  // //          std::unordered_set<h_t> sg_s2;
  // //          get_set(shash2desc, sg_s1, it1->first);
  // //          get_set(shash2desc, sg_s2, it2->first);
  // //          sg_s1.merge(sg_s2);
  // //          h_t sum = std::accumulate(sg_s1.begin(), sg_s1.end(), 0);
  // //          if (shash2desc.contains(sum)) {
  // //            if ((it1->second.a < shash2desc[sum].a) &&
  // //                (it2->second.a < shash2desc[sum].a)) {
  // //              shash2desc[shash2desc[sum].ss].ie--;
  // //              shash2desc[shash2desc[sum].tt].ie--;
  // //              it1->second.ie++;
  // //              it2->second.ie++;
  // //              shash2desc[sum].ss = it1->first;
  // //              shash2desc[sum].tt = it2->first;
  // //              std::cout << shash2desc[sum].a << "=" << it2->second.a <<
  // "+"
  // //                        << it1->second.a << std::endl;
  // //            }
  // //          }
  // //        }
  // //      }
  // //    }
  // //  }

  /* report_shashes(shash2desc, ghash2gID); */

  return 0;
}
