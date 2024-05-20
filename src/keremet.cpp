#include "keremet.h"
#include <cstdint>
#include <unordered_set>

#define FAST

struct esset {
  uint32_t ss;
  uint32_t tt;
  uint32_t c;
  uint32_t a;
  uint32_t ie = 0;
  uint32_t ke = 0;
};

typedef std::vector<std::vector<std::pair<encT, uint32_t>>> vvec_merrec;
typedef std::unordered_map<std::string, uint32_t> um_str_uint32;
typedef std::unordered_map<uint32_t, std::string> um_uint32_str;
typedef std::unordered_map<uint32_t, esset> um_uint32_esset;

void report_histogram(vvec_merrec &table) {
  std::map<size_t, uint32_t> hist_map;
  for (auto &row : table) {
    hist_map[row.size()]++;
  }
  std::cout << "Row size histogram:" << std::endl;
  for (auto kv : hist_map) {
    float rpercent =
        static_cast<float>(kv.second) / static_cast<float>(table.size());
    /* std::cout << "\t" << kv.first << " : " << rpercent << std::endl; */
    std::cout << "\t" << kv.first << " : " << kv.second << std::endl;
  }
}

void report_ghashes(um_str_uint32 &gID2ghash) {
  for (auto kv : gID2ghash) {
    std::cout << "Genome " << kv.first << " : " << kv.second << std::endl;
  }
  std::cout << "Number of genomes : " << gID2ghash.size() << std::endl;
}

void get_set(um_uint32_esset &shash2desc, std::unordered_set<uint32_t> &sg_s,
             uint32_t k) {
  sg_s.clear();
  auto el = shash2desc[k];
  if (el.a > 1) {
    std::unordered_set<uint32_t> ccs1;
    std::unordered_set<uint32_t> ccs2;
    get_set(shash2desc, ccs1, el.ss);
    get_set(shash2desc, ccs2, el.tt);
    sg_s.merge(ccs1);
    sg_s.merge(ccs2);
  } else {
    sg_s.insert(el.tt);
  }
}

void report_shashes(um_uint32_esset &shash2desc, um_uint32_str &ghash2gID) {
  for (auto kv : shash2desc) {
    std::vector<uint32_t> sg_v;
    std::cout << "ss:" << kv.second.ss << ", tt:" << kv.second.tt
              << ", ie:" << kv.second.ie << ", c:" << kv.second.c
              << ", a:" << kv.second.a << std::endl;
    auto el = kv.second;
    sg_v.push_back(el.tt);
    while (el.a > 1) {
      el = shash2desc[el.ss];
      sg_v.push_back(el.tt);
    }
    std::reverse(sg_v.begin(), sg_v.end());
    std::cout << kv.first << "(" << kv.second.c << ")"
              << " :";
    for (auto h : sg_v) {
      std::cout << " " << ghash2gID[h];
    }
    std::cout << std::endl;
  }
}

int main(int argc, char **argv) {
  PRINT_VERSION

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

  // read reference genome mapping and assign a hash to each genome
  std::unordered_map<std::string, std::string> gID2gpath;
  um_str_uint32 gID2ghash;
  um_uint32_str ghash2gID;
  read_refmap(gID2gpath, gID2ghash, ghash2gID, input_filepath);

  // get LSH positions
  vec_uint8 ppos;
  vec_uint8 npos;
  maskLSH lsh_vg;
  get_randpos(k, h, ppos, npos);
  lsh_vg = generateMaskLSH(ppos);

  uint32_t num_rows = pow(2, 2 * h);
  vvec_merrec table;
  table.resize(num_rows);

  um_uint32_esset shash2desc;

  std::vector<std::pair<std::string, std::string>> kvs{gID2gpath.begin(),
                                                       gID2gpath.end()};
#pragma omp parallel for schedule(dynamic) num_threads(16)
  for (unsigned int fix = 0; fix < gID2gpath.size(); ++fix) {
    std::string gID;
    gID = kvs[fix].first;
    uint32_t ghash = gID2ghash[gID];
    vec_str filepath_v;
    filepath_v.push_back(kvs[fix].second);
    inputHandler<encT> pI(filepath_v, k, w, h, &lsh_vg, &npos);
    uint64_t total_genome_len = pI.extractInput(1);
#pragma omp critical
    {
      for (auto &lsh_enc : pI.lsh_enc_vec) {
        bool seen = false;
        for (unsigned int i = 0; i < table[lsh_enc.first].size() && !seen;
             ++i) {
          if (lsh_enc.second == table[lsh_enc.first][i].first) {
            uint32_t shash = table[lsh_enc.first][i].second;
            uint32_t new_shash;
            /* uint64_t conc_sghash = 0; */
            /* conc_sghash += shash; */
            /* conc_sghash= (conc_sghash << 32) + ghash; */
            /* MurmurHash3_x86_32(&conc_sghash, 8, 0, &new_shash); */
            new_shash = shash + ghash;

            shash2desc[shash].c--;

            bool existing_hash = shash2desc.contains(new_shash);
            if (existing_hash) {
              shash2desc[new_shash].c++;
            } else {
              shash2desc[new_shash].tt = ghash;
              shash2desc[new_shash].ss = shash;
              shash2desc[new_shash].c = 1;
              shash2desc[new_shash].a = 1 + shash2desc[shash].a;
              shash2desc[shash].ie++;
              if (shash2desc.contains(ghash)) {
                shash2desc[ghash].ie++;
              } else {
                shash2desc[ghash].tt = ghash;
                shash2desc[ghash].ss = ghash;
                shash2desc[ghash].c = 1;
                shash2desc[ghash].a = 1;
                shash2desc[ghash].ie = 1;
              }
            }
            seen = true;
            table[lsh_enc.first][i].second = new_shash;
          }
        }
        if (!seen) {
          table[lsh_enc.first].push_back(std::make_pair(lsh_enc.second, ghash));
          if (shash2desc.contains(ghash)) {
            shash2desc[ghash].c++;
          } else {
            shash2desc[ghash].tt = ghash;
            shash2desc[ghash].ss = ghash;
            shash2desc[ghash].c = 1;
            shash2desc[ghash].a = 1;
            shash2desc[ghash].ie = 1;
          }
        }
      }
    }
  }

  /* report_histogram(table); */
  /* report_ghashes(gID2ghash); */
  /* report_shashes(shash2desc, ghash2gID); */

  // for (auto it1 = shash2desc.begin(); it1 != shash2desc.end(); ++it1) {
  //   if ((it1->second.a > 1)) {
  //     for (auto it2 = shash2desc.begin();
  //          (it2 != shash2desc.end()) && (it1->first > it2->first); ++it2) {
  //       if ((it2->second.a > 1) &&
  //           ((it1->second.c > 0) || (it2->second.c > 0)) &&
  //           shash2desc.contains(it1->first + it2->first)) {
  //         shash2desc[shash2desc[it1->first + it2->first].ss].ie--;
  //         shash2desc[shash2desc[it1->first + it2->first].tt].ie--;
  //         it1->second.ie++;
  //         it2->second.ie++;
  //         shash2desc[it1->first + it2->first].ss = it1->first;
  //         shash2desc[it1->first + it2->first].tt = it2->first;
  //         assert((it2->second.a + it1->second.a) ==
  //                shash2desc[it1->first + it2->first].a);
  //         std::cout << shash2desc[it1->first + it2->first].a << "="
  //                   << it2->second.a << "+" << it1->second.a << std::endl;
  //       }
  //     }
  //   }
  // }

  //  for (auto it1 = shash2desc.begin(); it1 != shash2desc.end(); ++it1) {
  //    if ((it1->second.c > 0) && (it1->second.a > 1)) {
  //      for (auto it2 = shash2desc.begin();
  //           (it2 != shash2desc.end()) && (it1->first > it2->first); ++it2) {
  //        if ((it2->second.c > 0) && (it2->second.a > 1)) {
  //          std::unordered_set<uint32_t> sg_s1;
  //          std::unordered_set<uint32_t> sg_s2;
  //          get_set(shash2desc, sg_s1, it1->first);
  //          get_set(shash2desc, sg_s2, it2->first);
  //          sg_s1.merge(sg_s2);
  //          uint32_t sum = std::accumulate(sg_s1.begin(), sg_s1.end(), 0);
  //          if (shash2desc.contains(sum)) {
  //            if ((it1->second.a < shash2desc[sum].a) &&
  //                (it2->second.a < shash2desc[sum].a)) {
  //              shash2desc[shash2desc[sum].ss].ie--;
  //              shash2desc[shash2desc[sum].tt].ie--;
  //              it1->second.ie++;
  //              it2->second.ie++;
  //              shash2desc[sum].ss = it1->first;
  //              shash2desc[sum].tt = it2->first;
  //              std::cout << shash2desc[sum].a << "=" << it2->second.a << "+"
  //                        << it1->second.a << std::endl;
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }

  uint32_t nonzero = 0;
  std::map<uint32_t, std::vector<uint32_t>, std::greater<int>> m_card_sh;
  for (auto it = shash2desc.begin(); it != shash2desc.end(); ++it) {
    m_card_sh[it->second.a].push_back(it->first);
    nonzero += (it->second.c > 0);
  }
  // for (auto kv : m_card_sh) {
  //   for (auto v : kv.second) {
  //     if (shash2desc[v].c > 0 || shash2desc[v].ke > 0) {
  //       shash2desc[shash2desc[v].ss].ke++;
  //       shash2desc[shash2desc[v].tt].ke++;
  //       shash2desc[v].ke++;
  //     }
  //   }
  // }

  // std::cout << "Number of subsets seen : " << shash2desc.size() << std::endl;
  // for (auto kv : m_card_sh) {
  //   for (auto v : kv.second) {
  //     if (shash2desc[v].ke == 0) {
  //       shash2desc[shash2desc[v].ss].ie--;
  //       shash2desc[shash2desc[v].tt].ie--;
  //       shash2desc.erase(v);
  //     }
  //   }
  // }
  std::cout << "Number of nonzero : " << nonzero << std::endl;
  std::cout << "Number of subsets kept : " << shash2desc.size() << std::endl;

  return 0;
}
