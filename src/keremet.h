#ifndef _KEREMET_H
#define _KEREMET_H

#include "CLI11.hpp"
#include "MurmurHash3.h"
#include "assess.h"
#include "common.h"
#include "encode.h"
#include "io.h"
#include "library.h"
#include "lsh.h"
#include "query.h"
#include "table.h"
#include "taxonomy.h"
#include "tree.h"
#include <cstdint>
#include <fstream>

void report_ghashes(um_str_h_t &gID2ghash) {
  for (auto kv : gID2ghash) {
    std::cout << "Genome " << kv.first << " : " << kv.second << std::endl;
  }
  std::cout << "Number of genomes : " << gID2ghash.size() << std::endl;
}

void get_set(um_h_esset &shash2desc,
             std::unordered_map<h_t, std::string> &shash2names,
             std::unordered_set<std::string> &sh_name, h_t k) {
  sh_name.clear();
  if (!shash2names.contains(k)) {
    auto &el = shash2desc[k];
    std::unordered_set<std::string> ccs1;
    std::unordered_set<std::string> ccs2;
    get_set(shash2desc, shash2names, ccs1, el.ss);
    get_set(shash2desc, shash2names, ccs2, k - el.ss);
    sh_name.merge(ccs1);
    sh_name.merge(ccs2);
  } else {
    sh_name.insert(shash2names[k]);
  }
}

void report_shashes(um_h_esset &shash2desc, um_uint32_h_t &ghash2gID) {
  for (auto kv : shash2desc) {
    std::cout << "ss:" << kv.second.ss << ", a:" << kv.second.a << std::endl;
    std::cout << kv.first << "(" << kv.second.c << ")"
              << " :";
    std::unordered_set<h_t> sg_s;
    /* get_set(shash2desc, sg_s, kv.first); */
    for (auto h : sg_s) {
      std::cout << " " << ghash2gID[h];
    }
    std::cout << std::endl;
  }
}

void read_refmap(std::unordered_map<std::string, std::string> &le_fpaths,
                 std::string input_filepath) {
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string fpath, le;
    if (!(std::getline(iss, le, '\t') && std::getline(iss, fpath, '\t'))) {
      std::cerr << "Failed to read file for taxon ID to input map!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    le_fpaths[le] = fpath;
  }
  input_file.close();
}

void get_randpos(uint8_t k, uint8_t h, vec_uint8 &ppos, vec_uint8 &npos) {
  uint8_t n;
  assert(h <= 16);
  assert(h < k);
  std::uniform_int_distribution<uint8_t> distrib(0, k - 1);
  for (uint8_t m = 0; m < h; m++) {
    n = distrib(gen);
    if (count(ppos.begin(), ppos.end(), n)) {
      m -= 1;
    } else {
      ppos.push_back(n);
    }
  }
  sort(ppos.begin(), ppos.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < k; ++i) {
    if (i != ppos[ix_pos])
      npos.push_back(i);
    else
      ix_pos++;
  }
}

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

#endif
