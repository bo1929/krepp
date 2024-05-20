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
#include <cstdint>
#include <fstream>

void read_refmap(std::unordered_map<std::string, std::string> &gID2gpath,
                 std::unordered_map<std::string, uint32_t> &gID2ghash,
                 std::unordered_map<uint32_t, std::string> &ghash2gID,
                 std::string input_filepath) {
  std::set<uint32_t> set_ghash;
  std::ifstream input_file(input_filepath);
  if (!input_file.good()) {
    std::cerr << "Error opening " << input_filepath << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  uint32_t ghash;

  while (std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string gpath, gID;
    if (!(std::getline(iss, gID, '\t') && std::getline(iss, gpath, '\t'))) {
      std::cerr << "Failed to read file for taxon ID to input map!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    gID2gpath[gID] = gpath;
    MurmurHash3_x86_32(gID.c_str(), gID.length(), 0, &ghash);
    uint32_t pk = 0;
    while (set_ghash.find(ghash) != set_ghash.end()) {
      pk++;
      ghash += pow(pk, 2);
    }
    gID2ghash[gID] = ghash;
  }
  input_file.close();
  for (auto &kv : gID2ghash) {
    ghash2gID[kv.second] = kv.first;
  }
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

#endif
