#ifndef _COMMON_H
#define _COMMON_H

#include "omp.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <curl/curl.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <locale>
#include <math.h>
#include <memory>
#include <numeric>
#include <ostream>
#include <parallel_hashmap/phmap.h>
#include <parallel_hashmap/btree.h>
#include <queue>
#include <random>
#include <regex>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <zlib.h>
#include <stdio.h>

#define VERSION "v0.0.3"

extern uint32_t num_threads;
extern std::string invocation;
extern std::string leave_out_ref;
extern thread_local std::random_device rd;
extern thread_local std::mt19937 gen;
extern const unsigned char seq_nt4_table[128];
extern const uint64_t nt4_lr_table[4];
extern const uint64_t nt4_bp_table[4];

typedef uint64_t sh_t;
typedef uint32_t se_t;
typedef uint64_t inc_t;
typedef uint32_t enc_t;
typedef uint32_t tuint_t;

typedef std::stringstream strstream;

template<typename T>
using vvec = std::vector<std::vector<T>>;
template<typename T>
using vec = std::vector<T>;

class RSeq;
class QSeq;
class QBatch;
class QMers;
class LSHF;
class Tree;
class Node;
class Subset;
class Record;
class CRecord;
class DynHT;
class FlatHT;
class Index;

typedef std::shared_ptr<RSeq> rseq_sptr_t;
typedef std::shared_ptr<QSeq> qseq_sptr_t;
typedef std::shared_ptr<QBatch> qbatch_sptr_t;
typedef std::shared_ptr<QMers> qmers_sptr_t;
typedef std::shared_ptr<LSHF> lshf_sptr_t;
typedef std::shared_ptr<Tree> tree_sptr_t;
typedef std::shared_ptr<Node> node_sptr_t;
typedef std::shared_ptr<Subset> subset_sptr_t;
typedef std::shared_ptr<Record> record_sptr_t;
typedef std::shared_ptr<CRecord> crecord_sptr_t;
typedef std::shared_ptr<DynHT> dynht_sptr_t;
typedef std::shared_ptr<FlatHT> flatht_sptr_t;
typedef std::shared_ptr<Index> index_sptr_t;
typedef std::pair<enc_t, se_t> cmer_t;
typedef std::pair<se_t, se_t> pse_t;

struct mer_t
{
  enc_t encoding = 0;
  sh_t sh = 0;
  mer_t() {}
  mer_t(enc_t encoding, sh_t sh)
    : encoding(encoding)
    , sh(sh)
  {}
};

static inline uint32_t gp_hash(const std::string& str)
{
  uint32_t b = 378551;
  uint32_t a = 63689;
  uint32_t h = 0;
  for (std::size_t i = 0; i < str.length(); i++) {
    h = h * a + str[i];
    a = a * b;
  }
  return (h & 0x7FFFFFFF);
}

static inline uint32_t xur32_hash(uint32_t h)
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

static inline uint64_t xur64_hash(uint64_t h)
{
  h ^= (h >> 33);
  h *= 0xff51afd7ed558ccdL;
  h ^= (h >> 33);
  h *= 0xc4ceb9fe1a85ec53L;
  h ^= (h >> 33);
  return h;
}

static inline uint32_t hdist_lr64(const uint64_t x, const uint64_t y)
{
  uint64_t z1 = x ^ y;
  return __builtin_popcount(z1 | (z1 >> 32));
}

static inline uint32_t zc_lr32(const uint32_t x, const uint32_t y)
{
  uint32_t z1 = x ^ y;
  return (z1 | (z1 >> 16)) & 0x0000ffff;
}

static inline uint32_t hdist_lr32(const uint32_t x, const uint32_t y)
{
  uint32_t z1 = x ^ y;
  return __builtin_popcount((z1 | (z1 >> 16)) & 0x0000ffff);
}

static inline uint32_t popcount_lr32(const uint32_t z)
{
  return __builtin_popcount((z | (z >> 16)) & 0x0000ffff);
}

static inline uint64_t revcomp_bp64(const uint64_t x, uint8_t k)
{
  uint64_t res = ~x;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  return (res >> (2 * (32 - k)));
}

static inline uint64_t rmoddp_bp64(uint64_t x)
{
  x = x & 0x5555555555555555;
  x = (x | (x >> 1)) & 0x3333333333333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff0000ffff;
  x = (x | (x >> 16)) & 0x00000000ffffffff;
  return x;
}

static inline uint64_t conv_bp64_lr64(uint64_t x)
{
  return (rmoddp_bp64(x >> 1) << 32) | rmoddp_bp64(x);
}

static inline void
compute_encoding(const char* s1, const char* s2, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr = 0;
  enc_bp = 0;
  for (; s1 < s2; ++s1) {
    enc_lr <<= 1;
    enc_bp <<= 2;
    enc_bp += nt4_bp_table[seq_nt4_table[*s1]];
    enc_lr += nt4_lr_table[seq_nt4_table[*s1]];
  }
}
static inline void update_encoding(const char* s1, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr <<= 1;
  enc_bp <<= 2;
  enc_lr &= 0xFFFFFFFEFFFFFFFE;
  enc_bp += nt4_bp_table[seq_nt4_table[*s1]];
  enc_lr += nt4_lr_table[seq_nt4_table[*s1]];
}

#define assertm(exp, msg) assert(((void)msg, exp))

#define EXTRAARGS                                                                                  \
  phmap::priv::hash_default_hash<K>, phmap::priv::hash_default_eq<K>,                              \
    std::allocator<std::pair<const K, V>>, 4

template<class K, class V>
using parallel_flat_phmap = phmap::parallel_flat_hash_map<K, V, EXTRAARGS, std::mutex>;

template<class K, class V>
using parallel_node_phmap = phmap::parallel_node_hash_map<K, V, EXTRAARGS, std::mutex>;

template<class K, class V>
using flat_phmap = phmap::flat_hash_map<K, V>;

template<class K, class V>
using btree_phmap = phmap::btree_map<K, V>;

template<class K, class V>
using node_phmap = phmap::node_hash_map<K, V>;

static inline double kendalls_tau(double* v1, double* v2, tuint_t len)
{
  assert(len > 1);
  uint32_t i, j;
  int nC = 0, nD = 0;
  uint32_t n1 = 0, n2 = 0;
  uint32_t n0 = len * (len - 1) / 2;
  for (i = 0; i < (len - 1); i++) {
    for (j = i + 1; j < len; j++) {
      if ((v1[i] > v1[j] && v2[i] > v2[j]) || (v1[i] < v1[j] && v2[i] < v2[j]))
        ++nC;
      if ((v1[i] > v1[j] && v2[i] < v2[j]) || (v1[i] < v1[j] && v2[i] > v2[j]))
        ++nD;
      if (v1[i] == v1[j])
        ++n1;
      if (v2[i] == v2[j])
        ++n2;
    }
  }
  return static_cast<double>(nC - nD) / sqrt((n0 - n1) * (n0 - n2));
}

#endif
