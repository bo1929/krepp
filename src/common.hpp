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
#include <iostream>
#include <locale>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <queue>
#include <random>
#include <regex>
#include <string>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <zlib.h>

extern uint32_t num_threads;
extern thread_local std::random_device rd;
extern thread_local std::mt19937 gen;
extern const unsigned char seq_nt4_table[128];
extern const uint64_t nt4_lr_table[4];
extern const uint64_t nt4_bp_table[4];

typedef uint64_t sh_t;
typedef uint64_t inc_t;
typedef uint32_t enc_t;
typedef uint_least32_t tuint_t;
typedef uint32_t se_t;

template<typename T>
using vvec = std::vector<std::vector<T>>;
template<typename T>
using vec = std::vector<T>;

class LSHF;
class Tree;
class Node;
class Subset;
class Record;
class CRecord;
class RSeq;
class QSeq;
class DynHT;
class FlatHT;
class QBatch;
class QMers;
class Library;

typedef std::shared_ptr<LSHF> lshf_sptr_t;
typedef std::shared_ptr<Tree> tree_sptr_t;
typedef std::shared_ptr<Node> node_sptr_t;
typedef std::shared_ptr<Subset> subset_sptr_t;
typedef std::shared_ptr<Record> record_sptr_t;
typedef std::shared_ptr<CRecord> crecord_sptr_t;
typedef std::shared_ptr<RSeq> rseq_sptr_t;
typedef std::shared_ptr<QSeq> qseq_sptr_t;
typedef std::shared_ptr<DynHT> dynht_sptr_t;
typedef std::shared_ptr<FlatHT> flatht_sptr_t;
typedef std::shared_ptr<QBatch> qbatch_sptr_t;
typedef std::shared_ptr<QMers> qmers_sptr_t;
typedef std::shared_ptr<Library> library_sptr_t;
typedef std::pair<enc_t, se_t> cmer_t;

struct mer_t
{
  enc_t encoding = 0;
  sh_t shash = 0;
  mer_t() {}
  mer_t(enc_t encoding, sh_t shash)
    : encoding(encoding)
    , shash(shash)
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
  uint32_t z2 = z1 >> 32;
  uint32_t zc = z1 | z2;
  return __builtin_popcount(zc);
}

static inline uint32_t hdist_lr32(const uint32_t x, const uint32_t y)
{
  uint32_t z1 = x ^ y;
  uint16_t z2 = z1 >> 16;
  uint16_t zc = z1 | z2;
  return __builtin_popcount(zc);
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

static inline void compute_encoding(char* s1, char* s2, uint64_t& enc_lr, uint64_t& enc_bp)
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
static inline void update_encoding(char* s1, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr <<= 1;
  enc_bp <<= 2;
  enc_lr &= 0xFFFFFFFEFFFFFFFE;
  enc_bp += nt4_bp_table[seq_nt4_table[*s1]];
  enc_lr += nt4_lr_table[seq_nt4_table[*s1]];
}

#define assertm(exp, msg) assert(((void)msg, exp))

#endif
