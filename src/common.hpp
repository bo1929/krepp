#ifndef _COMMON_H
#define _COMMON_H

#include "omp.h"
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <curl/curl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <random>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
// TODO: Organize all included headers across all headers.

extern uint32_t num_threads;
extern thread_local std::random_device rd;
extern thread_local std::mt19937 gen;
extern const unsigned char seq_nt4_table[128];
extern const uint64_t nt4_lr_table[4];
extern const uint64_t nt4_bp_table[4];

class LSHF;
class Tree;
class Node;
class Subset;
class Record;
class RefSeq;
class DynTable;

typedef std::shared_ptr<LSHF> lshf_sptr_t;
typedef std::shared_ptr<Tree> tree_sptr_t;
typedef std::shared_ptr<Node> node_sptr_t;
typedef std::shared_ptr<Subset> subset_sptr_t;
typedef std::shared_ptr<Record> record_sptr_t;
typedef std::shared_ptr<RefSeq> refseq_sptr_t;

typedef uint64_t sh_t;
typedef uint32_t enc_t;
typedef uint_least32_t tuint;

template<typename T>
using vvec = std::vector<std::vector<T>>;
template<typename T>
using vec = std::vector<T>;

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

static inline uint32_t hd_lr64(const uint64_t x, const uint64_t y)
{
  uint64_t z1 = x ^ y;
  uint32_t z2 = z1 >> 32;
  uint32_t zc = z1 | z2;
  return __builtin_popcount(zc);
}

static inline uint32_t hd_lr32(const uint32_t x, const uint32_t y)
{
  uint32_t z1 = x ^ y;
  uint16_t z2 = z1 >> 16;
  uint16_t zc = z1 | z2;
  return __builtin_popcount(zc);
}

static inline u_int64_t revcomp_b64(const u_int64_t& x, size_t k)
{
  uint64_t res = ~x;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  return (res >> (2 * (32 - k)));
}

static inline uint64_t rmoddp_b64(uint64_t x)
{
  x = x & 0x5555555555555555;
  x = (x | (x >> 1)) & 0x3333333333333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff0000ffff;
  x = (x | (x >> 16)) & 0x00000000ffffffff;
  return x;
}

static inline u_int64_t cast_bp64_lr64(u_int64_t x)
{
  return (rmoddp_b64(x >> 1) << 32) | rmoddp_b64(x);
}

struct mer_t
{
  enc_t encoding;
  sh_t shash;
  mer_t(enc_t encoding, sh_t shash)
    : encoding(encoding)
    , shash(shash)
  {}
};

#define assertm(exp, msg) assert(((void)msg, exp))

#endif
