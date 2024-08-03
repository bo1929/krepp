#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include <cstdint>
#include <unordered_map>

struct match_t
{
  node_sptr_t nd;
  uint32_t hdist;
  enc_t enc_lr;
  match_t(node_sptr_t nd, uint32_t hdist, enc_t enc_lr)
    : nd(nd)
    , hdist(hdist)
    , enc_lr(enc_lr)
  {}
};

struct sinfo_t
{
  float coverage_mer = 0;
  float coverage_pos = 0;
  uint32_t parsimony_score = 0;
  uint32_t min_hdist = std::numeric_limits<uint32_t>::max();
};

struct cdist_t
{
  uint32_t sub_hdist = std::numeric_limits<uint32_t>::max();
  uint32_t exc_hdist = std::numeric_limits<uint32_t>::max();
};

class QMers
{
  friend class QBatch;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t max_hdist = 3);
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void compute_coverage();
  void print_coverage();
  void compute_pseudoparsimony();
  std::pair<node_sptr_t, uint32_t> max_parsimonious_mer();

private:
  uint8_t k;
  uint32_t len;
  uint32_t max_hdist;
  tree_sptr_t tree;
  library_sptr_t library;
  crecord_sptr_t crecord;
  vec<vec<match_t>> match_vvec;
  std::unordered_map<node_sptr_t, sinfo_t> node_to_sinfo;
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs);
  void search_batch(uint32_t max_hdist);
  void search_mers(char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);

private:
  uint8_t k;
  uint32_t m;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t batch_size;
  lshf_sptr_t lshashf;
  library_sptr_t library;
  vec<std::string> name_batch;
  vec<std::string> seq_batch;
};

#endif
