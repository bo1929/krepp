#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

struct match_t
{
  enc_t enc_lr;
  uint32_t pos;
  uint32_t hdist;
  se_t se;
  match_t(enc_t enc_lr, se_t se, uint32_t pos, uint32_t hdist)
    : enc_lr(enc_lr)
    , se(se)
    , pos(pos)
    , hdist(hdist)
  {}
};

class QMers
{
  friend class QBatch;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t max_hdist = 3);
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void compute_coverage();
  bool check_better_coverage(qmers_sptr_t qmers);
  void print_info(std::string name);
  void place_baseline(); // TODO: Temporary, remove.

private:
  std::unordered_map<uint32_t, vec<match_t>> pos_to_matches;
  std::unordered_map<se_t, vec<match_t>> se_to_matches;
  std::unordered_map<se_t, float> se_to_coverage;
  float total_coverage = 0.0;
  uint32_t max_hdist;
  uint32_t len;
  uint8_t k;
  tree_sptr_t tree;
  library_sptr_t library;
  crecord_sptr_t crecord;
  /* std::unordered_map<uint32_t, float> hdist_to_coverage; */
  /* std::unordered_map<uint32_t, float> hdist_to_ratio; */
  /* std::unordered_map<uint32_t, uint32_t> hdist_to_nnodes; */
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
