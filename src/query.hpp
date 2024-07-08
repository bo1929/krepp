#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

struct match_t
{
  se_t se;
  enc_t enc_lr;
  uint32_t hdist;
  match_t(enc_t enc_lr, se_t se, uint32_t hdist)
    : enc_lr(enc_lr)
    , se(se)
    , hdist(hdist)
  {}
};

class QMers
{
public:
  QMers(library_sptr_t library, uint32_t max_hdist)
    : library(library)
    , max_hdist(max_hdist){};
  QMers(library_sptr_t library)
    : library(library)
    , max_hdist(5){};
  void add_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);

private:
  vec<match_t> match_v;
  library_sptr_t library;
  uint32_t max_hdist;
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
