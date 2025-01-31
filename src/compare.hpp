#ifndef _COMPARE_H
#define _COMPARE_H

#include "common.hpp"
#include "sketch.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"

namespace optimize {
  class HDistHistLLH;
}

typedef std::shared_ptr<CSummary> csummary_sptr_t;

class CSummary
{
  friend class CBatch;

public:
  CSummary(uint64_t enmers, uint32_t hdist_th)
    : hdist_th(hdist_th)
    , mismatch_count(enmers)
  {
    hdisthist_v.resize(hdist_th + 1, 0);
  }
  void optimize_likelihood(optimize::HDistHistLLH llhfunc, double rho);
  void add_matching_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr);

private:
  uint32_t hdist_th;
  double mismatch_count;
  double match_count = 0;
  double d_llh = std::numeric_limits<double>::quiet_NaN();
  double v_llh = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> hdisthist_v;
};

class CBatch
{
public:
  CBatch(sketch_sptr_t sketch, qseq_sptr_t qs, uint32_t hdist_th);
  void estimate_distances(std::ostream& output_stream);
  void search_mers(const char* seq, uint64_t len, CSummary& or_summary, CSummary& rc_summary);

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint32_t hdist_th;
  lshf_sptr_t lshf;
  sketch_sptr_t sketch;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t enmers;
  uint64_t batch_size;
  uint64_t bix;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  optimize::HDistHistLLH llhfunc;
  double rho;
};

#endif
