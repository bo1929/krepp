#ifndef _COMPARE_H
#define _COMPARE_H

#include "common.hpp"
#include "sketch.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"
#include <boost/math/tools/minima.hpp>

namespace optimize {
  class HDistHistLLH;
}

class CBatch
{
public:
  CBatch(sketch_sptr_t sketch, qseq_sptr_t qs, uint32_t hdist_th);
  void estimate_distances(std::ostream& output_stream);
  void search_mers(const char* seq, uint64_t len);
  void add_or_mer(uint32_t rix, enc_t enc_lr)
  {
    uint32_t hdist_curr;
    uint32_t hdist_min = hdist_th + 1;
    std::vector<enc_t>::const_iterator iter1 = sketch->bucket_start();
    std::vector<enc_t>::const_iterator iter2 = sketch->bucket_next();
    for (; iter1 < iter2; ++iter1) {
      hdist_curr = popcount_lr32((*iter1) ^ enc_lr);
      if (hdist_curr < hdist_min) {
        hdist_min = hdist_curr;
      }
    }
    if (hdist_min <= hdist_th) {
      or_mismatch_count--;
      or_match_count++;
      or_hdisthist_v[hdist_min]++;
    }
  }
  void add_rc_mer(uint32_t rix, enc_t enc_lr)
  {
    uint32_t hdist_curr;
    uint32_t hdist_min = hdist_th + 1;
    std::vector<enc_t>::const_iterator iter1 = sketch->bucket_start();
    std::vector<enc_t>::const_iterator iter2 = sketch->bucket_next();
    for (; iter1 < iter2; ++iter1) {
      hdist_curr = popcount_lr32((*iter1) ^ enc_lr);
      if (hdist_curr < hdist_min) {
        hdist_min = hdist_curr;
      }
    }
    if (hdist_min <= hdist_th) {
      rc_mismatch_count--;
      rc_match_count++;
      rc_hdisthist_v[hdist_min]++;
    }
  }

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint32_t hdist_th;
  lshf_sptr_t lshf;
  sketch_sptr_t sketch;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint32_t enmers;
  uint64_t batch_size;
  uint64_t bix;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  optimize::HDistHistLLH llhfunc;
  double rho;

  double or_match_count;
  double or_mismatch_count;
  std::vector<double> or_hdisthist_v;
  double or_d_llh = std::numeric_limits<double>::quiet_NaN();
  double or_v_llh = std::numeric_limits<double>::quiet_NaN();
  double rc_match_count;
  double rc_mismatch_count;
  std::vector<double> rc_hdisthist_v;
  double rc_d_llh = std::numeric_limits<double>::quiet_NaN();
  double rc_v_llh = std::numeric_limits<double>::quiet_NaN();
};

#endif
