#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"

namespace optimize {
  class HDistHistLLH;
}

class Minfo;
typedef std::shared_ptr<Minfo> minfo_sptr_t;
typedef std::unique_ptr<Minfo> minfo_uptr_t;
typedef std::stringstream strstream;

class QMers : public std::enable_shared_from_this<QMers>
{
  friend class QBatch;
  friend class Minfo;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th);
  qmers_sptr_t getptr() { return shared_from_this(); }
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);

private:
  uint32_t k;
  uint32_t h;
  uint32_t len;
  uint32_t hdist_th;
  uint32_t onmers;
  uint32_t enmers;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  library_sptr_t library = nullptr;
  uint32_t hdist_filt = std::numeric_limits<uint32_t>::max();
  parallel_flat_phmap<node_sptr_t, minfo_sptr_t> leaf_to_minfo = {};
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs, uint32_t hdist_th = 4);
  void search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  void summarize_minfo(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  void estimate_distances();
  void place_sequences();
  void report_distances(strstream& batch_report);
  void report_placement(strstream& batch_report);

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint32_t hdist_th;
  uint64_t bix;
  uint64_t batch_size;
  uint32_t onmers_or;
  uint32_t onmers_rc;
  uint32_t enmers;
  tree_sptr_t tree;
  lshf_sptr_t lshf;
  library_sptr_t library;
  node_sptr_t nd_pp = nullptr;
  minfo_sptr_t mi_pp = nullptr;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  optimize::HDistHistLLH llhfunc;

protected:
  parallel_flat_phmap<node_sptr_t, minfo_sptr_t> node_to_minfo = {};
};

class Minfo
{
  friend class QMers;
  friend class QBatch;

  struct match_t
  {
    enc_t enc_lr;
    uint32_t pos;
    uint32_t hdist;
    match_t(enc_t enc_lr, uint32_t pos, uint32_t hdist)
      : enc_lr(enc_lr)
      , pos(pos)
      , hdist(hdist)
    {}
  };

public:
  Minfo(uint32_t nmers, uint32_t hdist_th, double rho)
    : nmers(nmers)
    , hdist_th(hdist_th)
    , rho(rho)
  {
    rmatch_count = 1;
    mismatch_count = nmers;
    hdisthist_v.resize(hdist_th + 1, 0);
  }
  Minfo(uint32_t nmers, uint32_t hdist_th)
    : nmers(nmers)
    , hdist_th(hdist_th)
  {
    mismatch_count = nmers;
    hdisthist_v.resize(hdist_th + 1, 0);
  }
  Minfo(uint32_t hdist_th)
    : nmers(0)
    , hdist_th(hdist_th)
  {
    mismatch_count = nmers;
    hdisthist_v.resize(hdist_th + 1, 0);
  }
  void add(minfo_sptr_t minfo)
  {
    if (nmers) {
      gamma = (gamma + minfo->gamma) / 2.0;
      d_llh = (d_llh + minfo->d_llh) / 2.0;
      match_count = (match_count + minfo->match_count) / 2.0;
      mismatch_count = (mismatch_count + minfo->mismatch_count) / 2.0;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        hdisthist_v[x] = (hdisthist_v[x] + minfo->hdisthist_v[x]) / 2.0;
      }
    } else {
      gamma = minfo->gamma;
      d_llh = minfo->d_llh;
      match_count = minfo->match_count;
      mismatch_count = minfo->mismatch_count;
      for (uint32_t x = 0; x <= hdist_th; ++x) {
        hdisthist_v[x] = minfo->hdisthist_v[x];
      }
    }
    rmatch_count += minfo->rmatch_count;
    rho = std::max(rho, minfo->rho);
    hdist_min = std::min(hdist_min, minfo->hdist_min);
    nmers = std::max(nmers, minfo->nmers);
  }
  void update_match(enc_t enc_lr, uint32_t pos, uint32_t hdist_curr)
  {
    /* if (match_v.empty() || ((match_v.back()).pos != pos)) { */
    if (last_hdist == 0xFFFFFFFF || last_pos != pos) {
      match_count++;
      mismatch_count--;
      hdisthist_v[hdist_curr]++;
      /* match_v.emplace_back(enc_lr, pos, hdist_curr); */
      last_pos = pos;
      last_hdist = hdist_curr;
    } else {
      if (last_hdist > hdist_curr) {
        hdisthist_v[hdist_curr]++;
        /* hdisthist_v[(match_v.back()).hdist]--; */
        hdisthist_v[last_hdist]--;
        last_hdist = hdist_curr;
        /* (match_v.back()).enc_lr = enc_lr; */
        /* (match_v.back()).hdist = hdist_curr; */
      }
    }
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
    }
  }
  /* void compute_gamma(); */
  void optimize_likelihood(optimize::HDistHistLLH& llhfunc);
  double likelihood_ratio(double d, optimize::HDistHistLLH& llhfunc);

private:
  uint32_t hdist_th;
  double nmers;
  double mismatch_count;
  double match_count = 0;
  double gamma = 0.0;
  double rho = 0.0;
  uint32_t last_pos = 0;
  uint32_t last_hdist = 0xFFFFFFFF;
  uint32_t rmatch_count = 0;
  uint32_t hdist_min = std::numeric_limits<uint32_t>::max();
  double chisq = std::numeric_limits<double>::max();
  double d_llh = std::numeric_limits<double>::max();
  double v_llh = std::numeric_limits<double>::min();
  std::vector<double> hdisthist_v;
  /* std::vector<match_t> match_v; */
};

#endif
