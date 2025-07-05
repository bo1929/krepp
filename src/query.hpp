#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "index.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"

namespace optimize {
  class HDistHistLLH;
}

class Minfo;
class IMers;
typedef std::shared_ptr<Minfo> minfo_sptr_t;
typedef std::unique_ptr<Minfo> minfo_uptr_t;
typedef std::shared_ptr<IMers> imers_sptr_t;

class IMers : public std::enable_shared_from_this<IMers>
{
  friend class IBatch;

public:
  IMers(index_sptr_t index, uint64_t len, uint32_t hdist_th);
  imers_sptr_t getptr() { return shared_from_this(); }
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
  index_sptr_t index = nullptr;
  uint32_t hdist_filt = std::numeric_limits<uint32_t>::max();
  parallel_flat_phmap<node_sptr_t, minfo_sptr_t> leaf_to_minfo = {};
};

class IBatch
{
public:
  IBatch(index_sptr_t index,
         qseq_sptr_t qs,
         uint32_t hdist_th,
         double dist_max,
         uint32_t tau,
         bool no_filter,
         bool multi = false);
  void search_mers(const char* seq, uint64_t len, imers_sptr_t imers_or, imers_sptr_t imers_rc);
  void summarize_matches(imers_sptr_t imers_or, imers_sptr_t imers_rc);
  void estimate_distances(std::ostream& output_stream);
  void place_sequences(std::ostream& output_stream);
  void report_distances(strstream& batch_stream);
  void report_placement(strstream& batch_stream);

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint32_t hdist_th;
  double dist_max;
  bool no_filter;
  uint32_t tau;
  tree_sptr_t tree;
  lshf_sptr_t lshf;
  index_sptr_t index;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint32_t enmers;
  uint32_t onmers;
  uint32_t wnmers_or;
  uint32_t wnmers_rc;
  uint64_t batch_size;
  uint64_t bix;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  node_sptr_t nd_closest = nullptr;
  minfo_sptr_t mi_closest = nullptr;
  optimize::HDistHistLLH llhfunc;
  bool multi = false;

protected:
  parallel_flat_phmap<node_sptr_t, minfo_sptr_t> node_to_minfo = {};
};

class Minfo
{
  friend class IBatch;

  // struct match_t
  // {
  //   enc_t enc_lr;
  //   uint32_t pos;
  //   uint32_t hdist;
  //   match_t(enc_t enc_lr, uint32_t pos, uint32_t hdist)
  //     : enc_lr(enc_lr)
  //     , pos(pos)
  //     , hdist(hdist)
  //   {}
  // };

public:
  Minfo(uint32_t hdist_th, uint32_t nmers, double rho = 0.0)
    : nmers(nmers)
    , rho(rho)
  {
    rcard++;
    rmatch_count = rho > 0 ? 1 : 0;
    mismatch_count = nmers;
    hdisthist_v.resize(hdist_th + 1, 0);
  }
  Minfo(uint32_t hdist_th) { hdisthist_v.resize(hdist_th + 1, 0); }
  void join(minfo_sptr_t minfo)
  {
    double denom = nmers ? 0.5 : 1.0;
    /* gamma = (gamma + minfo->gamma) * denom; */
    match_count = (match_count + minfo->match_count) * denom;
    mismatch_count = (mismatch_count + minfo->mismatch_count) * denom;
    for (uint32_t x = 0; x < hdisthist_v.size(); ++x) {
      hdisthist_v[x] = (hdisthist_v[x] + minfo->hdisthist_v[x]) * denom;
    }
    hdist_min = std::min(hdist_min, minfo->hdist_min);
    nmers = std::max(nmers, minfo->nmers);
    rho = std::max(rho, minfo->rho);
    rcard += minfo->rcard;
    rmatch_count += minfo->rmatch_count;
  }
  void update_match(enc_t enc_lr, uint32_t pos, uint32_t hdist_curr)
  {
    /* if (match_v.empty() || ((match_v.back()).pos != pos)) { */
    if (last_hdist == 0xFFFFFFFF || last_pos != pos) {
      match_count++;
      mismatch_count--;
      hdisthist_v[hdist_curr]++;
      last_pos = pos;
      last_hdist = hdist_curr;
      /* match_v.emplace_back(enc_lr, pos, hdist_curr); */
    } else {
      if (last_hdist > hdist_curr) {
        hdisthist_v[hdist_curr]++;
        hdisthist_v[last_hdist]--;
        last_hdist = hdist_curr;
        /* hdisthist_v[(match_v.back()).hdist]--; */
        /* (match_v.back()).enc_lr = enc_lr; */
        /* (match_v.back()).hdist = hdist_curr; */
      }
    }
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
    }
  }
  double get_leq_tau(uint32_t tau)
  {
    double total_leq_tau = 0.0;
    for (uint32_t x = 0; x <= tau; ++x) {
      total_leq_tau += hdisthist_v[x];
    }
    return total_leq_tau;
  }
  /* void compute_gamma(); */
  void optimize_likelihood(optimize::HDistHistLLH& llhfunc);
  double likelihood_ratio(double d, optimize::HDistHistLLH& llhfunc);

#define PLACEMENT_FIELD(nd, mi)                                                                    \
  "\t\t\t\t[" << (nd->get_se() - 1) << ", 0, " << (nd->get_blen() / 2.0) << ", " << -mi->v_llh     \
              << ", " << exp(-mi->chisq / 2) << ", " << mi->d_llh << "]"

#define DISTANCE_FIELD(nd, mi) nd->get_name() << "\t" << mi->d_llh

private:
  double nmers = 0;
  double mismatch_count = 0;
  double match_count = 0;
  double rho = 0.0;
  /* double gamma = 0.0; */
  uint32_t rcard = 0;
  uint32_t rmatch_count = 0;
  uint32_t last_pos = 0;
  uint32_t last_hdist = 0xFFFFFFFF;
  uint32_t hdist_min = 0xFFFFFFFF;
  std::vector<double> hdisthist_v;
  double chisq = std::numeric_limits<double>::quiet_NaN();
  double v_llh = std::numeric_limits<double>::quiet_NaN();
  double d_llh = std::numeric_limits<double>::max();
  /* std::vector<match_t> match_v; */
};

#endif
