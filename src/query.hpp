#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "index.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"

#define CJC 4.0 / 3.0

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
         double chisq_value,
         double dist_max,
         uint32_t tau,
         bool no_filter,
         bool multi,
         bool summarize);
  void search_mers(const char* seq, uint64_t len, imers_sptr_t imers_or, imers_sptr_t imers_rc);
  void summarize_matches(imers_sptr_t imers_or, imers_sptr_t imers_rc);
  void estimate_distances(strstream& batch_stream);
  void report_distances(strstream& batch_stream);
  void place_sequences(strstream& batch_stream, bool tabular);
  void report_placement(strstream& batch_stream, bool tabular);
  const parallel_flat_phmap<node_sptr_t, double>& get_summary() { return node_to_wcount; }

private:
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint32_t hdist_th;
  double chisq_value;
  double dist_max;
  bool no_filter;
  bool summarize;
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
  parallel_flat_phmap<node_sptr_t, double> node_to_wcount = {};
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
  Minfo(uint32_t hdist_th, uint32_t nmers, double rho)
    : nmers(nmers)
    , rho(rho)
  {
    rmatch_count = 1;
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
    rmatch_count += minfo->rmatch_count;
  }
  void add(const minfo_sptr_t& minfo, double denom)
  {
    // gamma = gamma + minfo->gamma * denom;
    mismatch_count = nmers ? mismatch_count : minfo->nmers;
    match_count += minfo->match_count * denom;
    mismatch_count -= minfo->match_count * denom;
    for (uint32_t x = 0; x < hdisthist_v.size(); ++x) {
      hdisthist_v[x] = hdisthist_v[x] + minfo->hdisthist_v[x] * denom;
    }
    hdist_min = std::min(hdist_min, minfo->hdist_min);
    nmers = std::max(nmers, minfo->nmers);
    rho = std::max(rho, minfo->rho);
    rmatch_count++;
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

  std::string get_match_string() const {
    std::ostringstream oss;
    // oss << "[";
    for (size_t i = 0; i < hdisthist_v.size(); ++i) {
      if (i > 0) oss << "\t";
      oss << hdisthist_v[i];
    }
    // oss << "]";
    return oss.str();
  }

  double get_leq_tau(uint32_t tau)
  {
    double total_leq_tau = 0.0;
    for (uint32_t x = 0; x <= tau; ++x) {
      total_leq_tau += hdisthist_v[x];
    }
    return total_leq_tau;
  }
  double jukes_cantor_dist() { return -0.75 * log(1 - CJC * d_llh); }
  /* void compute_gamma(); */
  void optimize_likelihood(optimize::HDistHistLLH& llhfunc);
  double likelihood_ratio(double d, optimize::HDistHistLLH& llhfunc);

#define PLACEMENT_FIELD(nd, mi)                                                                                             \
  "[" << nd->get_en() << ", " << mi->jukes_cantor_dist() - nd->get_midpoint_pendant() << ", " << nd->get_midpoint_pendant() \
      << ", " << -mi->v_llh << ", " << mi->lwr << ", " << mi->d_llh << "]"

#define TABULAR_FIELD(nd, mi) nd->get_name(true) << "\t" << nd->get_en() << "\t" << mi->lwr << "\t" << mi->d_llh

#define MATCH_FIELD(nd, mi) nd->get_name() << "\t" << mi->get_match_string() << "\n"

#define DISTANCE_FIELD(nd, mi) nd->get_name() << "\t" << mi->d_llh << "\n"

private:
  double nmers = 0;
  double mismatch_count = 0;
  double match_count = 0;
  double rho = 0.0;
  /* double gamma = 0.0; */
  uint32_t rmatch_count = 0;
  uint32_t last_pos = 0;
  uint32_t last_hdist = 0xFFFFFFFF;
  uint32_t hdist_min = 0xFFFFFFFF;
  std::vector<double> hdisthist_v;
  double chisq = std::numeric_limits<double>::quiet_NaN();
  double lwr = 1; // std::numeric_limits<double>::quiet_NaN();
  double v_llh = std::numeric_limits<double>::quiet_NaN();
  double d_llh = std::numeric_limits<double>::max();
  /* std::vector<match_t> match_v; */
};

#endif
