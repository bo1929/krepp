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
  uint32_t onmers;
  uint32_t enmers;
  uint32_t hdist_th = 0;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  library_sptr_t library = nullptr;
  uint32_t filt_hdist = std::numeric_limits<uint32_t>::max();
  parallel_flat_phmap<node_sptr_t, minfo_uptr_t> node_to_minfo = {};
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs, uint32_t hdist_th = 3, double min_gamma = 0.5);
  void place_sequences();
  void estimate_distances();
  void search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  void summarize_minfo(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  node_sptr_t place_wrt_closest(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  node_sptr_t place_wrt_tau(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  void
  report_distances(std::stringstream& batch_report, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);

private:
  uint64_t bix;
  uint32_t k;
  uint32_t h;
  uint32_t m;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t batch_size;
  uint32_t hdist_th;
  double min_gamma;
  tree_sptr_t tree;
  lshf_sptr_t lshf;
  library_sptr_t library;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  optimize::HDistHistLLH llhfunc;
};

class Minfo
{
  friend class QMers;
  friend class QBatch;

public:
  Minfo(uint32_t hdist_th, double rho)
    : hdist_th(hdist_th)
    , rho(rho)
  {
    hdisthist_v.resize(hdist_th + 1, 0);
    last_pos = 0;
    /* match_v.clear(); */
  }
  void update_match(enc_t enc_lr, uint32_t pos, uint32_t curr_hdist)
  {
    /* if (match_v.empty() || ((match_v.back()).pos != pos)) { */
    if (last_hdist == 0xFFFFFFFF || last_pos != pos) {
      match_count++;
      hdisthist_v[curr_hdist]++;
      /* match_v.emplace_back(enc_lr, pos, curr_hdist); */
      last_pos = pos;
      last_hdist = curr_hdist;
    } else {
      if (last_hdist > curr_hdist) {
        hdisthist_v[curr_hdist]++;
        /* hdisthist_v[(match_v.back()).hdist]--; */
        hdisthist_v[last_hdist]--;
        last_hdist = curr_hdist;
        /* (match_v.back()).enc_lr = enc_lr; */
        /* (match_v.back()).hdist = curr_hdist; */
      }
    }
    if (curr_hdist < min_hdist) {
      min_hdist = curr_hdist;
    }
  }
  /* void compute_gamma(); */
  void optimize_likelihood(optimize::HDistHistLLH& llhfunc);

private:
  uint32_t hdist_th;
  double gamma = 0.0;
  double rho = 1;
  uint32_t last_pos = 0;
  uint32_t last_hdist = 0xFFFFFFFF;
  uint32_t match_count = 0;
  uint32_t min_hdist = std::numeric_limits<uint32_t>::max();
  double d_llh = std::numeric_limits<double>::max();
  /* std::vector<match_t> match_v; */
  std::vector<uint32_t> hdisthist_v;
};

namespace optimize {
  void simulate_hdhistllh(uint32_t k,
                          uint32_t h,
                          double rho,
                          uint32_t len,
                          uint32_t hdist_curr,
                          uint32_t hdist_th,
                          double min_gamma,
                          uint32_t num_replicates);
}

#endif
