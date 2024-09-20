#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

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
  uint32_t zc;
  uint32_t hdist;
  match_t(enc_t enc_lr, uint32_t pos, uint32_t zc, uint32_t hdist)
    : enc_lr(enc_lr)
    , pos(pos)
    , zc(zc)
    , hdist(hdist)
  {}
};

class QMers : public std::enable_shared_from_this<QMers>
{
  friend class QBatch;
  friend class Minfo;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th = 3, double min_covpos = 0.5);
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void summarize_mers();
  qmers_sptr_t getptr() { return shared_from_this(); }

private:
  uint32_t k;
  uint32_t h;
  uint32_t len;
  uint32_t nmers;
  uint32_t hdist_th = 0;
  double min_covpos = 0.0;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  library_sptr_t library = nullptr;
  parallel_flat_phmap<node_sptr_t, minfo_uptr_t> node_to_minfo = {};
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs);
  void search_batch(uint32_t hdist_th, double min_covpos);
  void print_summary(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix);
  void place_wrt_closest(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix);
  void place_wrt_tau(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix);
  double corr_dist_blen(node_sptr_t nd_c, qmers_sptr_t qmers_placement);
  void print_matches(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix);
  void search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);

private:
  uint32_t k;
  uint32_t m;
  uint64_t mask_bp;
  uint64_t mask_lr;
  uint64_t batch_size;
  tree_sptr_t tree;
  lshf_sptr_t lshf;
  library_sptr_t library;
  vec<std::string> name_batch;
  vec<std::string> seq_batch;
};

class Minfo
{
  friend class QBatch;

public:
  Minfo(qmers_sptr_t qmers, float rho)
    : qmers(qmers)
    , rho(rho)
  {
    homoc_v.resize(qmers->len, 0);
    subsc_v.resize(qmers->len, 0);
    obsbp_v.resize(qmers->len, 0);
    hdisthist_v.resize(qmers->hdist_th + 1, 0);
    match_v.clear();
  }
  void update_match(enc_t enc_lr, uint32_t pos, uint32_t zc, uint32_t curr_hdist)
  {
    if (match_v.empty() || ((match_v.back()).pos != pos)) {
      match_v.emplace_back(enc_lr, pos, zc, curr_hdist);
      match_count++;
    } else {
      if ((match_v.back()).hdist > curr_hdist) {
        (match_v.back()).enc_lr = enc_lr;
        (match_v.back()).zc = zc;
        (match_v.back()).hdist = curr_hdist;
      }
    }
  }
  void estimate_distance(optimize::HDistHistLLH& llhfunc);
  void summarize_matches();

private:
  float rho = 1;
  double covmer = 0.0;
  double covpos = 0.0;
  uint32_t match_count = 0;
  double d_was = std::numeric_limits<double>::max();
  double d_llh = std::numeric_limits<double>::max();
  double avg_hdist = std::numeric_limits<double>::max();
  uint32_t sup_hdist = std::numeric_limits<uint32_t>::max();
  std::vector<match_t> match_v;
  std::vector<uint32_t> homoc_v;
  std::vector<uint32_t> subsc_v;
  std::vector<uint32_t> obsbp_v;
  std::vector<uint32_t> hdisthist_v;
  qmers_sptr_t qmers = nullptr;
};

#endif
