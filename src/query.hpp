#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

namespace optimize {
  class HDistHistLLH;
  void simulate_hdhistllh();
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
  QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th = 3, double min_gamma = 0.5);
  qmers_sptr_t getptr() { return shared_from_this(); }
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void summarize_minfo();

private:
  uint32_t k;
  uint32_t h;
  uint32_t len;
  uint32_t nmers;
  uint32_t hdist_th = 0;
  double min_gamma = 0.0;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  library_sptr_t library = nullptr;
  parallel_flat_phmap<node_sptr_t, minfo_uptr_t> node_to_minfo = {};
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs);
  void place_batch(uint32_t hdist_th, double min_gamma);
  void search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  node_sptr_t place_wrt_closest(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);
  node_sptr_t place_wrt_tau(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);

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
  friend class QMers;
  friend class QBatch;

public:
  Minfo(qmers_sptr_t qmers, double rho)
    : qmers(qmers)
    , rho(rho)
  {
    hdisthist_v.resize(qmers->hdist_th + 1, 0);
    match_v.clear();
  }
  void update_match(enc_t enc_lr, uint32_t pos, uint32_t curr_hdist)
  {
    if (match_v.empty() || ((match_v.back()).pos != pos)) {
      match_count++;
      hdisthist_v[curr_hdist]++;
      match_v.emplace_back(enc_lr, pos, curr_hdist);
    } else {
      if ((match_v.back()).hdist > curr_hdist) {
        hdisthist_v[curr_hdist]++;
        hdisthist_v[(match_v.back()).hdist]--;
        (match_v.back()).enc_lr = enc_lr;
        (match_v.back()).hdist = curr_hdist;
      }
    }
  }
  void estimate_distance(optimize::HDistHistLLH& llhfunc);
  void compute_gamma();

private:
  double rho = 1;
  double gamma = 0.0;
  uint32_t match_count = 0;
  double d_llh = std::numeric_limits<double>::max();
  std::vector<match_t> match_v;
  std::vector<uint32_t> hdisthist_v;
  qmers_sptr_t qmers = nullptr;
};

#endif
