#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "library.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

struct minfo_t;
typedef std::shared_ptr<minfo_t> minfo_sptr_t;
typedef std::unique_ptr<minfo_t> minfo_uptr_t;

struct match_t
{
  uint32_t pos;
  uint32_t zc;
  uint32_t hdist;
  match_t(uint32_t pos, uint32_t zc, uint32_t hdist)
    : pos(pos)
    , zc(zc)
    , hdist(hdist)
  {}
};

struct minfo_t
{
  std::vector<match_t> match_v;
  std::vector<uint32_t> homoc_v;
  std::vector<uint32_t> subsc_v;
  std::vector<uint32_t> hdisthist_v;
  uint32_t match_count = 0;
  float covmer = 0.0;
  float covpos = 0.0;
  float wschdist = std::numeric_limits<float>::max();
  float avghdist = std::numeric_limits<float>::max();
  uint32_t maxhdist = std::numeric_limits<uint32_t>::max();
  minfo_t(uint32_t len)
  {
    homoc_v.resize(len, 0);
    subsc_v.resize(len, 0);
  }
  void update_match(uint32_t pos, uint32_t zc, uint32_t curr_hdist)
  {
    if (match_v.empty() || ((match_v.back()).pos != pos)) {
      match_v.emplace_back(pos, zc, curr_hdist);
      match_count++;
    } else {
      if ((match_v.back()).hdist > curr_hdist) {
        (match_v.back()).zc = zc;
        (match_v.back()).hdist = curr_hdist;
      }
    }
  }
  void print_info()
  {
    std::cout << covpos << "\t" << covmer << "\t" << wschdist << "\t" << avghdist << "\t"
              << match_count << std::endl;
  }
};

class QMers
{
  friend class QBatch;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th = 3, float min_covpos = 0.5);
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void summarize_matches();
  void print_matches(const std::string& name);
  void print_coverage(const std::string& name);
  void print_dist(const std::string& name);
  void print_summary(const std::string& name);

private:
  uint8_t k;
  uint32_t len;
  uint32_t hdist_th = 0;
  float min_covpos = 0.0;
  tree_sptr_t tree = nullptr;
  lshf_sptr_t lshf = nullptr;
  library_sptr_t library = nullptr;
  flat_phmap<node_sptr_t, minfo_uptr_t> node_to_minfo;
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs);
  void search_batch(uint32_t hdist_th, float min_covpos);
  void print_summary(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix);
  void search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc);

private:
  uint8_t k;
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

#endif
