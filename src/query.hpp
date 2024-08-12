#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"

struct minfo_t;
typedef std::shared_ptr<minfo_t> minfo_sptr_t;
struct ninfo_t;
typedef std::shared_ptr<ninfo_t> ninfo_sptr_t;

struct match_t
{
  uint32_t pos;
  enc_t enc_lr;
  uint32_t hdist;
  match_t(uint32_t pos, enc_t enc_lr, uint32_t hdist)
    : pos(pos)
    , enc_lr(enc_lr)
    , hdist(hdist)
  {}
};

struct minfo_t
{
  std::vector<match_t> match_v;
  uint32_t match_count = 0;
  float covmer = 0.0;
  float covpos = 0.0;
  float avghdist = std::numeric_limits<float>::max();
  minfo_t() {}
  void update_match(uint32_t pos, uint32_t curr_hdist, enc_t enc_lr)
  {
    if (match_v.empty() || (match_v.back()).pos != pos) {
      match_v.emplace_back(pos, enc_lr, curr_hdist);
      match_count++;
    } else {
      if ((match_v.back()).hdist > curr_hdist) {
        (match_v.back()).enc_lr = enc_lr;
        (match_v.back()).hdist = curr_hdist;
      }
    }
  }
  void print_info()
  {
    std::cout << covpos << "\t" << covmer << "\t" << avghdist << "\t" << 1 << "\t" << match_v.size()
              << std::endl;
  }
};

struct pinfo_t // TODO: May simplify, might be storing too much information.
{
  uint32_t infhdist = std::numeric_limits<uint32_t>::max();
  uint32_t suphdist = std::numeric_limits<uint32_t>::min();
  ;
  float avghdist = std::numeric_limits<float>::max();
  uint32_t match_count = 0;
};

struct ninfo_t
{
  std::vector<minfo_sptr_t> minfo_v;
  flat_phmap<uint32_t, pinfo_t>
    pos_to_pinfo; // TODO: You could potentially keep matches as they are instead of positional info and compute on demand.
  float score = 0;     // TODO: This is a score but which? Does this have to be float?
  float pavghdist = 0; // TODO: Remove?
  uint32_t taxa_count = 0;
  uint32_t match_count = 0;
  float covmer = 0.0;
  float covpos = 0.0;
  float avghdist = std::numeric_limits<float>::max();
  float exchdist = std::numeric_limits<float>::max();
  ninfo_t() {}
  ninfo_t(minfo_sptr_t mi)
  {
    taxa_count = 1;
    match_count = mi->match_count;
    avghdist = mi->avghdist;
    covmer = mi->covmer;
    covpos = mi->covpos;
    minfo_v.push_back(mi);
    map_pinfo(mi);
  }
  void add_minfo(minfo_sptr_t mi)
  {
    taxa_count++;
    avghdist = (avghdist * match_count + (mi->avghdist * mi->match_count));
    match_count += mi->match_count;
    avghdist /= match_count;
    covmer = std::max(mi->covmer, covmer);
    covpos = std::max(mi->covpos, covpos);
    minfo_v.push_back(mi);
    map_pinfo(mi);
  }
  void map_pinfo(minfo_sptr_t mi)
  {
    for (auto const& match : mi->match_v) {
      pinfo_t& pi = pos_to_pinfo[match.pos];
      pi.avghdist = (pi.avghdist * pi.match_count) + match.hdist;
      pi.match_count++;
      pi.avghdist /= pi.match_count;
      pi.infhdist = std::min(match.hdist, pi.infhdist);
      pi.suphdist = std::max(match.hdist, pi.suphdist);
    }
  }
  void print_info()
  {
    std::cout << covpos << "\t" << covmer << "\t" << avghdist << "\t" << taxa_count << "\t"
              << match_count << std::endl;
  }
};

class QMers
{
  friend class QBatch;

public:
  QMers(library_sptr_t library, uint64_t len, uint32_t max_hdist = 3, float min_covpos = 0.5);
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void print_matches();
  void compute_coverage();
  void print_coverage();
  void fill_ninfo();
  void compute_exchdist();
  float greedy_count_placement();
  float argmin_avghdist_placement();
  float argmin_diffhdist_placement();
  float countpos_avghdist_placement(); // TODO: Perhaps better naming;
  void display_placement();

private:
  uint8_t k;
  uint32_t len;
  uint32_t max_hdist;
  float min_covpos;
  tree_sptr_t tree;
  library_sptr_t library;
  crecord_sptr_t crecord;
  node_sptr_t placement;
  flat_phmap<node_sptr_t, minfo_sptr_t> node_to_minfo;
  flat_phmap<node_sptr_t, ninfo_sptr_t> node_to_ninfo;
};

class QBatch
{
public:
  QBatch(library_sptr_t library, qseq_sptr_t qs);
  void search_batch(uint32_t max_hdist, float min_covpos);
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
