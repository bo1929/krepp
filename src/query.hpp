#ifndef _QUERY_H
#define _QUERY_H

#include "common.hpp"
#include "index.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "table.hpp"
#include "hdhistllh.hpp"

#define CJC 4.0 / 3.0

inline double chisq_cdf(double chisq)
{
  constexpr double kStep = 0.05;
  // clang-format off
  static const double cdf_v[] = {
      0, 0.039877611676744917, 0.079655674554057948, 0.11923538474048505,
      0.15851941887820606, 0.19741265136584743, 0.23582284437790529, 0.27366130235123814,
      0.31084348322064831, 0.34728955942415984, 0.38292492254802613, 0.41768062642330722,
      0.45149376449985279, 0.48430777838827049, 0.51607269555385404, 0.54674529524626359,
      0.57628920283320673, 0.60467491375461524, 0.63187974930648094, 0.65788774738303646,
      0.68268949213708585, 0.70628188724820817, 0.72866787810723466, 0.74985612872569951,
      0.76986065955658356, 0.78870045266628941, 0.8063990308287794, 0.82298401712519587,
      0.83848668153245798, 0.85294148078070331, 0.86638559746228383, 0.87885848399588196,
      0.89040141660088401, 0.90105706393270379, 0.91086907448291399, 0.91988168627236577,
      0.92813936177414835, 0.93568645040877263, 0.94256688036799641, 0.94882388095672276,
      0.95449973610364158, 0.95963556918859116, 0.96427115887436687, 0.96844478521781896,
      0.97219310497300282, 0.97555105468991066, 0.97855177995664833, 0.98122658893032288,
      0.98360492815080769, 0.98571437852945709, 0.98758066934844768, 0.98922770809186655,
      0.99067762395256254, 0.99195082291448333, 0.99306605239391865, 0.99404047352989089,
      0.99488973933914415, 0.99562807709017354, 0.99626837339923191, 0.99682226070527025,
      0.99730020393673979, 0.99771158633795465, 0.99806479357356337, 0.99836729537434288,
      0.99862572412416828, 0.99884594991521847, 0.99903315171523244, 0.99919188439627193,
      0.99932614146864629, 0.99943941344636766, 0.99953474184192892, 0.99961476884872869,
      0.99968178281968489, 0.99973775969115897, 0.99978440053304518, 0.99982316542959837,
      0.99985530391214983, 0.99988188217516216, 0.99990380731196482, 0.99992184880680446,
      0.99993665751633376, 0.99994878236705187, 0.99995868498617491, 0.99996675247254063,
      0.99997330850196819, 0.99997862294845019, 0.99998292018905799, 0.99998638624680125,
      0.99998917491218453, 0.99999141297106009, 0.99999320465375052, 0.9999946354084408,
      0.99999577509059501, 0.99999668064971126, 0.99999739838509216, 0.99999796583351486,
      0.9999984133436961, 0.99999876538525601, 0.99999904163344677, 0.99999925786518407,
      0.99999942669685626, 0.99999955818993547, 0.99999966034651855, 0.99999973951354093,
      0.99999980071147365, 0.99999984790078966, 0.99999988419731922, 0.9999999120457681,
      0.99999993335910298, 0.99999994963017991, 0.99999996202087504, 0.9999999714330402,
      0.99999997856481948, 0.99999998395521628, 0.99999998801925716, 0.99999999107565507,
      0.99999999336850798, 0.99999999508426995, 0.99999999636498427, 0.99999999731857514,
      0.99999999802682471, 0.99999999855154165, 0.99999999893931535, 0.99999999922517047,
      0.99999999943536833, 0.99999999958954733, 0.99999999970235431, 0.99999999978468512,
      0.99999999984462307, 0.99999999988814992, 0.99999999991968003, 0.99999999994246291,
      0.99999999995888422, 0.99999999997069078, 0.99999999997915801, 0.99999999998521549,
      0.99999999998953815, 0.99999999999261502, 0.99999999999479972, 0.99999999999634714,
      0.99999999999744038, 0.99999999999821076, 0.99999999999875244, 0.99999999999913225,
      0.99999999999939782, 0.99999999999958322, 0.99999999999971223, 0.99999999999980171,
      0.99999999999986389, 0.99999999999990674, 0.99999999999993616, 0.99999999999995648,
      0.99999999999997047, 0.9999999999999799, 0.99999999999998646, 0.99999999999999079,
      0.99999999999999378, 0.99999999999999578, 0.99999999999999722, 0.99999999999999811,
      0.99999999999999878, 0.99999999999999911, 0.99999999999999944, 0.99999999999999956,
      0.99999999999999978, 0.99999999999999978, 0.99999999999999989, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      1
  };
  // clang-format on
  if (std::isnan(chisq)) {
    return chisq;
  }
  if (chisq <= 0.0) {
    return 0.0;
  }
  const double pos = std::sqrt(chisq) / kStep;
  constexpr size_t n = sizeof(cdf_v) / sizeof(cdf_v[0]);
  if (pos > static_cast<double>(n - 1)) {
    return 1.0;
  }
  size_t i = static_cast<size_t>(pos);
  if (i + 1 >= n) {
    i = n - 2;
  }
  return cdf_v[i] + (cdf_v[i + 1] - cdf_v[i]) * (pos - i);
}

namespace optimize {
  class HDistHistLLH;
}

class Minfo;
class IMers;
typedef std::shared_ptr<Minfo> minfo_sptr_t;
typedef std::unique_ptr<Minfo> minfo_uptr_t;
typedef std::shared_ptr<IMers> imers_sptr_t;

// A single placement of one query onto one edge of the tree.
struct placement_t
{
  se_t edge_num = 0;
  double pendant_length = 0;
  double distal_length = 0;
  double likelihood = 0;
  double like_weight_ratio = 0;
  double distance = 0;
  std::string distal_node;
  node_sptr_t node = nullptr;
};

class IMers : public std::enable_shared_from_this<IMers>
{
  friend class IBatch;

public:
  IMers(index_sptr_t index, uint64_t len, uint32_t hdist_th);
  imers_sptr_t getptr() { return shared_from_this(); }
  void add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  inline void
  add_matching_mer_view(uint32_t pos, uint32_t rix, enc_t enc_lr, const FlatHT* flatht, CRecord* crecord, uint32_t numerator);
  bool stage_mer(uint32_t pos, uint32_t rix, enc_t enc_lr);
  void flush_mers();

private:
  /* A queued lookup with its bucket already resolved. */
  struct PendingMer
  {
    uint32_t pos;
    enc_t enc_lr;
    const cmer_t* first;
    const cmer_t* last;
    CRecord* crecord;
  };

  void process_mer(const PendingMer& mer);
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
  vec<uint32_t> vnd_v;
  vec<se_t> se_v;
  uint32_t tix = 0;
  vec<PendingMer> pending_v;
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
  bool collect_placements(vec<placement_t>& placements);
  bool report_placement(strstream& batch_stream, bool tabular, bool has_previous);
  const parallel_flat_phmap<node_sptr_t, double>& get_summary() { return node_to_wcount; }

private:
  placement_t make_placement(const node_sptr_t& nd, const minfo_sptr_t& mi, const minfo_sptr_t& mi_parent);
  static void widen_hdist_filter(uint32_t& hdist_filt);
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

#define PP_JPLACE_FIELDS(pp)                                                                                                \
  "[" << (pp).edge_num << ", " << (pp).pendant_length << ", " << (pp).distal_length << ", " << (pp).likelihood << ", "      \
      << (pp).like_weight_ratio << ", " << (pp).distance << "]"

#define PP_TABULAR_FIELDS(pp)                                                                                               \
  (pp).distal_node << "\t" << (pp).edge_num << "\t" << (pp).like_weight_ratio << "\t" << (pp).distance

#define DISTANCE_FIELDS(nd, mi)                                                                                             \
  nd->get_name() << "\t" << mi->d_llh << "\t" << std::scientific << chisq_cdf(mi->chisq) << std::fixed << "\n"

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
