#include "seek.hpp"
#include <boost/math/tools/minima.hpp>
#include <iostream>
#include <vector>
#include <stack>
#include <climits>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <limits>
#include <utility>
#include <algorithm> // for std::upper_bound

using std::make_pair;
using std::make_tuple;
using std::max;
using std::min;
using std::pair;
using std::tuple;
using std::vector;

#define DD 0.1055
#define KK 10000

vector<long double> prefix_sum(const vector<double>& X)
{
  int n = X.size();
  vector<long double> P(n + 1, 0.0);
  for (int i = 0; i < n; ++i) {
    P[i + 1] = P[i] + (-X[i]);
    // std::cout << P[i + 1] << std::endl;
  }
  // for (int i = 0; i < n; ++i) {
  //   P[i] = std::lround(P[i] / 2) - std::lround(P[i] / 2) % 2;
  //   // std::cout << P[i + 1] << std::endl;
  // }
  return P;
}

vector<long double> prefix_maxima(const vector<long double>& P)
{
  int n = P.size();
  vector<long double> M(n);
  long double running_max = -std::numeric_limits<long double>::infinity();
  for (int i = 0; i < n; ++i) {
    running_max = max(running_max, P[i]);
    M[i] = running_max;
  }
  return M;
}

vector<long double> suffix_min(const vector<long double>& P)
{
  int n = P.size();
  vector<long double> S_min(n + 1, std::numeric_limits<long double>::infinity());
  for (int i = n - 1; i >= 0; --i)
    S_min[i] = min(P[i], S_min[i + 1]);
  return S_min;
}

pair<vector<tuple<int, int, int>>, vector<int>> find_monotone_segments(const vector<long double>& P)
{
  int n = P.size();
  vector<tuple<int, int, int>> segments;
  vector<int> labels(n, 0);
  int start = 0, dirn = 0;
  for (int i = 1; i < n; ++i) {
    int sgn = (P[i] > P[i - 1]) - (P[i] < P[i - 1]);
    if (dirn == 0 && sgn != 0)
      dirn = sgn;
    else if (sgn != 0 && sgn != dirn) {
      // std::cout << "SEGMENT SIZE: " << i - start + 1 << std::endl;
      segments.push_back(make_tuple(start, i - 1, dirn));
      start = i - 1;
      dirn = sgn;
    }
    labels[i] = segments.size();
  }
  segments.push_back(make_tuple(start, n - 1, dirn));
  return {segments, labels};
}

// vector<int> compute_jump_labels(const std::vector<int>& pm,
//                                 const std::vector<int>& labels,
//                                 const vector<tuple<int, int, int>>& segments,
//                                 const std::vector<long double>& P)
// {
//   int S = segments.size();
//   std::vector<double> keys(S);
//   for (int i = 0; i < S; ++i) {
//     int a, b, d;
//     std::tie(a, b, d) = segments[i];
//     keys[i] = (d == 1) ? a : b;
//   }

//   std::vector<int> jl;
//   for (int a : pm) {
//     int seg_label = labels[a];

//     auto it = std::upper_bound(keys.begin() + seg_label, keys.end(), P[a]);
//     if (it == keys.begin() + seg_label) {
//       jl.push_back(seg_label);
//     } else {
//       jl.push_back(std::distance(keys.begin(), it) - 1);
//     }
//   }
//   return jl;
// }

vector<int> compute_jump_labels(const vector<int>& pm,
                                const vector<int>& labels,
                                const vector<tuple<int, int, int>>& segments,
                                const vector<long double>& P)
{
  vector<int> jl;
  int S = segments.size();
  // std::cout << "NUM SEGMENTS=" << S << std::endl;
  for (int a : pm) {
    int left = labels[a], right = S - 1, best_idx = left;
    while (left <= right) {
      int mid = (left + right) / 2;
      int i, j, d;
      std::tie(i, j, d) = segments[mid];
      long double Ps = P[i], Pe = P[j];
      if (d == 1) {
        if (Ps <= P[a]) {
          best_idx = mid;
          left = mid + 1;
        } else
          right = mid - 1;
      } else {
        if (Pe <= P[a]) {
          best_idx = mid;
          left = mid + 1;
        } else
          right = mid - 1;
      }
    }
    jl.push_back(best_idx);
  }
  return jl;
}

vector<vector<pair<long double, int>>>
build_segment_arrays(const vector<tuple<int, int, int>>& segments, const vector<long double>& P)
{
  vector<vector<pair<long double, int>>> segment_data;
  for (auto& t : segments) {
    int s, e, d;
    std::tie(s, e, d) = t;
    if (d == 1) {
      vector<pair<long double, int>> vals;
      for (int i = s; i <= e; ++i)
        vals.push_back(make_pair(P[i], i));
      segment_data.push_back(vals);
    } else {
      pair<long double, int> min_val = make_pair(P[s], s);
      for (int i = s; i <= e; ++i)
        if (P[i] < min_val.first) min_val = make_pair(P[i], i);
      segment_data.push_back({min_val});
    }
  }
  return segment_data;
}

// int binary_search_in_segment(const vector<long double>& P,
//                              int a,
//                              int K,
//                              int seg_start,
//                              int seg_end,
//                              const vector<pair<long double, int>>& seg_vals)
// {
//   int lo = max(a + K, seg_start);
//   int hi = seg_end;
//   if (seg_start > hi) return -1;
//   if (seg_vals.size() == 1 && seg_vals[0].second >= lo) {
//     if (seg_vals[0].first <= P[a] && seg_vals[0].second >= lo)
//       return seg_vals[0].second;
//     else
//       return -1;
//   }
//   int l = 0, r = seg_vals.size() - 1, pos = -1;
//   while (l <= r) {
//     int m = (l + r) / 2;
//     if (seg_vals[m].first <= P[a]) {
//       pos = m;
//       l = m + 1;
//     } else
//       r = m - 1;
//   }
//   if (pos == -1) return -1;
//   while (pos >= 0 && seg_vals[pos].second < lo)
//     pos--;
//   return pos >= 0 ? seg_vals[pos].second : -1;
// }

int binary_search_in_segment(const vector<long double>& P,
                             int a,
                             int K,
                             int seg_start,
                             int seg_end,
                             const vector<pair<long double, int>>& seg_vals)
{
  int lo = std::max(a + K, seg_start);
  if (seg_start > seg_end) return -1;

  // Decreasing or constant segment (single min value)
  if (seg_vals.size() == 1) {
    int idx = seg_vals[0].second;
    return (seg_vals[0].first <= P[a] && idx >= lo) ? idx : -1;
  }

  // Increasing segment: use upper_bound on value
  auto it = std::upper_bound(seg_vals.begin(), seg_vals.end(), make_pair(P[a], std::numeric_limits<int>::max()));
  if (it == seg_vals.begin()) return -1; // no value <= P[a]

  --it; // last element <= P[a]
  while (it >= seg_vals.begin() && it->second < lo)
    --it; // ensure index >= lo
  return (it >= seg_vals.begin() && it->second >= lo) ? it->second : -1;
}

vector<pair<int, int>> monotonic_segment_optimized_intervals(const vector<double>& X, int K)
{
  int n = X.size();
  vector<long double> P = prefix_sum(X);
  vector<long double> M = prefix_maxima(P);
  vector<long double> S_min = suffix_min(P);
  auto [segments, labels] = find_monotone_segments(P);
  vector<vector<pair<long double, int>>> segment_data = build_segment_arrays(segments, P);
  vector<int> pm;
  long double cur_max = -std::numeric_limits<long double>::infinity();
  for (int i = 0; i < P.size(); ++i)
    if (P[i] >= cur_max) {
      pm.push_back(i);
      cur_max = P[i];
    }
  vector<int> jl = compute_jump_labels(pm, labels, segments, P);
  vector<pair<int, int>> intervals;
  for (int idx = 0; idx < pm.size(); ++idx) {
    int a = pm[idx];
    if (a >= n) continue;
    int seg_idx = jl[idx];
    int seg_start, seg_end, d;
    std::tie(seg_start, seg_end, d) = segments[seg_idx];
    int b = binary_search_in_segment(P, a, K, seg_start, seg_end, segment_data[seg_idx]);
    if (b == -1) continue;
    if (a > 0 && P[b] <= M[a - 1]) continue;
    if (b < n && S_min[b + 1] <= P[a]) continue;
    if (P[a] >= P[b]) intervals.push_back(make_pair(a, b));
  }
  for (auto& p : intervals) {
    std::cout << p.first << "," << p.second << "\t" << P[p.first] - P[p.second + 1] << "\t" << P[p.first] - P[p.second]
              << "\n";
  }
  return intervals;
}

SBatch::SBatch(sketch_sptr_t sketch, qseq_sptr_t qs, uint32_t hdist_th)
  : sketch(sketch)
  , hdist_th(hdist_th)
{
  lshf = sketch->get_lshf();
  k = lshf->get_k();
  h = lshf->get_h();
  m = lshf->get_m();
  batch_size = qs->cbatch_size;
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->identifer_batch, identifer_batch);
  llhfunc = optimize::HDistHistLLH(h, k, hdist_th);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  rho = sketch->get_rho();
  mp = llhfunc.mp(DD, rho);
}

void SBatch::seek_sequences(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    x_v.clear();
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    enmers = len - k + 1;
    onmers = 0;

    SSummary or_summary(enmers, hdist_th);
    SSummary rc_summary(enmers, hdist_th);
    search_mers(seq, len, or_summary, rc_summary);

    vector<pair<int, int>> res = monotonic_segment_optimized_intervals(x_v, KK);
    //for (uint32_t i = 0; i < x_v.size(); ++i) {
    // #pragma omp critical
    //   std::cout << x_v[i] << std::endl;
    // }
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void SBatch::search_mers(const char* seq, uint64_t len, SSummary& or_summary, SSummary& rc_summary)
{
  uint32_t i, l;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (i = l = 0; i < len;) {
    if (seq_nt4_table[seq[i]] >= 4) {
      l = 0, i++;
      continue;
    }
    l++, i++;
    if (l < k) {
      continue;
    }
    if (l == k) {
      compute_encoding(seq + i - k, seq + i, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(seq + i - 1, orenc64_lr, orenc64_bp);
    }
    orenc64_bp = orenc64_bp & mask_bp;
    orenc64_lr = orenc64_lr & mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    onmers++; // TODO: Incorporate missing partial/fraction?
    uint32_t hdist_curr = std::numeric_limits<uint32_t>::max();
    uint32_t hdist_or = hdist_curr, hdist_rc = hdist_curr;
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (sketch->check_partial(orrix)) {
        hdist_or = or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (sketch->check_partial(rcrix)) {
        hdist_rc = rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (sketch->check_partial(orrix)) {
      hdist_or = or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (sketch->check_partial(rcrix)) {
      hdist_rc = rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
    }
#endif /* CANONICAL */
    hdist_curr = std::min(hdist_or, hdist_rc);
    if (hdist_curr <= hdist_th) {
      x_v.push_back((static_cast<double>(hdist_curr) - k * DD) / (DD * (1 - DD)));
    } else {
      x_v.push_back(mp);
    }
  }
}

uint32_t SSummary::add_matching_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr)
{
  uint32_t hdist_curr;
  uint32_t hdist_min = hdist_th + 1;
  std::pair<vec_enc_it, vec_enc_it> indices = sketch->bucket_indices(rix);
  for (; indices.first < indices.second; ++indices.first) {
    hdist_curr = popcount_lr32((*indices.first) ^ enc_lr);
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
    }
  }
  if (hdist_min <= hdist_th) {
    mismatch_count--;
    match_count++;
    hdisthist_v[hdist_min]++;
  }
  return hdist_min;
}
void SSummary::optimize_likelihood(optimize::HDistHistLLH llhfunc, double rho)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
  d_llh = sol_r.first;
  v_llh = sol_r.second;
}
