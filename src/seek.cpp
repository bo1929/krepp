#include "seek.hpp"
#include <boost/math/tools/minima.hpp>

#define DD 0.1
#define KK 1000

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
  mp = llhfunc.mp(DD, 0.2);
}

void SBatch::seek_sequences(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    enmers = len - k + 1;
    onmers = 0;

    SSummary or_summary(enmers, hdist_th);
    SSummary rc_summary(enmers, hdist_th);
    search_mers(seq, len, or_summary, rc_summary);

    // if (or_summary.match_count + rc_summary.match_count) {
    //   or_summary.mismatch_count = onmers - or_summary.match_count;
    //   rc_summary.mismatch_count = onmers - rc_summary.match_count;
    //   or_summary.optimize_likelihood(llhfunc, rho);
    //   rc_summary.optimize_likelihood(llhfunc, rho);
    //   if (or_summary.d_llh < rc_summary.d_llh) {
    //     batch_stream << identifer_batch[bix] << "\t" << or_summary.d_llh << "\n";
    //   } else {
    //     batch_stream << identifer_batch[bix] << "\t" << rc_summary.d_llh << "\n";
    //   }
    // } else {
    //   batch_stream << identifer_batch[bix] << "\tNaN\n";
    // }
    // for (uint32_t i = 0; i < x_v.size(); ++i) {
    // #pragma omp critical
    //   // std::cout << x_v[i] << std::endl;
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
