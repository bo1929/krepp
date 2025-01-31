#include "compare.hpp"
#include <boost/math/tools/minima.hpp>

CBatch::CBatch(sketch_sptr_t sketch, qseq_sptr_t qs, uint32_t hdist_th)
  : sketch(sketch)
  , hdist_th(hdist_th)
{
  lshf = sketch->get_lshf();
  k = lshf->get_k();
  h = lshf->get_h();
  m = lshf->get_m();
  batch_size = qs->batch_size;
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->identifer_batch, identifer_batch);
  llhfunc = optimize::HDistHistLLH(h, k, hdist_th);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  rho = sketch->get_rho();
}

void CBatch::estimate_distances(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    enmers = len - k + 1;

    CSummary or_summary(enmers, hdist_th);
    CSummary rc_summary(enmers, hdist_th);
    search_mers(seq, len, or_summary, rc_summary);

    if (or_summary.match_count + rc_summary.match_count) {
      or_summary.optimize_likelihood(llhfunc, rho);
      rc_summary.optimize_likelihood(llhfunc, rho);
      if (or_summary.d_llh < rc_summary.d_llh) {
        batch_stream << identifer_batch[bix] << "\t" << or_summary.d_llh << "\n";
      } else {
        batch_stream << identifer_batch[bix] << "\t" << rc_summary.d_llh << "\n";
      }
    } else {
      batch_stream << identifer_batch[bix] << "\tNaN\n";
    }
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void CBatch::search_mers(const char* seq, uint64_t len, CSummary& or_summary, CSummary& rc_summary)
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
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (sketch->set_partial(orrix)) {
        or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (sketch->set_partial(rcrix)) {
        rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (sketch->set_partial(orrix)) {
      or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (sketch->set_partial(rcrix)) {
      rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
    }
#endif /* CANONICAL */
  }
}

void CSummary::add_matching_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr)
{
  uint32_t hdist_curr;
  uint32_t hdist_min = hdist_th + 1;
  std::vector<enc_t>::const_iterator iter1 = sketch->bucket_start();
  std::vector<enc_t>::const_iterator iter2 = sketch->bucket_next();
  for (; iter1 < iter2; ++iter1) {
    hdist_curr = popcount_lr32((*iter1) ^ enc_lr);
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
    }
  }
  if (hdist_min <= hdist_th) {
    mismatch_count--;
    match_count++;
    hdisthist_v[hdist_min]++;
  }
}
void CSummary::optimize_likelihood(optimize::HDistHistLLH llhfunc, double rho)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
  d_llh = sol_r.first;
  v_llh = sol_r.second;
}
