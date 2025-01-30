#include "compare.hpp"

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
  or_hdisthist_v.resize(hdist_th + 1, 0);
  rc_hdisthist_v.resize(hdist_th + 1, 0);
  rho = sketch->get_rho();
}

void CBatch::estimate_distances(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    enmers = len - k + 1;

    or_match_count = 0;
    or_mismatch_count = enmers;
    rc_match_count = 0;
    rc_mismatch_count = enmers;
    std::fill(or_hdisthist_v.begin(), or_hdisthist_v.end(), 0);
    std::fill(rc_hdisthist_v.begin(), rc_hdisthist_v.end(), 0);
    search_mers(seq, len);

    if (or_match_count + rc_match_count) {
      llhfunc.set_parameters(or_hdisthist_v.data(), or_mismatch_count, rho);
      std::pair<double, double> sol_r;
      sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
      or_d_llh = sol_r.first;
      or_v_llh = sol_r.second;
      llhfunc.set_parameters(rc_hdisthist_v.data(), rc_mismatch_count, rho);
      sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
      rc_d_llh = sol_r.first;
      rc_v_llh = sol_r.second;
      if (or_d_llh < rc_d_llh) {
        batch_stream << identifer_batch[bix] << "\t" << or_d_llh << "\n";
      } else {
        batch_stream << identifer_batch[bix] << "\t" << rc_d_llh << "\n";
      }
    } else {
      batch_stream << identifer_batch[bix] << "\tNaN\n";
    }
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void CBatch::search_mers(const char* seq, uint64_t len)
{
  uint32_t i, l;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  uint32_t curr_hdist;
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
        add_or_mer(orrix, lshf->drop_ppos_lr(orenc64_lr));
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (sketch->set_partial(rcrix)) {
        add_rc_mer(rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (sketch->set_partial(orrix)) {
      add_or_mer(orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (sketch->set_partial(rcrix)) {
      add_rc_mer(rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
    }
#endif /* CANONICAL */
  }
}
