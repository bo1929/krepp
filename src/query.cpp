#include "query.hpp"
#include "common.hpp"
#include "hdhistllh.hpp"

QMers::QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th, double min_gamma)
  : library(library)
  , len(len)
  , nmers(0)
  , hdist_th(hdist_th)
  , min_gamma(min_gamma)
{
  lshf = library->get_lshf();
  tree = library->get_tree();
  k = lshf->get_k();
  h = lshf->get_h();
}

void QBatch::search_batch(uint32_t hdist_th, double min_gamma)
{
#pragma omp parallel for num_threads(num_threads)
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, hdist_th, min_gamma);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, hdist_th, min_gamma);

    search_mers(seq, len, qmers_or, qmers_rc);

    qmers_or->summarize_minfo();
    qmers_rc->summarize_minfo();

#pragma omp critical
    {
      add_to_report(qmers_or, qmers_rc, bix);
    }
  }
  std::cout << report;
}

void Minfo::compute_gamma()
{
  uint32_t s;
  uint32_t i, j, k;
  for (i = 0; i < match_v.size(); ++i) {
    if (i == 0) {
      gamma = qmers->k;
      continue;
    }
    if (match_v[i].pos > match_v[i - 1].pos) {
      s = (match_v[i].pos - match_v[i - 1].pos);
    } else {
      s = (match_v[i - 1].pos - match_v[i].pos);
    }
    if (s > qmers->k) {
      gamma += qmers->k;
    } else {
      gamma += s;
    }
  }
  gamma /= qmers->len;
}

void Minfo::estimate_distance(optimize::HDistHistLLH& llhfunc)
{
  llhfunc.set_mc(hdisthist_v.data());
  llhfunc.set_ro(rho);
  optimize::Lbfgsb solver;
  optimize::State state =
    solver.minimize(llhfunc, optimize::d_init, optimize::d_lb, optimize::d_ub);
  d_llh = (state.x().transpose())(0);
}

void QMers::summarize_minfo()
{
  optimize::HDistHistLLH llhfunc(h, k, hdist_th, nmers);
  for (auto& [nd, mi] : node_to_minfo) {
    /* mi->compute_gamma(); */
    mi->gamma = static_cast<double>(mi->match_count) / static_cast<double>(len - k + 1);
    mi->estimate_distance(llhfunc);
  }
}

QBatch::QBatch(library_sptr_t library, qseq_sptr_t qs)
  : library(library)
{
  lshf = library->get_lshf();
  tree = library->get_tree();
  k = lshf->get_k();
  m = lshf->get_m();
  batch_size = qs->batch_size;
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->name_batch, name_batch);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
}

void QBatch::add_to_report(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint32_t bix)
{
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    if (mi->gamma > qmers_or->min_gamma &&
        !(qmers_rc->node_to_minfo.contains(nd) && mi->d_llh > qmers_rc->node_to_minfo[nd]->d_llh)) {
      report += name_batch[bix] + "\t" + "or" + "\t" + std::to_string(mi->match_count) + "\t" +
                std::to_string(qmers_or->nmers) + "\t" + std::to_string(mi->rho) + "\t" +
                nd->get_name() + "\t" + std::to_string(mi->d_llh) + "\t" +
                std::to_string(mi->gamma);
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        report += "\t" + std::to_string(mi->hdisthist_v[i]);
      }
      report += "\n";
    }
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->gamma > qmers_rc->min_gamma &&
        !(qmers_or->node_to_minfo.contains(nd) && mi->d_llh > qmers_or->node_to_minfo[nd]->d_llh)) {
      report += name_batch[bix] + "\t" + "rc" + "\t" + std::to_string(mi->match_count) + "\t" +
                std::to_string(qmers_or->nmers) + "\t" + std::to_string(mi->rho) + "\t" +
                nd->get_name() + "\t" + std::to_string(mi->d_llh) + "\t" +
                std::to_string(mi->gamma);
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        report += "\t" + std::to_string(mi->hdisthist_v[i]);
      }
      report += "\n";
    }
  }
}

void QBatch::search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
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
      if (library->check_partial(orrix)) {
        qmers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (library->check_partial(rcrix)) {
        qmers_or->add_matching_mer(i - k, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (library->check_partial(orrix)) {
      qmers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (library->check_partial(rcrix)) {
      qmers_rc->add_matching_mer(len - i, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
    }
#endif /* CANONICAL */
  }
}

void QMers::add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  se_t se;
  node_sptr_t nd;
  uint32_t ix;
  uint32_t curr_hdist;
  std::queue<se_t> subset_q;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = library->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->get_next(rix);
  crecord_sptr_t crecord = library->get_crecord(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist > hdist_th) {
      continue;
    }
    subset_q.push(iter1->second);
    while (!subset_q.empty()) {
      se = subset_q.front();
      subset_q.pop();
      if (!crecord->check_node(se)) {
        pse = crecord->get_pse(se);
        subset_q.push(pse.first);
        subset_q.push(pse.second);
        continue;
      }
      nd = crecord->get_node(se);
      if (nd->check_leaf() && (nd->get_name() != leave_out_ref)) { // TODO: Remove testing.
        if (!node_to_minfo.contains(nd)) {
          node_to_minfo[nd] = std::make_unique<Minfo>(getptr(), crecord->get_rho(se));
        }
        node_to_minfo[nd]->update_match(iter1->first, pos, curr_hdist);
      } else { // TODO: This might not be needed.
        for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
          subset_q.push((*std::next(nd->get_children(), i))->get_se());
        }
      }
    }
  }
  nmers++;
}
