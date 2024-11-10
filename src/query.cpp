#include "query.hpp"
#include <boost/math/tools/minima.hpp>
#include <string>

#define CHISQ_THRESHOLD 2.706

QBatch::QBatch(index_sptr_t index, qseq_sptr_t qs, uint32_t hdist_th)
  : index(index)
  , hdist_th(hdist_th)
{
  lshf = index->get_lshf();
  tree = index->get_tree();
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
}

void QBatch::search_mers(const char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  enmers = len - k + 1;
  onmers_or = 0;
  onmers_rc = 0;
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
      if (index->check_partial(orrix)) {
        qmers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
        onmers_or++;
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (index->check_partial(rcrix)) {
        qmers_or->add_matching_mer(i - k, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
        onmers_rc++;
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (index->check_partial(orrix)) {
      qmers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
      onmers_or++;
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (index->check_partial(rcrix)) {
      qmers_rc->add_matching_mer(len - i, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      onmers_rc++;
    }
#endif /* CANONICAL */
  }
}

void QBatch::summarize_minfo(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  nd_pp = tree->get_root();
  mi_pp = std::make_shared<Minfo>(hdist_th);
  node_to_minfo.clear();
  qmers_or->hdist_filt = 2 * qmers_or->hdist_filt + 1;
  qmers_rc->hdist_filt = 2 * qmers_rc->hdist_filt + 1;
  for (auto [nd, mi] : qmers_or->leaf_to_minfo) {
    // mi->compute_gamma();
    if (mi->hdist_min > qmers_or->hdist_filt) {
      continue;
    }
    mi->optimize_likelihood(llhfunc);
    if (mi->d_llh < mi_pp->d_llh) {
      nd_pp = nd;
      mi_pp = mi;
    }
    node_to_minfo[nd] = mi;
  }
  for (auto [nd, mi] : qmers_rc->leaf_to_minfo) {
    // mi->compute_gamma();
    if (mi->hdist_min > qmers_rc->hdist_filt) {
      continue;
    }
    mi->optimize_likelihood(llhfunc);
    if (mi->d_llh < mi_pp->d_llh) {
      nd_pp = nd;
      mi_pp = mi;
    }
    node_to_minfo[nd] = mi;
    // If both in reverse-complement and the original sequence, decide:
    if ((qmers_or->leaf_to_minfo).contains(nd) &&
        mi->d_llh > (qmers_or->leaf_to_minfo)[nd]->d_llh) {
      node_to_minfo[nd] = (qmers_or->leaf_to_minfo)[nd];
    }
  }
}

void QBatch::estimate_distances()
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(index, len, hdist_th);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(index, len, hdist_th);

    search_mers(seq, len, qmers_or, qmers_rc);
    summarize_minfo(qmers_or, qmers_rc);
    report_distances(batch_stream);
  }
#pragma omp critical
  std::cout << batch_stream.rdbuf();
}

void QBatch::report_distances(strstream& batch_stream)
{
  for (auto& [nd, mi] : node_to_minfo) {
    // TODO: Add distance threshold for reporting or likelihood ratio test.
    batch_stream << identifer_batch[bix] << "\t" << nd->get_name() << "\t" << mi->d_llh << "\n";
  }
}

void QBatch::place_sequences()
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(index, len, hdist_th);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(index, len, hdist_th);

    search_mers(seq, len, qmers_or, qmers_rc);
    summarize_minfo(qmers_or, qmers_rc);
    report_placement(batch_stream);
  }
#pragma omp critical
  std::cout << batch_stream.rdbuf();
}

void QBatch::report_placement(strstream& batch_stream)
{
  if (node_to_minfo.size() == 0) {
    batch_stream << "\t\t\t{\"p\" : [ ], \"n \": [\"" + identifer_batch[bix] + "\"]},\n";
    return;
  } else if (node_to_minfo.size() > 1) { // OR place on the closest leaf.
    node_sptr_t nd_curr = nullptr;
    minfo_sptr_t mi_curr = nullptr;

    /*
    // Place on the LCA of all leaves that pass llh-ratio test:
    mi_curr = std::make_shared<Minfo>(hdist_th);
    for (auto const& [nd, mi] : node_to_minfo) {
      mi->chisq = mi->likelihood_ratio(mi_pp->d_llh, llhfunc);
      if (mi->chisq < CHISQ_THRESHOLD) {
        nd_curr = Tree::compute_lca(nd_curr, nd);
        mi_curr.add(mi);
      }
    }
    mi_curr->chisq = mi_curr->likelihood_ratio(mi_pp->d_llh, llhfunc);
    nd_pp = nd_curr;
    mi_pp = mi_curr;
    */
    vec<node_sptr_t> nd_v;
    while (nd_curr = tree->next_post_order(nd_curr)) {
      if (nd_curr->check_leaf()) {
        if (!node_to_minfo.contains(nd_curr)) {
          node_to_minfo[nd_curr] = std::make_shared<Minfo>(enmers, hdist_th);
        } else {
          mi_curr = node_to_minfo[nd_curr];
          mi_curr->chisq = mi_curr->likelihood_ratio(mi_pp->d_llh, llhfunc);
          if (mi_curr->chisq < CHISQ_THRESHOLD) {
            nd_v.push_back(nd_curr);
          }
        }
      } else {
        mi_curr = std::make_shared<Minfo>(hdist_th);
        for (uint32_t nix = 0; nix < nd_curr->get_nchildren(); ++nix) {
          mi_curr->add(node_to_minfo[*std::next(nd_curr->get_children(), nix)]);
        }
        if (mi_curr->rmatch_count) {
          mi_curr->optimize_likelihood(llhfunc);
          mi_curr->chisq = mi_curr->likelihood_ratio(mi_pp->d_llh, llhfunc);
          if (mi_curr->chisq < CHISQ_THRESHOLD) {
            mi_curr->optimize_likelihood(llhfunc);
            nd_v.push_back(nd_curr);
          }
        }
        node_to_minfo[nd_curr] = mi_curr;
      }
    }

    std::sort(nd_v.begin(), nd_v.end(), [&](node_sptr_t lhs, node_sptr_t rhs) {
      if (lhs->get_card() == rhs->get_card())
        return node_to_minfo[lhs]->d_llh > node_to_minfo[rhs]->d_llh;
      else
        return lhs->get_card() < rhs->get_card();
    });
    nd_pp = *(nd_v.rbegin());
    mi_pp = node_to_minfo[nd_pp];
  } else {
    mi_pp->chisq = 0;
  }
  batch_stream << "\t\t\t{\"p\" : [[" + std::to_string(nd_pp->get_se() - 1) + ", " +
                    std::to_string(mi_pp->chisq) + ", " +
                    /* std::to_string(mi_pp->d_llh) + ", " + */
                    std::to_string(mi_pp->v_llh) + ", " + std::to_string(mi_pp->d_llh) + ", " +
                    std::to_string(1e-5) + ", " +
                    /* std::to_string(mi_pp->rmatch_count) + ", " + */
                    std::to_string(nd_pp->get_blen() / 2.0) + "]], \"n\" : [\"" +
                    identifer_batch[bix] + "\"]},\n";
}

QMers::QMers(index_sptr_t index, uint64_t len, uint32_t hdist_th)
  : index(index)
  , len(len)
  , hdist_th(hdist_th)
  , onmers(0)
{
  lshf = index->get_lshf();
  tree = index->get_tree();
  k = lshf->get_k();
  h = lshf->get_h();
  enmers = len - k + 1;
}

void QMers::add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{ // TODO: Potential acceleration using SIMD?
  se_t se;
  node_sptr_t nd;
  uint32_t ix;
  uint32_t hdist_curr;
  std::queue<se_t> se_q;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = index->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = index->get_next(rix);
  crecord_sptr_t crecord = index->get_crecord(rix);
  for (; iter1 < iter2; ++iter1) {
    hdist_curr = popcount_lr32(iter1->first ^ enc_lr);
    if (hdist_curr > hdist_th) {
      continue;
    }
    if (hdist_curr < hdist_filt) {
      hdist_filt = hdist_curr;
    }
    se_q.push(iter1->second);
    while (!se_q.empty()) {
      se = se_q.front();
      se_q.pop();
      if (!tree->check_node(se)) {
        pse = crecord->get_pse(se);
        se_q.push(pse.first);
        se_q.push(pse.second);
        continue;
      }
      nd = tree->get_node(se);
      if (nd->check_leaf()) {
        if (!leaf_to_minfo.contains(nd)) {
          leaf_to_minfo[nd] = std::make_shared<Minfo>(enmers, hdist_th, crecord->get_rho(se));
        }
        leaf_to_minfo[nd]->update_match(iter1->first, pos, hdist_curr);
      } else {
        for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
          se_q.push((*std::next(nd->get_children(), i))->get_se());
        }
      }
    }
  }
  onmers++;
}

/* void Minfo::compute_gamma() */
/* { */
/*   // Alternative 1: number of k-mers covered. */
/*   gamma = static_cast<double> match_count / static_cast<double>(nmers); */
/*   // Alternative 2: number of positions covered. */
/*   // Requires matches to be collected. */
/*   uint32_t s; */
/*   uint32_t i, j, k; */
/*   uint32_t ugamma = 0; */
/*   for (i = 0; i < match_v.size(); ++i) { */
/*     if (i == 0) { */
/*       ugamma = qmers->k; */
/*       continue; */
/*     } */
/*     if (match_v[i].pos > match_v[i - 1].pos) { */
/*       s = (match_v[i].pos - match_v[i - 1].pos); */
/*     } else { */
/*       s = (match_v[i - 1].pos - match_v[i].pos); */
/*     } */
/*     if (s > qmers->k) { */
/*       ugamma += qmers->k; */
/*     } else { */
/*       ugamma += s; */
/*     } */
/*   } */
/*   gamma = ugamma / (nmers + qmers->k - 1); */
/* } */

double Minfo::likelihood_ratio(double d, optimize::HDistHistLLH& llhfunc)
{ // TODO: What about testing this the other way around?
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  return 2 * (llhfunc(d) - v_llh);
}

void Minfo::optimize_likelihood(optimize::HDistHistLLH& llhfunc)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);

  // Locating Function Minima using Brent's algorithm, depends on boost::math.
  std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
  d_llh = sol_r.first;
  v_llh = sol_r.second;

  // Optimizing with L-BFGS-B, requires Eigen and l-bfgs-b libraries.
  /* optimize::Lbfgsb solver; */
  /* optimize::State state = */
  /*   solver.minimize(llhfunc, optimize::d_init, optimize::d_lb, optimize::d_ub); */
  /* d_llh = (state.x().transpose())(0); */
}
