#include "query.hpp"
#include <boost/math/tools/minima.hpp>

/* #define CHISQ_THRESHOLD 3.841 */
/* #define CHISQ_THRESHOLD 2.706 */
/* #define CHISQ_THRESHOLD 1.642 */

IBatch::IBatch(index_sptr_t index,
               qseq_sptr_t qs,
               uint32_t hdist_th,
               double chisq_value,
               double dist_max,
               uint32_t tau,
               bool no_filter,
               bool multi,
               bool summarize)
  : index(index)
  , hdist_th(hdist_th)
  , chisq_value(chisq_value)
  , dist_max(dist_max)
  , tau(tau)
  , no_filter(no_filter)
  , multi(multi)
  , summarize(summarize)
{
  lshf = index->get_lshf();
  tree = index->get_tree();
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
}

void IBatch::search_mers(const char* seq, uint64_t len, imers_sptr_t imers_or, imers_sptr_t imers_rc)
{
  enmers = len - k + 1;
  onmers = 0;
  wnmers_or = 0;
  wnmers_rc = 0;
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
    onmers++; // TODO: Incorporate missing fraction/partial?
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (index->check_partial(orrix)) {
        imers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
        wnmers_or++;
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (index->check_partial(rcrix)) {
        imers_rc->add_matching_mer(i - k, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
        wnmers_rc++;
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (index->check_partial(orrix)) {
      imers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
      wnmers_or++;
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (index->check_partial(rcrix)) {
      imers_rc->add_matching_mer(len - i, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)));
      wnmers_rc++;
    }
#endif /* CANONICAL */
  }
}

void IBatch::summarize_matches(imers_sptr_t imers_or, imers_sptr_t imers_rc)
{
  nd_closest = tree->get_root();
  mi_closest = std::make_shared<Minfo>(hdist_th);
  node_to_minfo.clear();
  imers_or->hdist_filt = 2 * imers_or->hdist_filt + 1;
  imers_rc->hdist_filt = 2 * imers_rc->hdist_filt + 1;
  for (auto [nd, mi] : imers_or->leaf_to_minfo) {
    mi->mismatch_count = onmers - mi->match_count;
    // mi->compute_gamma();
    if (mi->hdist_min > imers_or->hdist_filt) {
      continue;
    }
    mi->optimize_likelihood(llhfunc);
    if (mi->d_llh <= mi_closest->d_llh) {
      nd_closest = nd;
      mi_closest = mi;
    }
    node_to_minfo.emplace(nd, mi);
  }
  for (auto [nd, mi] : imers_rc->leaf_to_minfo) {
    mi->mismatch_count = onmers - mi->match_count;
    // mi->compute_gamma();
    if (mi->hdist_min > imers_rc->hdist_filt) {
      continue;
    }
    mi->optimize_likelihood(llhfunc);
    if (mi->d_llh <= mi_closest->d_llh) {
      nd_closest = nd;
      mi_closest = mi;
    }
    node_to_minfo[nd] = mi;
    // If both in reverse-complement and the original sequence, decide:
    if ((imers_or->leaf_to_minfo).contains(nd)) {
      minfo_sptr_t mi_or = (imers_or->leaf_to_minfo)[nd];
      if ((mi->d_llh > mi_or->d_llh) || ((mi->d_llh == mi_or->d_llh) && (mi->match_count < mi_or->match_count))) {
        node_to_minfo[nd] = mi_or;
      }
    }
  }
  if (nd_closest != tree->get_root()) {
    node_to_minfo[nd_closest] = mi_closest;
  }
}

void IBatch::estimate_distances(strstream& batch_stream)
{
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    imers_sptr_t imers_or = std::make_shared<IMers>(index, len, hdist_th);
    imers_sptr_t imers_rc = std::make_shared<IMers>(index, len, hdist_th);

    search_mers(seq, len, imers_or, imers_rc);
    summarize_matches(imers_or, imers_rc);
    batch_stream.precision(STRSTREAM_PRECISION);
    batch_stream << std::fixed;
    report_distances(batch_stream);
  }
}

void IBatch::report_distances(strstream& batch_stream)
{
  if (summarize) {
    vec<node_sptr_t> nd_v;
    nd_v.reserve(node_to_minfo.size());
    for (auto& [nd, mi] : node_to_minfo) {
      mi->chisq = mi_closest->likelihood_ratio(mi->d_llh, llhfunc);
      if (mi->chisq < chisq_value && mi->d_llh < dist_max) {
        nd_v.push_back(nd);
      }
    }
    for (node_sptr_t& nd : nd_v) {
      node_to_wcount[nd] += 1.0 / nd_v.size();
    }
  } else {
    if (node_to_minfo.empty() || (!no_filter && (mi_closest->d_llh > dist_max))) {
      batch_stream << identifer_batch[bix] << "\tNaN\tNaN\n";
      return;
    }
    if (multi) {
      if (no_filter) {
        for (const auto& [nd, mi] : node_to_minfo) {
          batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd, mi);
        }
      } else {
        for (auto& [nd, mi] : node_to_minfo) {
          mi->chisq = mi_closest->likelihood_ratio(mi->d_llh, llhfunc);
          if (mi->chisq < chisq_value && mi->d_llh < dist_max) {
            batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd, mi);
          }
        }
      }
    } else {
      batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd_closest, mi_closest);
    }
  }
}

void IBatch::place_sequences(strstream& batch_stream, bool tabular)
{
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    imers_sptr_t imers_or = std::make_shared<IMers>(index, len, hdist_th);
    imers_sptr_t imers_rc = std::make_shared<IMers>(index, len, hdist_th);

    search_mers(seq, len, imers_or, imers_rc);
    summarize_matches(imers_or, imers_rc);
    batch_stream.precision(STRSTREAM_PRECISION);
    batch_stream << std::fixed;
    report_placement(batch_stream, tabular);
  }
}

void IBatch::report_placement(strstream& batch_stream, bool tabular)
{
  if (node_to_minfo.size() == 0 || !(no_filter || (mi_closest->get_leq_tau(tau) > 1.0))) {
    return;
  }
  node_sptr_t nd_pp = nd_closest;
  minfo_sptr_t mi_pp = mi_closest;
  mi_pp->chisq = 0;

  if (!tabular) batch_stream << "\t\t\t{\"n\" : [\"" << identifer_batch[bix] << "\"], \"p\" : [";
  if (node_to_minfo.size() == 1) {
    if (summarize) {
      node_to_wcount[nd_pp] += 1.0;
    } else {
      if (tabular) {
        batch_stream << identifer_batch[bix] << "\t" << TABULAR_FIELD(nd_pp, mi_pp) << "\n";
      } else {
        batch_stream << PLACEMENT_FIELD(nd_pp, mi_pp) << "]},\n";
      }
    }
    return;
  }

  vec<node_sptr_t> nd_v;
  nd_v.reserve(node_to_minfo.size());
  parallel_flat_phmap<node_sptr_t, minfo_sptr_t> pp_map;

  for (auto& [nd_curr, mi_curr] : node_to_minfo) {
    pp_map[nd_curr] = mi_curr;
    double denom = 1.0;
    node_sptr_t nd_parent = nd_curr;
    // assert(!(std::isnan(mi_curr->d_llh) || std::isnan(mi_curr->v_llh)));
    while ((nd_parent = nd_parent->get_parent())) {
      if (nd_parent->check_taxon() && nd_curr->check_taxon()) {
        denom = 1.0;
        // denom /= nd_parent->get_nchildren();
      } else {
        denom /= nd_parent->get_nchildren();
      }
      if (!pp_map.contains(nd_parent)) {
        pp_map[nd_parent] = std::make_shared<Minfo>(hdist_th);
      }
      pp_map[nd_parent]->add(mi_curr, denom);
    }
  }

  // Collect candidate placements.
  for (auto& [nd_curr, mi_curr] : pp_map) {
    if (no_filter || (mi_curr->get_leq_tau(tau) > 1.0)) {
      if (!nd_curr->check_leaf()) {
        mi_curr->optimize_likelihood(llhfunc);
      }
      mi_curr->chisq = mi_closest->likelihood_ratio(mi_curr->d_llh, llhfunc);
      if ((mi_curr->chisq < chisq_value) && nd_curr->get_parent()) {
        nd_v.push_back(nd_curr);
      }
    }
  }
  // assert(nd_v.size() > 0);

  double total_lwr = 0;
  for (uint32_t i = 0; i < nd_v.size(); ++i) {
    nd_pp = nd_v[i];
    mi_pp = pp_map[nd_pp];
    mi_pp->lwr = exp(-mi_pp->chisq / 2);
    total_lwr = total_lwr + mi_pp->lwr;
  }

  if (multi) {
    for (uint32_t i = 0; i < nd_v.size(); ++i) {
      nd_pp = nd_v[i];
      mi_pp = pp_map[nd_pp];
      mi_pp->lwr = mi_pp->lwr / total_lwr;
      if (summarize) {
        node_to_wcount[nd_pp] += 1.0 / nd_v.size();
      } else {
        if (i > 0 && !tabular) batch_stream << ",";
        if (tabular) {
          batch_stream << identifer_batch[bix] << "\t" << TABULAR_FIELD(nd_pp, mi_pp) << "\n";
        } else {
          batch_stream << "\n\t\t\t\t" << PLACEMENT_FIELD(nd_pp, mi_pp);
        }
      }
    }
    if (!summarize && !tabular) {
      batch_stream << "]\n\t\t\t},\n";
    }
  } else {
    if (nd_v.size() > 1) {
      // Sort: prefer higher card, then lower d_llh
      std::sort(nd_v.begin(), nd_v.end(), [&](node_sptr_t lhs, node_sptr_t rhs) {
        return (lhs->get_card() == rhs->get_card()) ? pp_map[lhs]->d_llh > pp_map[rhs]->d_llh
                                                    : lhs->get_card() < rhs->get_card();
      });
    }
    nd_pp = nd_v.back();
    mi_pp = pp_map[nd_pp];
    mi_pp->lwr = mi_pp->lwr / total_lwr;
    if (summarize) {
      node_to_wcount[nd_pp] += 1.0;
    } else {
      if (tabular) {
        batch_stream << identifer_batch[bix] << "\t" << TABULAR_FIELD(nd_pp, mi_pp) << "\n";
      } else {
        batch_stream << PLACEMENT_FIELD(nd_pp, mi_pp) << "]},\n";
      }
    }
  }
}

IMers::IMers(index_sptr_t index, uint64_t len, uint32_t hdist_th)
  : index(index)
  , len(len)
  , hdist_th(hdist_th)
  , onmers(0)
{
  lshf = index->get_lshf();
  tree = index->get_tree();
  k = lshf->get_k();
  h = lshf->get_h();
  if (len) {
    enmers = len - k + 1;
  } else {
    enmers = 0;
  }
}

void IMers::add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  se_t se;
  pse_t pse;
  node_sptr_t nd;
  uint32_t hdist_curr;
  std::queue<se_t> se_q;
  std::pair<vec_cmer_it, vec_cmer_it> indices = index->bucket_indices(rix);
  crecord_sptr_t crecord = index->get_crecord(rix);
  for (; indices.first < indices.second; ++indices.first) {
    hdist_curr = popcount_lr32(indices.first->first ^ enc_lr);
    if (hdist_curr > hdist_th) {
      continue;
    }
    if (hdist_curr < hdist_filt) {
      hdist_filt = hdist_curr;
    }
    se_q.push(indices.first->second);
    while (!se_q.empty()) {
      se = se_q.front();
      se_q.pop();
      if (tree->check_node(se)) {
        if (!(nd = tree->get_node(se))) {
          continue;
        } else if (nd->check_leaf()) {
          if (!leaf_to_minfo.contains(nd)) {
            leaf_to_minfo[nd] = std::make_shared<Minfo>(hdist_th, enmers, crecord->get_rho(se));
          }
          leaf_to_minfo[nd]->update_match(indices.first->first, pos, hdist_curr);
          continue;
        }
      }
      pse = crecord->get_pse(se);
      se_q.push(pse.first);
      se_q.push(pse.second);
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
/*       ugamma = imers->k; */
/*       continue; */
/*     } */
/*     if (match_v[i].pos > match_v[i - 1].pos) { */
/*       s = (match_v[i].pos - match_v[i - 1].pos); */
/*     } else { */
/*       s = (match_v[i - 1].pos - match_v[i].pos); */
/*     } */
/*     if (s > imers->k) { */
/*       ugamma += imers->k; */
/*     } else { */
/*       ugamma += s; */
/*     } */
/*   } */
/*   gamma = ugamma / (nmers + imers->k - 1); */
/* } */

double Minfo::likelihood_ratio(double d, optimize::HDistHistLLH& llhfunc)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  return 2 * (llhfunc(d) - v_llh);
}

void Minfo::optimize_likelihood(optimize::HDistHistLLH& llhfunc)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  // Locating Function Minima using Brent's algorithm, depends on boost::math.
  std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.33, 13);
  d_llh = sol_r.first;
  v_llh = sol_r.second;
}
