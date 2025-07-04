#include "query.hpp"
#include <boost/math/tools/minima.hpp>

/* #define CHISQ_THRESHOLD 3.841 */
#define CHISQ_THRESHOLD 2.706
/* #define CHISQ_THRESHOLD 1.642 */

IBatch::IBatch(index_sptr_t index,
               qseq_sptr_t qs,
               uint32_t hdist_th,
               double dist_max,
               uint32_t tau,
               bool no_filter,
               bool multi)
  : index(index)
  , hdist_th(hdist_th)
  , dist_max(dist_max)
  , tau(tau)
  , no_filter(no_filter)
  , multi(multi)
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
    if (mi->d_llh < mi_closest->d_llh) {
      nd_closest = nd;
      mi_closest = mi;
    }
    node_to_minfo[nd] = mi;
  }
  for (auto [nd, mi] : imers_rc->leaf_to_minfo) {
    mi->mismatch_count = onmers - mi->match_count;
    // mi->compute_gamma();
    if (mi->hdist_min > imers_rc->hdist_filt) {
      continue;
    }
    mi->optimize_likelihood(llhfunc);
    if (mi->d_llh < mi_closest->d_llh) {
      nd_closest = nd;
      mi_closest = mi;
    }
    node_to_minfo[nd] = mi;
    // If both in reverse-complement and the original sequence, decide:
    if ((imers_or->leaf_to_minfo).contains(nd) &&
        mi->d_llh > (imers_or->leaf_to_minfo)[nd]->d_llh) {
      node_to_minfo[nd] = (imers_or->leaf_to_minfo)[nd];
    }
  }
}

void IBatch::estimate_distances(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    imers_sptr_t imers_or = std::make_shared<IMers>(index, len, hdist_th);
    imers_sptr_t imers_rc = std::make_shared<IMers>(index, len, hdist_th);

    search_mers(seq, len, imers_or, imers_rc);
    summarize_matches(imers_or, imers_rc);
    report_distances(batch_stream);
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void IBatch::report_distances(strstream& batch_stream)
{
  if (node_to_minfo.empty() || (!no_filter && (mi_closest->d_llh > dist_max))) {
    batch_stream << identifer_batch[bix] << "\tNaN\tNaN\n";
  } else if (!multi) {
    batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd_closest, mi_closest) << "\n";
  } else if (no_filter) {
    for (auto& [nd, mi] : node_to_minfo) {
      batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd, mi) << "\n";
    }
  } else {
    for (auto& [nd, mi] : node_to_minfo) {
      mi->chisq = mi_closest->likelihood_ratio(mi->d_llh, llhfunc);
      if (mi->chisq < CHISQ_THRESHOLD && mi->d_llh < dist_max) {
        batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELD(nd, mi) << "\n";
      }
    }
  }
}

void IBatch::place_sequences(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    imers_sptr_t imers_or = std::make_shared<IMers>(index, len, hdist_th);
    imers_sptr_t imers_rc = std::make_shared<IMers>(index, len, hdist_th);

    search_mers(seq, len, imers_or, imers_rc);
    summarize_matches(imers_or, imers_rc);
    report_placement(batch_stream);
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void IBatch::report_placement(strstream& batch_stream)
{
  if (node_to_minfo.size() == 0 || !(no_filter || mi_closest->get_leq_tau(tau) > 1.0)) {
    return;
  }

  node_sptr_t nd_pp = nd_closest;
  minfo_sptr_t mi_pp = mi_closest;
  mi_pp->chisq = 0;
  vec<node_sptr_t> nd_v = {};

  if (node_to_minfo.size() > 1) {
    node_sptr_t nd_curr = nullptr;
    minfo_sptr_t mi_curr = nullptr;

    // Traverse tree in post-order, collect candidate placements
    while ((nd_curr = tree->next_post_order(nd_curr))) {
      if (nd_curr->check_leaf()) {
        if (!node_to_minfo.contains(nd_curr)) {
          node_to_minfo[nd_curr] = std::make_shared<Minfo>(hdist_th, enmers, 0.0);
        } else {
          mi_curr = node_to_minfo[nd_curr];
          mi_curr->chisq = mi_closest->likelihood_ratio(mi_curr->d_llh, llhfunc);
          if (mi_curr->chisq < CHISQ_THRESHOLD && (no_filter || mi_curr->get_leq_tau(tau) > 1.0)) {
            nd_v.push_back(nd_curr);
          }
        }
      } else {
        mi_curr = std::make_shared<Minfo>(hdist_th);
        for (uint32_t nix = 0; nix < nd_curr->get_nchildren(); ++nix) {
          mi_curr->join(node_to_minfo[*std::next(nd_curr->get_children(), nix)]);
        }
        if (mi_curr->rmatch_count > 1.0 && (no_filter || mi_curr->get_leq_tau(tau) > 1.0)) {
          mi_curr->optimize_likelihood(llhfunc);
          mi_curr->chisq = mi_closest->likelihood_ratio(mi_curr->d_llh, llhfunc);
          if (mi_curr->chisq < CHISQ_THRESHOLD) {
            nd_v.push_back(nd_curr);
          }
        }
        node_to_minfo[nd_curr] = mi_curr;
      }
    }

    if (!multi) {
      // Sort: prefer higher card, then higher d_llh
      std::sort(nd_v.begin(), nd_v.end(), [&](node_sptr_t lhs, node_sptr_t rhs) {
        return (lhs->get_card() == rhs->get_card())
                 ? node_to_minfo[lhs]->d_llh > node_to_minfo[rhs]->d_llh
                 : lhs->get_card() < rhs->get_card();
      });
      nd_pp = nd_v.back();
      mi_pp = node_to_minfo[nd_pp];

      // Don't place on top if all matches,
      // choose either of the children based on the total branch length.
      // TODO: Do something more elegant for such cases...
      if ((mi_pp->rcard == mi_pp->rmatch_count) && !nd_pp->check_leaf()) {
        nd_curr = *std::next(nd_pp->get_children(), 0);
        mi_curr = node_to_minfo[nd_curr];
        for (uint32_t nix = 1; nix < nd_pp->get_nchildren(); ++nix) {
          if (nd_curr->total_blen < (*std::next(nd_pp->get_children(), nix))->total_blen) {
            nd_curr = *std::next(nd_pp->get_children(), nix);
            mi_curr = node_to_minfo[nd_curr];
          }
        }
        nd_pp = nd_curr;
        mi_pp = mi_curr;
        if (std::isnan(mi_pp->chisq)) {
          mi_pp->optimize_likelihood(llhfunc);
          mi_pp->chisq = mi_closest->likelihood_ratio(mi_pp->d_llh, llhfunc);
        }
      }
    }
  }
  batch_stream << "\t\t\t{\"n\" : [\"" << identifer_batch[bix] << "\"], \"p\" : [\n";
  if (multi) {
    for (uint32_t i = 0; i < nd_v.size(); ++i) {
      nd_pp = nd_v[i];
      mi_pp = node_to_minfo[nd_pp];
      if (i > 0) batch_stream << ",\n";
      batch_stream << PLACEMENT_FIELD(nd_pp, mi_pp);
    }
  } else {
    batch_stream << PLACEMENT_FIELD(nd_pp, mi_pp);
  }
  batch_stream << "]\n\t\t\t},\n";
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
  node_sptr_t nd;
  uint32_t hdist_curr;
  std::queue<se_t> se_q;
  std::pair<se_t, se_t> pse;
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
      if (!tree->check_node(se)) {
        pse = crecord->get_pse(se);
        se_q.push(pse.first);
        se_q.push(pse.second);
        continue;
      }
      nd = tree->get_node(se);
      if (nd->check_leaf()) {
        if (!leaf_to_minfo.contains(nd)) {
          leaf_to_minfo[nd] = std::make_shared<Minfo>(hdist_th, enmers, crecord->get_rho(se));
        }
        leaf_to_minfo[nd]->update_match(indices.first->first, pos, hdist_curr);
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
}
