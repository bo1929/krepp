#include "query.hpp"
#include <boost/math/tools/minima.hpp>

/* #define CHISQ_THRESHOLD 3.841 */
/* #define CHISQ_THRESHOLD 2.706 */
/* #define CHISQ_THRESHOLD 1.642 */

namespace {
  constexpr size_t kPendingBlock = 64;

  inline void prefetch(const void* ptr)
  {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_prefetch(ptr);
#else
    (void)ptr;
#endif
  }
} // namespace

IBatch::IBatch(index_sptr_t index,
               qseq_sptr_t qs,
               uint32_t hdist_th,
               double chisq_value,
               double dist_max,
               uint32_t tau,
               bool no_filter,
               bool multi,
               bool summarize,
               bool p_value)
  : index(index)
  , hdist_th(hdist_th)
  , chisq_value(chisq_value)
  , dist_max(dist_max)
  , tau(tau)
  , no_filter(no_filter)
  , multi(multi)
  , summarize(summarize)
  , p_value(p_value)
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
    if (imers_or->stage_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr))) {
      wnmers_or++;
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (imers_rc->stage_mer(len - i, rcrix, lshf->drop_ppos_lr(conv_bp64_lr64(rcenc64_bp)))) {
      wnmers_rc++;
    }
#endif /* CANONICAL */
  }
  imers_or->flush_mers();
  imers_rc->flush_mers();
}

void IBatch::widen_hdist_filter(uint32_t& hdist_filt)
{
  const uint32_t unset = std::numeric_limits<uint32_t>::max();
  if (hdist_filt != unset) {
    hdist_filt = 2 * hdist_filt + 1;
  }
}

void IBatch::summarize_matches(imers_sptr_t imers_or, imers_sptr_t imers_rc)
{
  nd_closest = tree->get_root();
  mi_closest = std::make_shared<Minfo>(hdist_th);
  node_to_minfo.clear();
  widen_hdist_filter(imers_or->hdist_filt);
  widen_hdist_filter(imers_rc->hdist_filt);
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
    batch_stream.precision(5);
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
      if (mi->chisq < chisq_value && (std::isnan(dist_max) || mi->d_llh < dist_max)) {
        nd_v.push_back(nd);
      }
    }
    for (node_sptr_t& nd : nd_v) {
      node_to_wcount[nd] += 1.0 / nd_v.size();
    }
  } else {
    if (node_to_minfo.empty() || (!std::isnan(dist_max) && (mi_closest->d_llh > dist_max))) {
      batch_stream << identifer_batch[bix] << "\tNA\tNaN";
      if (p_value) {
        batch_stream << "\tNaN";
      }
      batch_stream << "\n";
      return;
    }
    if (multi) {
      vec<const std::pair<const node_sptr_t, minfo_sptr_t>*> rows;
      rows.reserve(node_to_minfo.size());
      for (const auto& entry : node_to_minfo) {
        if (p_value || !no_filter) {
          entry.second->chisq = mi_closest->likelihood_ratio(entry.second->d_llh, llhfunc);
        }
        if (no_filter || entry.second->chisq < chisq_value) {
          if (std::isnan(dist_max) || entry.second->d_llh < dist_max) {
            rows.push_back(&entry);
          }
        }
      }
      std::sort(rows.begin(), rows.end(), [](const auto* lhs, const auto* rhs) {
        if (lhs->second->d_llh != rhs->second->d_llh) return lhs->second->d_llh < rhs->second->d_llh;
        return lhs->first->get_se() < rhs->first->get_se();
      });
      for (const auto* entry : rows) {
        batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELDS(entry->first, entry->second);
        if (p_value) {
          append_p_value(batch_stream, entry->second->chisq);
        }
        batch_stream << "\n";
      }
    } else {
      if (p_value) {
        mi_closest->chisq = mi_closest->likelihood_ratio(mi_closest->d_llh, llhfunc);
      }
      batch_stream << identifer_batch[bix] << "\t" << DISTANCE_FIELDS(nd_closest, mi_closest);
      if (p_value) {
        append_p_value(batch_stream, mi_closest->chisq);
      }
      batch_stream << "\n";
    }
  }
}

void IBatch::place_sequences(strstream& batch_stream, bool tabular)
{
  bool has_previous = false;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    imers_sptr_t imers_or = std::make_shared<IMers>(index, len, hdist_th);
    imers_sptr_t imers_rc = std::make_shared<IMers>(index, len, hdist_th);

    search_mers(seq, len, imers_or, imers_rc);
    summarize_matches(imers_or, imers_rc);
    batch_stream.precision(5);
    batch_stream << std::fixed;
    if (report_placement(batch_stream, tabular, has_previous) && !summarize && !tabular) {
      has_previous = true;
    }
  }
}

placement_t IBatch::make_placement(const node_sptr_t& nd, const minfo_sptr_t& mi, const minfo_sptr_t& mi_parent)
{
  placement_t pp;
  pp.node = nd;
  pp.edge_num = nd->get_en();
  double t = std::isnan(nd->get_blen()) ? 0.0 : nd->get_blen();
  double d_c = mi->jukes_cantor_dist();
  double d_p = mi_parent ? mi_parent->jukes_cantor_dist() : std::numeric_limits<double>::quiet_NaN();
  double x, y;
  if (std::isfinite(d_p) && (d_p > 0.0)) {
    double eps = std::numeric_limits<double>::epsilon() * t;
    x = (t > 0.0) ? std::clamp(t * d_c / d_p, eps, t - eps) : 0.0;
    y = std::max(0.0, std::min(d_c - x, d_p - t + x));
  } else {
    x = std::min(t / 2.0, d_c);
    y = std::max(0.0, d_c - x);
  }
  pp.distal_length = x;
  pp.pendant_length = y;
  pp.likelihood = -mi->v_llh;
  pp.like_weight_ratio = mi->lwr;
  pp.distance = mi->d_llh;
  pp.distal_node = nd->get_name(true);
  return pp;
}

bool IBatch::collect_placements(vec<placement_t>& placements)
{
  placements.clear();
  if (node_to_minfo.size() == 0 || !(no_filter || (mi_closest->get_leq_tau(tau) > 1.0))) {
    return false;
  }
  node_sptr_t nd_pp = nd_closest;
  minfo_sptr_t mi_pp = mi_closest;
  mi_pp->chisq = 0;

  if (node_to_minfo.size() == 1) {
    placements.push_back(make_placement(nd_pp, mi_pp, nullptr));
    return true;
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
        denom /= nd_parent->get_eff_nchildren();
      }
      if (!pp_map.contains(nd_parent)) {
        pp_map[nd_parent] = std::make_shared<Minfo>(hdist_th);
      }
      pp_map[nd_parent]->add(mi_curr, denom);
    }
  }

  auto parent_minfo = [&pp_map](const node_sptr_t& nd) -> minfo_sptr_t {
    node_sptr_t nd_parent = nd->get_parent();
    if (!nd_parent || !pp_map.contains(nd_parent)) {
      return nullptr;
    }
    return pp_map[nd_parent];
  };

  // Collect candidate placements.
  for (auto& [nd_curr, mi_curr] : pp_map) {
    if (nd_curr->get_nchildren() != nd_curr->get_eff_nchildren() || nd_curr->get_nchildren() == 1) {
      continue;
    }
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
    placements.reserve(nd_v.size());
    for (uint32_t i = 0; i < nd_v.size(); ++i) {
      nd_pp = nd_v[i];
      mi_pp = pp_map[nd_pp];
      mi_pp->lwr = mi_pp->lwr / total_lwr;
      placements.push_back(make_placement(nd_pp, mi_pp, parent_minfo(nd_pp)));
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
    placements.push_back(make_placement(nd_pp, mi_pp, parent_minfo(nd_pp)));
  }
  return true;
}

bool IBatch::report_placement(strstream& batch_stream, bool tabular, bool has_previous)
{
  vec<placement_t> placements;
  if (!collect_placements(placements)) {
    return false;
  }

  if (summarize) {
    for (const placement_t& pp : placements) {
      node_to_wcount[pp.node] += 1.0 / placements.size();
    }
    return true;
  }

  if (tabular) {
    for (const placement_t& pp : placements) {
      batch_stream << identifer_batch[bix] << "\t" << PP_TABULAR_FIELDS(pp) << "\n";
    }
    return true;
  }

  if (has_previous) batch_stream << ",\n";
  batch_stream << "\t\t\t{\"n\" : [\"" << identifer_batch[bix] << "\"], \"p\" : [";
  // Layout follows the shape of the search, not the number of placements: a
  // single matched node, or single-best mode, prints inline. Multi mode with
  // more than one match breaks across lines even if only one placement
  // survives. collect_placements does not insert into or erase from
  // node_to_minfo, so the branch it took can be re-derived here; both inline
  // paths push exactly one placement, so front() is safe.
  if (node_to_minfo.size() == 1 || !multi) {
    batch_stream << PP_JPLACE_FIELDS(placements.front()) << "]}";
  } else {
    for (uint32_t i = 0; i < placements.size(); ++i) {
      if (i > 0) batch_stream << ",";
      batch_stream << "\n\t\t\t\t" << PP_JPLACE_FIELDS(placements[i]);
    }
    batch_stream << "]\n\t\t\t}";
  }
  return true;
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
  pending_v.reserve(kPendingBlock);
  if (len) {
    enmers = len - k + 1;
  } else {
    enmers = 0;
  }
}

void IMers::add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  add_matching_mer_view(
    pos, rix, enc_lr, index->get_flatht_view(rix), index->get_crecord_view(rix), index->get_numerator_view(rix));
}

bool IMers::stage_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  if (!index->check_partial_view(rix)) {
    return false;
  }
  const FlatHT* flatht = index->get_flatht_view(rix);
  const uint32_t m = lshf->get_m();
  uint32_t row = rix / m;
  const uint32_t numerator = index->get_numerator_view(rix);
  if (numerator > 1) {
    row = row * numerator + (rix % m);
  }
  pending_v.push_back({pos, enc_lr, flatht->bucket_data(row), flatht->bucket_data(row + 1), index->get_crecord_view(rix)});
  const PendingMer& mer = pending_v.back();
  prefetch(mer.first);
  prefetch(mer.first + 8);
  prefetch(mer.crecord);
  if (pending_v.size() >= kPendingBlock) {
    flush_mers();
  }
  return true;
}

void IMers::flush_mers()
{
  for (const PendingMer& mer : pending_v) {
    process_mer(mer);
  }
  pending_v.clear();
}

void IMers::add_matching_mer_view(uint32_t pos,
                                  uint32_t rix,
                                  enc_t enc_lr,
                                  const FlatHT* flatht,
                                  CRecord* crecord,
                                  uint32_t numerator)
{
  const uint32_t m = lshf->get_m();
  uint32_t row = rix / m;
  if (numerator > 1) {
    row = row * numerator + (rix % m);
  }
  const PendingMer mer{pos, enc_lr, flatht->bucket_data(row), flatht->bucket_data(row + 1), crecord};
  process_mer(mer);
}

void IMers::process_mer(const PendingMer& mer)
{
  se_t se;
  pse_t pse;
  node_sptr_t nd;
  uint32_t hdist_curr;
  const enc_t enc_lr = mer.enc_lr;
  const uint32_t pos = mer.pos;
  CRecord* crecord = mer.crecord;
  const se_t nsubsets = crecord->get_nsubsets();
  if (vnd_v.size() < nsubsets) {
    vnd_v.resize(nsubsets, 0);
  }
  for (const cmer_t* first = mer.first; first < mer.last; ++first) {
    hdist_curr = popcount_lr32(first->first ^ enc_lr);
    if (hdist_curr > hdist_th) {
      continue;
    }
    if (hdist_curr < hdist_filt) {
      hdist_filt = hdist_curr;
    }
    if (++tix == 0) {
      std::fill(vnd_v.begin(), vnd_v.end(), 0);
      ++tix;
    }
    se_v.push_back(first->second);
    while (!se_v.empty()) {
      se = se_v.back();
      se_v.pop_back();
      if (se >= nsubsets) {
        error_exit("Invalid ID " + std::to_string(se) + " in the index record.");
      }
      if (vnd_v[se] == tix) {
        continue;
      }
      vnd_v[se] = tix;
      if (tree->check_node(se)) {
        if (!(nd = tree->get_node(se))) {
          continue;
        } else if (nd->check_leaf()) {
          if (!leaf_to_minfo.contains(nd)) {
            leaf_to_minfo[nd] = std::make_shared<Minfo>(hdist_th, enmers, crecord->get_rho(se));
          }
          leaf_to_minfo[nd]->update_match(first->first, pos, hdist_curr);
          continue;
        }
      }
      pse = crecord->get_pse(se);
      if (pse.first >= nsubsets || pse.second >= nsubsets) {
        error_exit("Invalid parent ID in the index record.");
      }
      se_v.push_back(pse.first);
      se_v.push_back(pse.second);
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
  std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16);
  d_llh = sol_r.first;
  v_llh = sol_r.second;
}
