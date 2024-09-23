#include "query.hpp"
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

void QBatch::place_batch(uint32_t hdist_th, double min_gamma)
{
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, hdist_th, min_gamma);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, hdist_th, min_gamma);

    search_mers(seq, len, qmers_or, qmers_rc);

    qmers_or->summarize_minfo();
    qmers_rc->summarize_minfo();

    node_sptr_t nd_placement = place_wrt_tau(qmers_or, qmers_rc);
    /* nd_placement = place_wrt_closest(qmers_or, qmers_rc); */
  }
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
  std::queue<se_t> se_q;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = library->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->get_next(rix);
  crecord_sptr_t crecord = library->get_crecord(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist > hdist_th) {
      continue;
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
      if (nd->check_leaf() && (nd->get_name() != leave_out_ref)) { // TODO: Remove testing.
        if (!node_to_minfo.contains(nd)) {
          node_to_minfo[nd] = std::make_unique<Minfo>(getptr(), crecord->get_rho(se));
        }
        node_to_minfo[nd]->update_match(iter1->first, pos, curr_hdist);
      } else { // TODO: This might not be needed.
        for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
          se_q.push((*std::next(nd->get_children(), i))->get_se());
        }
      }
    }
  }
  nmers++;
}

node_sptr_t QBatch::place_wrt_closest(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  node_sptr_t nd_placement = nullptr;
  double dmin_llh = std::numeric_limits<double>::max();
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    if (mi->gamma > qmers_or->min_gamma && mi->d_llh < dmin_llh) {
      nd_placement = nd;
      dmin_llh = mi->d_llh;
    }
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->gamma > qmers_or->min_gamma && mi->d_llh < dmin_llh) {
      nd_placement = nd;
      dmin_llh = mi->d_llh;
    }
  }
  return nd_placement;
}

node_sptr_t QBatch::place_wrt_tau(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  btree_phmap<se_t, double> se_to_dmin;
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    if (mi->gamma < qmers_or->min_gamma) {
      continue;
    }
    se_to_dmin[nd->get_se()] = mi->d_llh;
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->gamma < qmers_or->min_gamma) {
      continue;
    }
    if (!(se_to_dmin.contains(nd->get_se()) && se_to_dmin[nd->get_se()] < mi->d_llh)) {
      se_to_dmin[nd->get_se()] = mi->d_llh;
    }
  }

  if (se_to_dmin.size() < 2) {
    return tree->get_node((se_to_dmin.begin())->first);
  }

  node_sptr_t subtree_root = Tree::compute_lca(tree->get_node((se_to_dmin.begin())->first),
                                               tree->get_node((se_to_dmin.end())->first));
  se_t se_max = subtree_root->get_se();
  se_t se_min = (se_to_dmin.begin())->first;

  double tau;
  double blen_curr, blen_ddiff, blen_sdiff;
  vec<bool> subtree_v(se_max + 1, false);
  node_sptr_t nd_curr, nd_prev;

  vec<double> dmin_v;
  vec<double> blen_v;
  vec<se_t> se_v;
  dmin_v.reserve(se_to_dmin.size());
  blen_v.reserve(se_to_dmin.size());
  se_v.reserve(se_to_dmin.size());

  for (tuint_t se = se_min; se < se_max; ++se) {
    nd_curr = tree->get_node(se);
    if (se_to_dmin.contains(se)) {
      subtree_v[se] = true;
      subtree_v[nd_curr->get_parent()->get_se()] = true;
    }
    if (!subtree_v[se]) {
      continue;
    }
    if (se == se_min) {
      for (const auto& kv : se_to_dmin) {
        dmin_v.push_back(kv.second);
        blen_curr = Tree::compute_distance(tree->get_node(kv.first), nd_curr->get_parent());
        blen_curr = kv.first == se_min ? blen_curr - nd_curr->get_blen() / 2.0
                                       : blen_curr + nd_curr->get_blen() / 2.0;
        blen_v.push_back(blen_curr);
        se_v.push_back(se);
      }
    } else {
      if (nd_curr == nd_prev->get_parent()) {
        blen_ddiff = nd_curr->get_blen() / 2.0 + nd_prev->get_blen() / 2.0;
        blen_sdiff = -blen_ddiff;
      } else {
        blen_ddiff = nd_prev->get_blen() / 2.0 - nd_curr->get_blen() / 2.0;
        blen_sdiff = -nd_prev->get_blen() / 2.0 - nd_curr->get_blen() / 2.0;
      }
      for (uint32_t i = 0; i < se_v.size(); ++i) {
        if (se_v[i] == nd_prev->get_se()) {
          blen_v[i] -= blen_sdiff;
        } else {
          blen_v[i] -= blen_ddiff;
        }
      }
      nd_prev = nd_curr;
    }
    tau = kendalls_tau(dmin_v.data(), blen_v.data(), se_v.size());
  }

  node_sptr_t nd_placement = nullptr;

  return nd_placement;
}

void optimize::simulate_hdhistllh()
{
  const uint32_t num_replicates = 1000;
  const uint32_t hdist_th = 4;
  const double min_gamma = 0.0;
  const double dist_max = 0.25;
  const uint32_t len = 150;
  const double rho = 1.0;
  const uint32_t k = 29;
  const uint32_t h = 13;
  uint32_t hdist_max = (len * dist_max);

#pragma omp parallel for num_threads(num_threads)
  for (uint32_t hdist_curr = 0; hdist_curr < hdist_max; ++hdist_curr) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dindexed(rho);

    uint32_t hdist_seq;
    double d_llh;

    vec<bool> subs_v(len, false);
    vec<uint32_t> hdisthist_v(hdist_th + 1, 0);
    vec<uint32_t> pos_v(len);
    std::iota(std::begin(pos_v), std::end(pos_v), 0);

    optimize::HDistHistLLH llhfunc(h, k, hdist_th, (len - k + 1));
    for (uint64_t bix = 0; bix < num_replicates; ++bix) {
      std::random_shuffle(pos_v.begin(), pos_v.end());
      for (uint32_t pix = 0; pix < hdist_curr; ++pix) {
        subs_v[pos_v[pix]] = true;
      }

      for (uint32_t six = 0; six < (len - k + 1); ++six) {
        hdist_seq = 0;
        for (uint32_t mix = 0; mix < k; ++mix) {
          if (subs_v[six + mix]) {
            hdist_seq++;
          }
        }
        std::bernoulli_distribution dcoll(llhfunc.prob_collide(hdist_seq));
        if ((hdist_seq <= hdist_th) && dcoll(gen) && dindexed(gen)) {
          hdisthist_v[hdist_seq]++;
        }
      }

      llhfunc.set_mc(hdisthist_v.data());
      llhfunc.set_ro(rho);
      optimize::Lbfgsb solver;
      optimize::State state =
        solver.minimize(llhfunc, optimize::d_init, optimize::d_lb, optimize::d_ub);
      d_llh = (state.x().transpose())(0);

      std::fill(subs_v.begin(), subs_v.end(), false);
      std::fill(hdisthist_v.begin(), hdisthist_v.end(), 0);
    }
  }
}
