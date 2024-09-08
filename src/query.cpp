#include "query.hpp"
#include "hdhistllh.hpp"

QMers::QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th, double min_covpos)
  : library(library)
  , len(len)
  , num_kmers(0)
  , hdist_th(hdist_th)
  , min_covpos(min_covpos)
{
  lshf = library->get_lshf();
  tree = library->get_tree();
  k = lshf->get_k();
  h = lshf->get_h();
}

void QBatch::search_batch(uint32_t hdist_th, double min_covpos)
{
#pragma omp parallel for num_threads(num_threads)
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, hdist_th, min_covpos);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, hdist_th, min_covpos);

    search_mers(seq, len, qmers_or, qmers_rc);

    qmers_or->summarize_mers();
    qmers_rc->summarize_mers();

#pragma omp critical
    {
      print_summary(qmers_or, qmers_rc, bix);
      /* print_matches(qmers_or, qmers_rc, bix); */
      /* place_wrt_closest(qmers_or, qmers_rc, bix); */
    }
  }
}

void Minfo::summarize_matches()
{
  avg_hdist = 0;
  sup_hdist = 0;
  uint32_t s;
  uint32_t i, j, k;
  uint32_t pos_diff;
  for (i = 0; i < match_v.size(); ++i) {
    for (j = 0; j < qmers->k; ++j) {
      homoc_v[match_v[i].pos + j]++;
    }
    for (j = 0, pos_diff = 0; j < match_v[i].hdist; ++j) {
      k = qmers->lshf->get_npos_accdiff(match_v[i].zc, pos_diff);
      /* obsbp_v[match_v[i].pos + k] += */
      /*   (match_v[i].enc_lr & (0x00010001 << pos_diff)) >> (pos_diff - subsc_v[match_v[i].pos + k]); */
      subsc_v[match_v[i].pos + k]++;
    }
    hdisthist_v[match_v[i].hdist]++;
    avg_hdist += match_v[i].hdist;
    sup_hdist = std::max(sup_hdist, match_v[i].hdist);
    if (i == 0) {
      covpos = qmers->k;
      covmer = match_count;
      continue;
    }
    if (match_v[i].pos > match_v[i - 1].pos) {
      s = (match_v[i].pos - match_v[i - 1].pos);
    } else {
      s = (match_v[i - 1].pos - match_v[i].pos);
    }
    if (s > qmers->k) {
      covpos += qmers->k;
    } else {
      covpos += s;
    }
  }
  covpos /= qmers->len;
  covmer /= (qmers->len - qmers->k + 1);
  avg_hdist /= match_count;
}

void Minfo::estimate_distance(optimize::HDistHistLLH& llhfunc)
{
  /* #pragma omp critical */
  /*   { */
  /*     for (uint32_t i = 0; i < obsbp_v.size(); ++i) { */
  /*       if (subsc_v[i] > 1) { */
  /*         uint32_t bp = 0x00010001 & obsbp_v[i]; */
  /*         for (uint32_t j = 1; j < subsc_v[i]; ++j) { */
  /*           if (bp != (0x00010001 & (obsbp_v[i] >> 1))) { */
  /*             std::cout << std::bitset<32>(obsbp_v[i]) << std::endl; */
  /*             std::cout << "i: " << i << std::endl; */
  /*             std::cout << static_cast<double>(subsc_v[i]) << " / " << static_cast<double>(homoc_v[i]) */
  /*                       << std::endl; */
  /*           } */
  /*           bp = (0x00010001 & (obsbp_v[i] >> 1)); */
  /*         } */
  /*       } */
  /*     } */
  /*   } */
  d_was = 0;
  d_llh = 0;
  llhfunc.set_mc(hdisthist_v.data());
  llhfunc.set_ro(static_cast<double>(ro)); // TODO: fix this
  optimize::Lbfgsb solver;
  optimize::State state =
    solver.minimize(llhfunc, optimize::d_init, optimize::d_lb, optimize::d_ub);
  d_llh = (state.x().transpose())(0);
  uint32_t ddd = 0;
  for (uint32_t i = 0; i < qmers->len; ++i) {
    if (homoc_v[i] > 0) {
      d_was += (static_cast<double>(subsc_v[i]) / static_cast<double>(homoc_v[i]));
      ddd++;
    } else {
      /* d_was += d_llh; */
    }
  }
  d_was = d_was / ddd;
}

void QBatch::place_wrt_closest(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix)
{
  node_sptr_t nd_placement = nullptr;
  double dmin_was = std::numeric_limits<double>::max();
  std::cout << name_batch[bix] << "\t";
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    if (mi->covpos > qmers_or->min_covpos && mi->d_was < dmin_was) {
      nd_placement = nd;
      dmin_was = mi->d_was;
    }
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->covpos > qmers_or->min_covpos && mi->d_was < dmin_was) {
      nd_placement = nd;
      dmin_was = mi->d_was;
    }
  }
  if (nd_placement) {
    nd_placement = nd_placement->get_parent();
    std::cout << nd_placement->get_name() << "\n";
  } else {
    std::cout << "UP"
              << "\n";
  }
}

void QMers::summarize_mers()
{
  optimize::HDistHistLLH llhfunc(h, k, hdist_th, num_kmers);
  for (auto& [nd, mi] : node_to_minfo) {
    mi->summarize_matches();
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
  uint32_t curr_hdist, curr_zc;
  std::queue<se_t> qsubset;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = library->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->get_next(rix);
  crecord_sptr_t crecord = library->get_crecord(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_zc = zc_lr32(iter1->first, enc_lr);
    curr_hdist = __builtin_popcount(curr_zc);
    if (curr_hdist > hdist_th) {
      continue;
    }
    qsubset.push(iter1->second);
    while (!qsubset.empty()) {
      se = qsubset.front();
      qsubset.pop();
      if (!crecord->check_node(se)) {
        pse = crecord->get_pse(se);
        qsubset.push(pse.first);
        qsubset.push(pse.second);
        continue;
      }
      nd = crecord->get_node(se);
      if (nd->check_leaf() && (nd->get_name() != leave_out_ref)) { // TODO: Remove testing.
        if (!node_to_minfo.contains(nd)) {
          node_to_minfo[nd] = std::make_unique<Minfo>(getptr(), crecord->get_wdensity(se));
        }
        node_to_minfo[nd]->update_match(iter1->first, pos, curr_zc, curr_hdist);
      } else { // TODO: this might not be needed.
        for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
          qsubset.push((*std::next(nd->get_children(), i))->get_se());
        }
      }
    }
  }
  num_kmers++;
}

void QBatch::print_summary(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix)
{
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    /* if (mi->covpos > qmers_or->min_covpos && !(qmers_rc->node_to_minfo.contains(nd) && */
    /*                                            mi->covpos < qmers_rc->node_to_minfo[nd]->covpos)) { */
    std::cout << name_batch[bix] << "\t"
              << "or"
              << "\t" << mi->match_count << "\t" << mi->qmers->len << "\t" << mi->ro << "\t"
              << nd->get_name() << "\t" << mi->d_was << "\t" << mi->d_llh << "\t" << mi->covpos
              << "\t" << mi->covmer;
    for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
      std::cout << "\t" << mi->hdisthist_v[i];
    }
    std::cout << "\n";
    /* } */
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    /* if (mi->covpos > qmers_rc->min_covpos && !(qmers_or->node_to_minfo.contains(nd) && */
    /*                                            mi->covpos < qmers_or->node_to_minfo[nd]->covpos)) { */
    std::cout << name_batch[bix] << "\t"
              << "rc"
              << "\t" << mi->match_count << "\t" << mi->qmers->len << "\t" << mi->ro << "\t"
              << nd->get_name() << "\t" << mi->d_was << "\t" << mi->d_llh << "\t" << mi->covpos
              << "\t" << mi->covmer;
    for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
      std::cout << "\t" << mi->hdisthist_v[i];
    }
    std::cout << "\n";
    /* } */
  }
}

void QBatch::print_matches(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix)
{
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    // TODO: Remove leave_out_ref.
    if (mi->covpos > qmers_or->min_covpos && !(qmers_rc->node_to_minfo.contains(nd) &&
                                               mi->covpos < qmers_rc->node_to_minfo[nd]->covpos)) {
      std::cout << leave_out_ref << "\t" << gp_hash(name_batch[bix]) * bix << "\t"
                << nd->get_name();
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        std::cout << "\t" << mi->hdisthist_v[i];
      }
      std::cout << "\n";
    }
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->covpos > qmers_rc->min_covpos && !(qmers_or->node_to_minfo.contains(nd) &&
                                               mi->covpos < qmers_or->node_to_minfo[nd]->covpos)) {
      std::cout << leave_out_ref << "\t" << gp_hash(name_batch[bix]) * bix << "\t"
                << nd->get_name();
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        std::cout << "\t" << mi->hdisthist_v[i];
      }
      std::cout << "\n";
    }
  }
}
