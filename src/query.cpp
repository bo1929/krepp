#include "query.hpp"

QMers::QMers(library_sptr_t library, uint64_t len, uint32_t max_hdist)
  : library(library)
  , len(len)
  , max_hdist(max_hdist)
{
  crecord = library->get_crecord();
  tree = library->get_tree();
  k = library->get_lshashf()->get_k();
}

void QMers::place_baseline() // TODO: Temporary, remove.
{
  std::unordered_map<std::string, float> node_to_scores = {};
  for (auto const& [pos, match_v] : pos_to_matches) {
    std::unordered_map<std::string, uint32_t> node_to_hdist = {};
    for (auto const& match : match_v) {
      node_sptr_t node_curr = crecord->get_node(match.se);
      if ("G001536045" == node_curr->get_name())
        continue;
      node_to_hdist[node_curr->get_name()] = match.hdist;
      while (node_curr = node_curr->get_parent()) {
        if (node_to_hdist.find(node_curr->get_name()) == node_to_hdist.end())
          node_to_hdist[node_curr->get_name()] = match.hdist;
        if (node_to_hdist[node_curr->get_name()] > match.hdist)
          node_to_hdist[node_curr->get_name()] = match.hdist;
      }
    }
    for (auto const& kv : node_to_hdist) {
      node_to_scores[kv.first] += pow((1.0 - static_cast<float>(kv.second) / k), k);
    }
  }
  node_sptr_t node_curr = tree->get_root();
  float threshold_score = node_to_scores[node_curr->get_name()] / 2;
  while (node_to_scores[node_curr->get_name()] > threshold_score && !node_curr->check_leaf()) {
    /* while (!node_curr->check_leaf()) { */
    vec<node_sptr_t> children_nd_v = node_curr->get_children();
    float max_child_score = 0;
    for (auto child : children_nd_v) {
      if (node_to_scores[child->get_name()] > max_child_score) {
        max_child_score = node_to_scores[child->get_name()];
        node_curr = child;
      }
    }
    /* std::cout << node_curr->get_name() << ", " << max_child_score << std::endl; */
  }
  std::cout << node_curr->get_name() << std::endl;
}

void QBatch::search_batch(uint32_t max_hdist)
{
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    qmers_sptr_t qmers;
    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, max_hdist);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, max_hdist);
    search_mers(seq, len, qmers_or, qmers_rc);
    qmers_or->compute_coverage();
    qmers_rc->compute_coverage();
    /* if (qmers_or->check_better_coverage(qmers_rc)) { */
    if (qmers_or->total_coverage > qmers_rc->total_coverage) {
      qmers = qmers_or;
    } else {
      qmers = qmers_rc;
    }
    qmers->place_baseline();
    qmers->print_info(name_batch[bix]);
  }
}

void QMers::print_info(std::string name)
{
  /* for (uint32_t hdist = 0; hdist <= max_hdist; ++hdist) { */
  /*   std::cout << std::setprecision(4) << name << "\t" << hdist_to_nnodes[hdist] << "\t" << hdist */
  /*             << "\t" << hdist_to_coverage[hdist] << "\t" << hdist_to_ratio[hdist] << std::endl; */
  /* } */
  uint32_t ri = rand() * name.size();
  for (auto const& [se, coverage] : se_to_coverage) {
    if (coverage > 0.5)
      std::cout << ri << "\t" << name << "\t" << crecord->get_node(se)->get_name() << "\t"
                << coverage << std::endl;
  }
}

bool QMers::check_better_coverage(qmers_sptr_t qmers)
{
  bool is_better = true;
  /* for (uint32_t hdist = 0; hdist <= max_hdist; ++hdist) { */
  /*   if (hdist_to_coverage[hdist] > qmers->hdist_to_coverage[hdist]) { */
  /*     is_better = true; */
  /*     break; */
  /*   } else if (hdist_to_coverage[hdist] < qmers->hdist_to_coverage[hdist]) { */
  /*     is_better = false; */
  /*     break; */
  /*   } else { */
  /*     continue; */
  /*   } */
  /* } */
  return is_better;
}

void QMers::compute_coverage()
{
  /* uint64_t nk = len - k + 1; */
  float ir = 1.0 / len;
  /* float ik = 1.0 / nk; */
  for (auto const& [se, matches] : se_to_matches) {
    vec<uint32_t> min_hdist_v(len, max_hdist + 1);
    for (auto& m : matches) {
      for (uint32_t i = 0; i < k; ++i) {
        if (min_hdist_v[m.pos + i] > m.hdist) {
          min_hdist_v[m.pos + i] = m.hdist;
        }
      }
    }
    for (uint32_t hdist : min_hdist_v) {
      if (hdist <= max_hdist) {
        se_to_coverage[se] += ir;
      }
      total_coverage += se_to_coverage[se];
    }
  }
  // std::unordered_map<uint32_t, uint32_t> frag_to_hdist, pos_to_hdist;
  // std::unordered_map<uint32_t, std::unordered_set<se_t>> hdist_to_setse;
  // for (auto const& [se, matches] : se_to_matches) {
  //   /* tuint_t card = crecord->get_node(se)->get_card(); */
  //   for (auto& m : matches) {
  //     hdist_to_setse[m.hdist].insert(se);
  //     if (pos_to_hdist.find(m.pos) == pos_to_hdist.end()) {
  //       pos_to_hdist[m.pos] = m.hdist;
  //     } else if (pos_to_hdist[m.pos] > m.hdist) {
  //       pos_to_hdist[m.pos] = m.hdist;
  //     } else {
  //       continue;
  //     }
  //     for (uint32_t i = 0; i < k; ++i) {
  //       if (frag_to_hdist.find(i + m.pos) == frag_to_hdist.end()) {
  //         frag_to_hdist[i + m.pos] = m.hdist;
  //       } else if (frag_to_hdist[i + m.pos] > m.hdist) {
  //         frag_to_hdist[i + m.pos] = m.hdist;
  //       } else {
  //         continue;
  //       }
  //     }
  //   }
  //}
  // for (auto const& [hdist, setse] : hdist_to_setse) {
  //   hdist_to_nnodes[hdist] = 0;
  //   for (auto const& se : setse) {
  //     if (se_to_matches[se].size() > (hdist / 2)) {
  //       hdist_to_nnodes[hdist] += 1;
  //     }
  //   }
  //   /* hdist_to_nnodes[hdist] = setse.size(); */
  // }
  // for (auto const& [pos, hdist] : pos_to_hdist) {
  //   for (uint32_t i = hdist; i <= max_hdist; ++i) {
  //     hdist_to_ratio[i] += ik;
  //   }
  // }
  // for (auto const& [pos, hdist] : frag_to_hdist) {
  //   for (uint32_t i = hdist; i <= max_hdist; ++i) {
  //     hdist_to_coverage[i] += ir;
  //   }
  // }
}

void QBatch::search_mers(char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  uint32_t i, l;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp, rcenc64_lr;
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
    rcenc64_lr = conv_bp64_lr64(rcenc64_bp);
    orrix = lshashf->compute_hash(orenc64_bp);
    rcrix = lshashf->compute_hash(rcenc64_bp);
    if (library->check_partial(orrix)) {
      qmers_or->add_matching_mer(i - k, orrix, lshashf->drop_ppos_lr(orenc64_lr));
    }
    if (library->check_partial(rcrix)) {
      qmers_rc->add_matching_mer(len - i, rcrix, lshashf->drop_ppos_lr(rcenc64_lr));
    }
  }
}

QBatch::QBatch(library_sptr_t library, qseq_sptr_t qs)
  : library(library)
{
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  lshashf = library->get_lshashf();
  k = lshashf->get_k();
  m = lshashf->get_m();
  mask_bp = u64m >> ((32 - k) * 2);
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->name_batch, name_batch);
  batch_size = qs->batch_size;
}

void QMers::add_matching_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  se_t se;
  uint32_t curr_hdist;
  std::queue<se_t> q;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = library->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->get_next(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist > max_hdist) {
      continue;
    }
    q.push(iter1->second);
    while (!q.empty()) {
      se = q.front();
      q.pop();
      if (crecord->check_node(se)) {
        if (crecord->get_node(se)->check_leaf()) {
          if (se_to_matches[se].empty()) {
            se_to_matches[se].emplace_back(iter1->first, se, pos, curr_hdist);
            pos_to_matches[pos].emplace_back(iter1->first, se, pos, curr_hdist);
          } else if (se_to_matches[se].back().pos != pos) {
            se_to_matches[se].emplace_back(iter1->first, se, pos, curr_hdist);
            pos_to_matches[pos].emplace_back(iter1->first, se, pos, curr_hdist);
          } else if (se_to_matches[se].back().hdist > curr_hdist) {
            se_to_matches[se].back().hdist = curr_hdist;
            pos_to_matches[pos].back().hdist = curr_hdist;
            se_to_matches[se].back().enc_lr = iter1->first;
            pos_to_matches[pos].back().enc_lr = iter1->first;
          } else {
            continue;
          }
        } else { // Alternatively, remove this and use the internal directly.
          tuint_t nchildren = crecord->get_node(se)->get_nchildren();
          for (tuint_t i = 1; i <= nchildren; ++i) {
            q.push(se - i);
          }
        }
      } else {
        pse = crecord->get_pse(se);
        q.push(pse.first);
        q.push(pse.second);
      }
    }
  }
}
