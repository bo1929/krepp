#include "common.hpp"
#include "query.hpp"
#include <cstdint>
#include <functional>
#include <unordered_map>

QMers::QMers(library_sptr_t library, uint64_t len, uint32_t max_hdist)
  : library(library)
  , len(len)
  , max_hdist(max_hdist)
{
  crecord = library->get_crecord();
  tree = library->get_tree();
  k = library->get_lshashf()->get_k();
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
    std::pair<node_sptr_t, uint32_t> pmer_or = qmers_or->max_parsimonious_mer();
    std::pair<node_sptr_t, uint32_t> pmer_rc = qmers_rc->max_parsimonious_mer();
    if (pmer_or.second > pmer_rc.second) {
      std::cout << pmer_or.first->get_name() << std::endl;
    } else {
      std::cout << pmer_rc.first->get_name() << std::endl;
    }
  }
}

std::pair<node_sptr_t, uint32_t> QMers::max_parsimonious_mer()
{

  node_sptr_t nd_m = tree->get_root();
  uint32_t max_score;
  while (!nd_m->check_leaf()) {
    std::cout << nd_m->get_name() << "/" << node_to_score[nd_m] << std::endl;
    vec<node_sptr_t> children = nd_m->get_children();
    max_score = 0;
    bool sc = false;
    bool fs = false;
    for (auto nd_c : children) {
      if (node_to_score.find(nd_c) != node_to_score.end()) {
        std::cout << node_to_score[nd_c] << ",";
      } else {
        std::cout << 0 << ",";
      }
    }
    std::cout << std::endl;
    for (auto nd_c : children) {
      if (node_to_score.find(nd_c) != node_to_score.end()) {
        if (sc && node_to_score[nd_c] == max_score) {
          nd_m = nd_m->get_parent();
          fs = true;
          break;
        }
        if (node_to_score[nd_c] > max_score) {
          nd_m = nd_c;
          max_score = node_to_score[nd_c];
        }
        sc = true;
      }
    }
    if (fs || max_score == 0)
      break;
  }
  return std::make_pair(nd_m, max_score);
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
  node_sptr_t nd;
  uint32_t curr_hdist;
  std::queue<se_t> q;
  std::pair<se_t, se_t> pse;
  std::set<se_t> leaf_se_set, node_se_set;
  std::map<node_sptr_t, uint32_t> node_to_dval, node_to_cval;

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
          leaf_se_set.insert(se);
        } else {
          node_se_set.insert(se);
        }
        nd = crecord->get_node(se);
        node_to_dval[nd] = node_to_dval.find(nd) == node_to_dval.end()
                             ? curr_hdist
                             : std::min(node_to_dval[nd], curr_hdist);
      } else {
        pse = crecord->get_pse(se);
        q.push(pse.first);
        q.push(pse.second);
      }
    }
  }

  for (const se_t se : node_se_set) {
    nd = crecord->get_node(se);
    curr_hdist = node_to_dval[nd];
    tree->set_subtree(nd);
    while (nd = tree->next_post_order()) {
      node_to_dval[nd] = node_to_dval.find(nd) == node_to_dval.end()
                           ? curr_hdist
                           : std::min(node_to_dval[nd], curr_hdist);
    }
    tree->reset_traversal();
  }

  for (const se_t se : leaf_se_set) {
    nd = crecord->get_node(se);
    curr_hdist = node_to_dval[nd];
    nd = nd->get_parent();
    while (nd) {
      node_to_dval[nd] = node_to_dval.find(nd) == node_to_dval.end()
                           ? curr_hdist
                           : std::min(node_to_dval[nd], curr_hdist);
      node_se_set.insert(crecord->get_se(nd->get_shash()));
      nd = nd->get_parent();
    }
  }

  for (auto it = node_se_set.rbegin(); it != node_se_set.rend(); ++it) {
    nd = crecord->get_node(*it);
    if (nd->get_parent()) {
      node_to_cval[nd] = node_to_cval[nd->get_parent()];
      for (auto const& nd_c : nd->get_parent()->get_children()) {
        if (nd_c != nd && node_to_dval.find(nd_c) != node_to_dval.end()) {
          node_to_cval[nd] = std::min(node_to_cval[nd], node_to_dval[nd_c]);
        }
      }
    } else {
      node_to_cval[nd] = k;
    }
  }

  for (auto const& [nd, cval] : node_to_cval) {
    /* std::cout << nd->get_name() << ", c: " << node_to_cval[nd] << ", d: " << node_to_dval[nd] */
    /*           << std::endl; */
    if (cval > node_to_dval[nd] && nd != tree->get_root()) {
      node_to_score[nd]++;
    }
  }
  /* std::cout << std::endl; */
  node_to_score[tree->get_root()] = 0;
}
