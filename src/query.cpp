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
  match_vvec.resize(len - k + 1);
}

void QBatch::search_batch(uint32_t max_hdist)
{
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, max_hdist);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, max_hdist);

    search_mers(seq, len, qmers_or, qmers_rc);

    qmers_or->compute_coverage();
    /* qmers_or->print_coverage(); */
    qmers_or->compute_pseudoparsimony();
    std::pair<node_sptr_t, uint32_t> pmer_or = qmers_or->max_parsimonious_mer();
    /* std::cout << name_batch[bix] << "\t" << pmer_or.first->get_name() << std::endl; */

    qmers_rc->compute_coverage();
    /* qmers_rc->print_coverage(); */
    qmers_rc->compute_pseudoparsimony();
    std::pair<node_sptr_t, uint32_t> pmer_rc = qmers_rc->max_parsimonious_mer();
    /* std::cout << name_batch[bix] << "\t" << pmer_rc.first->get_name() << std::endl; */

    /* std::cout << pmer_or.second << ", " << pmer_rc.second; */
    if (pmer_or.second > pmer_rc.second) {
      std::cout << name_batch[bix] << "\t" << pmer_or.first->get_name() << std::endl;
    } else {
      std::cout << name_batch[bix] << "\t" << pmer_rc.first->get_name() << std::endl;
    }
  }
}

std::pair<node_sptr_t, uint32_t> QMers::max_parsimonious_mer()
{
  node_sptr_t nd_c, nd_p, nd_m;
  uint32_t max_score;
  nd_m = tree->get_root();
  /* for (auto kv : node_to_sinfo) { */
  /*   std::cout << kv.first->get_name() << ": " << (kv.second).parsimony_score << std::endl; */
  /* } */
  auto it_child = node_to_sinfo.end();
  while (!nd_m->check_leaf()) {
    max_score = 0;
    nd_p = nd_m;
    for (tuint_t i = 0; i < nd_p->get_nchildren(); ++i) {
      nd_c = *std::next(nd_p->get_children(), i);
      it_child = node_to_sinfo.find(nd_c);
      if (it_child == node_to_sinfo.end()) {
        continue;
      }
      if (max_score != 0 && (it_child->second).parsimony_score == max_score) {
        nd_m = nd_m->get_parent();
        break;
      }
      if ((it_child->second).parsimony_score > max_score) {
        nd_m = nd_c;
        max_score = (it_child->second).parsimony_score;
      }
    }
    if (nd_m == nd_p)
      break;
  }
  return std::make_pair(nd_m, max_score);
}

void QMers::compute_pseudoparsimony()
{
  node_sptr_t nd_p;
  for (uint32_t i = 0; i < match_vvec.size(); ++i) {
    std::unordered_map<node_sptr_t, cdist_t> node_to_score;
    auto it_sibling = node_to_score.end();
    auto it_curr = node_to_score.end();
    std::set<se_t> se_set;
    for (const match_t m : match_vvec[i]) {
      if (node_to_sinfo[m.nd].coverage_pos < 0.33)
        continue;
      se_set.insert(m.nd->get_senc());
      node_to_score[m.nd].sub_hdist = m.hdist;
      nd_p = m.nd->get_parent();
      while (nd_p) {
        if (m.hdist < node_to_score[nd_p].sub_hdist) {
          node_to_score[nd_p].sub_hdist = m.hdist;
        }
        se_set.insert(nd_p->get_senc());
        nd_p = nd_p->get_parent();
      }
    }
    for (auto it_se = se_set.rbegin(); it_se != se_set.rend(); ++it_se) {
      it_curr = node_to_score.find(crecord->get_node(*it_se));
      nd_p = it_curr->first->get_parent();
      if (!nd_p || it_curr->first->check_leaf())
        continue;
      (it_curr->second).exc_hdist = node_to_score[nd_p].exc_hdist;
      for (tuint_t i = 0; i < nd_p->get_nchildren(); ++i) {
        it_sibling = node_to_score.find(*std::next(nd_p->get_children(), i));
        if (it_sibling != node_to_score.end() && it_sibling->first != it_curr->first) {
          (it_curr->second).exc_hdist =
            std::min((it_curr->second).exc_hdist, (it_sibling->second).sub_hdist);
        }
      }
    }
    for (it_curr = node_to_score.begin(); it_curr != node_to_score.end(); ++it_curr) {
      /* std::cout << (it_curr->second).sub_hdist << ", " << (it_curr->second).exc_hdist << std::endl; */
      if ((it_curr->second).sub_hdist < (it_curr->second).exc_hdist) {
        node_to_sinfo[it_curr->first].parsimony_score++;
      }
    }
  }
}

void QMers::print_coverage()
{
  for (auto const& [nd, ri] : node_to_sinfo) {
    std::cout << nd->get_name() << ": " << ri.coverage_mer << ", " << ri.coverage_pos << ", "
              << ri.min_hdist << std::endl;
  }
}

void QMers::compute_coverage()
{
  std::unordered_map<node_sptr_t, uint32_t> node_to_pos;
  float irk = 1.0 / (len - k + 1);
  float irp = 1.0 / (len);
  uint32_t pos;
  for (uint32_t i = 0; i < match_vvec.size(); ++i) {
    for (const match_t m : match_vvec[i]) {
      if (node_to_sinfo[m.nd].min_hdist == m.hdist) {
        node_to_sinfo[m.nd].coverage_mer += irk;
        node_to_pos.try_emplace(m.nd, i);
        pos = node_to_pos[m.nd];
        if ((pos == i) || i > (pos + k)) {
          node_to_sinfo[m.nd].coverage_pos += irp * k;
        } else {
          node_to_sinfo[m.nd].coverage_pos += irp * (i - pos);
        }
        node_to_pos[m.nd] = i;
      }
    }
  }
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
  uint32_t ix;
  uint32_t curr_hdist;
  std::queue<se_t> qsubset;
  std::pair<se_t, se_t> pse;
  std::unordered_map<node_sptr_t, uint32_t> node_to_ix;
  std::vector<cmer_t>::const_iterator iter1 = library->get_first(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->get_next(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist > max_hdist) {
      continue;
    }
    qsubset.push(iter1->second);
    while (!qsubset.empty()) {
      se = qsubset.front();
      qsubset.pop();
      if (crecord->check_node(se)) {
        nd = crecord->get_node(se);
        if (nd->check_leaf()) {
          /* if (nd->get_name() == "G001042715") */
          /*   continue; */
          node_to_ix.try_emplace(nd, match_vvec[pos].size());
          ix = node_to_ix[nd];
          if (node_to_ix[nd] == ix) {
            match_vvec[pos].emplace_back(nd, curr_hdist, enc_lr);
          } else if (curr_hdist < match_vvec[pos][ix].hdist) {
            match_vvec[pos][ix].hdist = curr_hdist;
            match_vvec[pos][ix].enc_lr = enc_lr;
          } else {
            continue;
          }
          if (node_to_sinfo[nd].min_hdist > curr_hdist) {
            node_to_sinfo[nd].min_hdist = curr_hdist;
          }
        } else {
          for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
            qsubset.push((*std::next(nd->get_children(), i))->get_senc());
          }
        }
      } else {
        pse = crecord->get_psenc(se);
        qsubset.push(pse.first);
        qsubset.push(pse.second);
      }
    }
  }
}
