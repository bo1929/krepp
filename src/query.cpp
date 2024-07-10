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
    if (qmers_or->check_better_coverage(qmers_rc)) {
      qmers = qmers_or;
    } else {
      qmers = qmers_rc;
    }
    qmers->print_info(name_batch[bix]);
  }
}

void QMers::print_info(std::string name)
{
  for (uint32_t hdist = 0; hdist <= max_hdist; ++hdist) {
    std::cout << std::setprecision(4) << name << "\t" << hdist_to_nnodes[hdist] << "\t" << hdist
              << "\t" << hdist_to_coverage[hdist] << "\t" << hdist_to_ratio[hdist] << std::endl;
  }
}

bool QMers::check_better_coverage(qmers_sptr_t qmers)
{
  bool is_better = true;
  for (uint32_t hdist = 0; hdist <= max_hdist; ++hdist) {
    if (hdist_to_coverage[hdist] > qmers->hdist_to_coverage[hdist]) {
      is_better = true;
      break;
    } else if (hdist_to_coverage[hdist] < qmers->hdist_to_coverage[hdist]) {
      is_better = false;
      break;
    } else {
      continue;
    }
  }
  return is_better;
}

void QMers::compute_coverage()
{
  float ir = 1.0 / len;
  float ik = 1.0 / (len - k + 1);
  std::unordered_map<uint32_t, uint32_t> frag_to_hdist, pos_to_hdist;
  std::unordered_map<uint32_t, std::unordered_set<se_t>> hdist_to_setse;
  for (auto const& [se, matches] : se_to_matches) {
    /* tuint_t card = crecord->get_node(se)->get_card(); */
    for (auto& m : matches) {
      hdist_to_setse[m.hdist].insert(se);
      if (pos_to_hdist.find(m.pos) == pos_to_hdist.end()) {
        pos_to_hdist[m.pos] = m.hdist;
      } else if (pos_to_hdist[m.pos] > m.hdist) {
        pos_to_hdist[m.pos] = m.hdist;
      } else {
        continue;
      }
      for (uint32_t i = 0; i < k; ++i) {
        if (frag_to_hdist.find(i + m.pos) == frag_to_hdist.end()) {
          frag_to_hdist[i + m.pos] = m.hdist;
        } else if (frag_to_hdist[i + m.pos] > m.hdist) {
          frag_to_hdist[i + m.pos] = m.hdist;
        } else {
          continue;
        }
      }
    }
  }
  for (auto const& [hdist, setse] : hdist_to_setse) {
    hdist_to_nnodes[hdist] = setse.size();
  }
  for (auto const& [pos, hdist] : pos_to_hdist) {
    for (uint32_t i = hdist; i <= max_hdist; ++i) {
      hdist_to_ratio[i] += ik;
    }
  }
  for (auto const& [pos, hdist] : frag_to_hdist) {
    for (uint32_t i = hdist; i <= max_hdist; ++i) {
      hdist_to_coverage[i] += ir;
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
      qmers_or->add_mer(i - k, orrix, lshashf->drop_ppos_lr(orenc64_lr));
    }
    if (library->check_partial(rcrix)) {
      qmers_rc->add_mer(len - i, rcrix, lshashf->drop_ppos_lr(rcenc64_lr));
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

void QMers::add_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  se_t se;
  uint32_t curr_hdist;
  std::queue<se_t> q;
  std::pair<se_t, se_t> pse;
  std::vector<cmer_t>::const_iterator iter1 = library->begin(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->end(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist > max_hdist) {
      continue;
    }
    q.push(iter1->second);
    while (!q.empty()) {
      se = q.front();
      q.pop();
      if (crecord->check_node(se)) { // && crecord->get_node(se)->check_leaf()) {
        if (se_to_matches[se].empty()) {
          se_to_matches[se].emplace_back(iter1->first, pos, curr_hdist);
        } else if (se_to_matches[se].back().pos != pos) {
          se_to_matches[se].emplace_back(iter1->first, pos, curr_hdist);
        } else if (se_to_matches[se].back().hdist > curr_hdist) {
          se_to_matches[se].back().hdist = curr_hdist;
          se_to_matches[se].back().enc_lr = iter1->first;
        } else {
          continue;
        }
      } else {
        pse = crecord->get_pse(se);
        q.push(pse.first);
        q.push(pse.second);
      }
    }
  }
}
