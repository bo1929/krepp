#include "query.hpp"

QMers::QMers(library_sptr_t library, uint64_t len, uint32_t hdist_th, float min_covpos)
  : library(library)
  , len(len)
  , hdist_th(hdist_th)
  , min_covpos(min_covpos)
{
  lshf = library->get_lshf();
  tree = library->get_tree();
  k = lshf->get_k();
}

void QBatch::search_batch(uint32_t hdist_th, float min_covpos)
{
#pragma omp parallel for num_threads(num_threads)
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();

    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, len, hdist_th, min_covpos);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, len, hdist_th, min_covpos);

    search_mers(seq, len, qmers_or, qmers_rc);

    qmers_or->summarize_matches();
    qmers_rc->summarize_matches();

#pragma omp critical
    {
      print_summary(qmers_or, qmers_rc, bix);
    }
  }
}

void QMers::summarize_matches()
{
  float irk = 1.0 / (len - k + 1);
  float irp = 1.0 / (len);
  uint32_t ix1, ix2;
  uint32_t itp;
  for (auto const& [nd, mi] : node_to_minfo) {
    mi->hdisthist_v.resize(hdist_th + 1, 0);
    mi->avghdist = 0;
    mi->maxhdist = 0;
    for (uint32_t i = 0; i < mi->match_v.size(); ++i) {
      for (uint32_t j = 0; j < k; ++j) {
        mi->homoc_v[mi->match_v[i].pos + j]++;
      }
      for (uint32_t j = 0; j < mi->match_v[i].hdist; ++j) {
        mi->subsc_v[mi->match_v[i].pos + lshf->get_ppos_diff(mi->match_v[i].zc)]++;
      }
      mi->hdisthist_v[mi->match_v[i].hdist]++;
      mi->avghdist += mi->match_v[i].hdist;
      mi->maxhdist = std::max(mi->maxhdist, mi->match_v[i].hdist);
      if (i == 0) {
        mi->covmer += irk * mi->match_v.size();
        mi->covpos += irp * k;
        continue;
      }
      if ((mi->match_v[i].pos > mi->match_v[i - 1].pos)) {
        ix1 = i;
        ix2 = i - 1;
      } else {
        ix2 = i;
        ix1 = i - 1;
      }
      if (mi->match_v[ix1].pos > (mi->match_v[ix2].pos + k)) {
        mi->covpos += irp * k;
      } else {
        mi->covpos += irp * (mi->match_v[ix1].pos - mi->match_v[ix2].pos);
      }
    }
    mi->avghdist /= mi->match_v.size();
    itp = 0;
    mi->wschdist = 0;
    for (uint32_t i = 0; i < len; ++i) {
      if (mi->homoc_v[i] > 0) {
        mi->wschdist +=
          (static_cast<float>(mi->subsc_v[i]) / static_cast<float>(mi->homoc_v[i])) >= 0.5;
        itp++;
      } else {
        /* mi->wschdist += (mi->maxhdist + 1.0) / k; */
      }
    }
    mi->wschdist = mi->wschdist / itp;
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
    orrix = lshf->compute_hash(orenc64_bp);
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (library->check_partial(orrix)) {
      qmers_or->add_matching_mer(i - k, orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    if (library->check_partial(rcrix)) {
      qmers_rc->add_matching_mer(len - i, rcrix, lshf->drop_ppos_lr(rcenc64_lr));
    }
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
          node_to_minfo[nd] = std::make_unique<minfo_t>(len);
        }
        node_to_minfo[nd]->update_match(pos, curr_zc, curr_hdist);
      } else { // TODO: this might be not needed.
        for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
          qsubset.push((*std::next(nd->get_children(), i))->get_se());
        }
      }
    }
  }
}

void QBatch::print_summary(qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc, uint64_t bix)
{
  for (auto const& [nd, mi] : qmers_or->node_to_minfo) {
    if (mi->covpos > qmers_or->min_covpos && !(qmers_rc->node_to_minfo.contains(nd) &&
                                               mi->covpos < qmers_rc->node_to_minfo[nd]->covpos)) {
      std::cout << name_batch[bix] << "\t"
                << "or"
                << "\t" << nd->get_name() << "\t" << mi->wschdist << "\t" << mi->avghdist << "\t"
                << mi->covpos << "\t" << mi->covmer << "\t";
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        std::cout << "\t" << mi->hdisthist_v[i];
      }
      std::cout << "\n";
    }
  }
  for (auto const& [nd, mi] : qmers_rc->node_to_minfo) {
    if (mi->covpos > qmers_rc->min_covpos && !(qmers_or->node_to_minfo.contains(nd) &&
                                               mi->covpos < qmers_or->node_to_minfo[nd]->covpos)) {
      std::cout << name_batch[bix] << "\t"
                << "or"
                << "\t" << nd->get_name() << "\t" << mi->wschdist << "\t" << mi->avghdist << "\t"
                << mi->covpos << "\t" << mi->covmer << "\t";
      for (uint32_t i = 0; i < mi->hdisthist_v.size(); ++i) {
        std::cout << "\t" << mi->hdisthist_v[i];
      }
      std::cout << "\n";
    }
  }
}

void QMers::print_summary(const std::string& name)
{
  for (auto const& [nd, mi] : node_to_minfo) {
    std::cout << nd->get_name() << "\t" << name << "\t" << mi->wschdist << "\t" << mi->avghdist
              << "\t" << mi->covpos << "\t" << mi->covmer << "\n";
  }
}

void QMers::print_dist(const std::string& name)
{
  for (auto const& [nd, mi] : node_to_minfo) {
    std::cout << nd->get_name() << "\t" << name << "\t" << mi->avghdist << "\t" << mi->wschdist
              << std::endl;
  }
}

void QMers::print_coverage(const std::string& name)
{
  for (auto const& [nd, mi] : node_to_minfo) {
    std::cout << name << "\t" << nd->get_name() << "\t" << mi->covmer << "\t" << mi->covpos
              << std::endl;
  }
}

void QMers::print_matches(const std::string& name)
{
  for (auto const& [nd, mi] : node_to_minfo) {
    for (auto const& match : mi->match_v) {
      std::cout << name << "\t" << nd->get_name() << "\t" << match.pos << "\t" << match.hdist
                << std::endl;
    }
  }
}
