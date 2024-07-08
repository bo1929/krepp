#include "query.hpp"

void QBatch::search_batch(uint32_t max_hdist)
{
  for (uint64_t bix = 0; bix < batch_size; ++bix) {
    char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    qmers_sptr_t qmers_or = std::make_shared<QMers>(library, max_hdist);
    qmers_sptr_t qmers_rc = std::make_shared<QMers>(library, max_hdist);
    search_mers(seq, len, qmers_or, qmers_rc);
  }
}

void QBatch::search_mers(char* seq, uint64_t len, qmers_sptr_t qmers_or, qmers_sptr_t qmers_rc)
{
  uint32_t i, l;
  uint32_t orrix, rcrix;
  uint64_t enc64_bp, enc64_lr, orenc64_bp, orenc64_lr, rcenc64_bp, rcenc64_lr;
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
      compute_encoding(seq + i - k, seq, enc64_lr, enc64_bp);
    } else {
      update_encoding(seq + i - 1, enc64_lr, enc64_bp);
    }
    orenc64_bp = enc64_bp & mask_bp;
    orenc64_lr = enc64_lr & mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    rcenc64_lr = conv_bp64_lr64(rcenc64_bp);
    orrix = lshashf->compute_hash(orenc64_bp);
    rcrix = lshashf->compute_hash(rcenc64_bp);
    if (library->partial_exists(orrix)) {
      qmers_or->add_mer(i - k, orrix, lshashf->drop_ppos_lr(orenc64_lr));
    }
    if (library->partial_exists(rcrix)) {
      qmers_rc->add_mer(i - k, rcrix, lshashf->drop_ppos_lr(rcenc64_lr));
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
  mask_bp = u64m >> (32 - k) * 2;
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->name_batch, name_batch);
  batch_size = qs->batch_size;
}

void QMers::add_mer(uint32_t pos, uint32_t rix, enc_t enc_lr)
{
  uint32_t curr_hdist;
  std::vector<cmer_t>::const_iterator iter1 = library->begin(rix);
  std::vector<cmer_t>::const_iterator iter2 = library->end(rix);
  for (; iter1 < iter2; ++iter1) {
    curr_hdist = hdist_lr32(iter1->first, enc_lr);
    if (curr_hdist <= max_hdist) {
      match_v.emplace_back(iter1->first, iter1->second, curr_hdist);
    }
  }
}
