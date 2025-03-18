#include "rqseq.hpp"

RSeq::RSeq(std::string input, lshf_sptr_t lshf, uint8_t w, uint32_t r, bool frac)
  : w(w)
  , r(r)
  , frac(frac)
  , lshf(lshf)
{
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  k = lshf->get_k();
  m = lshf->get_m();
  mask_bp = u64m >> ((32 - k) * 2);
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  is_url = std::regex_match(input, url_regexp);
  if (is_url) {
#if defined _WLCURL && _WLCURL == 1
    input_path = download_url(input);
#else
    std::cerr << "Failed to download from URL, compiled without libcurl!" << std::endl;
#endif
  } else {
    input_path = input;
  }
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    std::cerr << "Failed to open thefile at " << input_path << std::endl;
    exit(EXIT_FAILURE);
  }
  kseq = kseq_init(gfile);
}

RSeq::~RSeq()
{
  kseq_destroy(kseq);
  gzclose(gfile);
  if (is_url) {
    std::filesystem::remove(input_path);
  }
}

template<typename T>
void RSeq::extract_mers(vvec<T>& table, sh_t sh)
{
  uint32_t i, l;
  uint32_t rix, rix_res;
  uint8_t ldiff;
  if (w > k) {
    ldiff = w - k + 1;
  } else {
    ldiff = 1;
    w = k;
  }
  uint64_t kix = 0;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  std::vector<std::pair<uint64_t, uint64_t>> lsh_enc_win(ldiff);
  std::pair<uint64_t, uint64_t> cminimizer, pminimizer;
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
    lsh_enc_win[kix % ldiff].first = orenc64_bp & mask_bp;
    lsh_enc_win[kix % ldiff].second = orenc64_lr & mask_lr;
    kix++;
    if ((l < w) && (i != len)) {
      continue;
    }
    cminimizer =
      *std::min_element(lsh_enc_win.begin(),
                        lsh_enc_win.end(),
                        [](std::pair<uint64_t, uint64_t> lhs, std::pair<uint64_t, uint64_t> rhs) {
                          return xur64_hash(lhs.second) < xur64_hash(rhs.second);
                        });
#ifdef CANONICAL
    rcenc64_bp = revcomp_bp64(cminimizer.first, k);
    if (cminimizer.first < rcenc64_bp) {
      cminimizer.first = rcenc64_bp;
      cminimizer.second = conv_bp64_lr64(rcenc64_bp);
    }
#endif /* CANONICAL */
    rix = lshf->compute_hash(cminimizer.first);
    rix_res = rix % m;
    if (frac ? rix_res <= r : rix_res == r) {
      rix = frac ? rix / m * (r + 1) + rix_res : rix / m;
      if constexpr (std::is_same_v<T, mer_t>) {
        table[rix].emplace_back(lshf->drop_ppos_lr(cminimizer.second), sh);
      } else {
        table[rix].push_back(lshf->drop_ppos_lr(cminimizer.second));
      }
      wnix++;
      if (cminimizer.first != pminimizer.first) {
        wcix++;
      }
      pminimizer = cminimizer;
    }
  }
}

QSeq::~QSeq()
{
  kseq_destroy(kseq);
  gzclose(gfile);
}

void QSeq::clear_curr_batch()
{
  seq_batch.clear();
  identifer_batch.clear();
}

QSeq::QSeq(std::string input)
{
  is_url = std::regex_match(input, url_regexp);
  if (is_url) {
#if defined _WLCURL && _WLCURL == 1
    input_path = download_url(input);
#else
    std::cerr << "Failed to download from URL, compiled without libcurl!" << std::endl;
#endif
  } else {
    input_path = input;
  }
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    std::cerr << "Failed to open thefile at " << input_path << std::endl;
    exit(EXIT_FAILURE);
  }
  kseq = kseq_init(gfile);
}

bool QSeq::read_next_batch()
{
  seq_batch.clear();
  identifer_batch.clear();
  seq_batch.reserve(rbatch_size);
  identifer_batch.reserve(rbatch_size);
  bool cont_reading = false;
  uint64_t ix = 0, bpc = 0;
  // Alternatively, use (ix < rbatch_size).
  while ((bpc < bpc_limit) && (cont_reading = kseq_read(kseq) >= 0)) {
    bpc += kseq->seq.l;
    seq_batch.emplace_back(kseq->seq.s);
    identifer_batch.emplace_back(kseq->name.s);
    ix++;
  }
  cbatch_size = ix;
  return cont_reading;
}

bool QSeq::is_batch_finished()
{
  assert(seq_batch.size() == identifer_batch.size());
  return seq_batch.empty() && identifer_batch.empty();
}

template void RSeq::extract_mers(vvec<mer_t>& table, sh_t sh);
template void RSeq::extract_mers(vvec<enc_t>& table, sh_t sh);
