#include "rqseq.hpp"

RSeq::RSeq(uint8_t w, uint32_t r, bool frac, sh_t shash, lshf_sptr_t lshashf, std::string input)
  : w(w)
  , r(r)
  , frac(frac)
  , shash(shash)
  , lshashf(lshashf)
{
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  k = lshashf->get_k();
  m = lshashf->get_m();
  mask_bp = u64m >> (32 - k) * 2;
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  is_url = std::regex_match(input, url_regexp);
  if (is_url) {
    input_path = download_url(input);
  } else {
    input_path = input;
  }
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    std::cerr << "Failed to open gfile at " << input_path << std::endl;
    exit(1);
  }
  kseq = kseq_init(gfile);
}

std::string RSeq::download_url(std::string url)
{
  char tmp_path[FILENAME_MAX] = "/tmp/seq";
  const char* sx = std::to_string(gp_hash(url)).c_str();
  strcat(tmp_path, sx);
  strcat(tmp_path, ".XXXXXX");
  int tmp_fd = mkstemp(tmp_path);
  CURL* curl;
  FILE* fp;
  CURLcode resb;
  curl = curl_easy_init();
  if (curl) {
    fp = fopen(tmp_path, "wb");
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    resb = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    fclose(fp);
  }
  return tmp_path;
}

RSeq::~RSeq()
{
  kseq_destroy(kseq);
  gzclose(gfile);
  if (is_url) {
    std::filesystem::remove(input_path);
  }
}

void RSeq::extract_mers(vvec<mer_t>& table)
{
  uint32_t i, l;
  uint32_t rix_res;
  uint64_t wix = 0, kix = 0;
  uint8_t ldiff = w - k + 1;
  uint64_t enc64_bp, enc64_lr, orenc64_bp, orenc64_lr, rcenc64_bp;
  std::vector<std::pair<uint32_t, enc_t>> lsh_enc_win(ldiff);
  std::pair<uint32_t, enc_t> minimizer;
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
    if (orenc64_bp < rcenc64_bp) {
      orenc64_bp = rcenc64_bp;
      orenc64_lr = conv_bp64_lr64(rcenc64_bp);
    }
    lsh_enc_win[kix % ldiff].first = lshashf->compute_hash(orenc64_bp);
    lsh_enc_win[kix % ldiff].second = lshashf->drop_ppos_lr(orenc64_lr);
    if (l >= w || (i == len)) {
      minimizer =
        *std::min_element(lsh_enc_win.begin(),
                          lsh_enc_win.end(),
                          [](std::pair<uint32_t, enc_t> lhs, std::pair<uint32_t, enc_t> rhs) {
                            return xur32_hash(lhs.second) < xur32_hash(rhs.second);
                          });
      rix_res = minimizer.first % m;
      if (frac ? rix_res <= r : rix_res == r) {
        minimizer.first = (minimizer.first / m);
        assert(minimizer.first < table.size());
        table[minimizer.first].emplace_back(minimizer.second, shash);
        wix++;
      }
    }
    kix++;
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
  name_batch.clear();
}

QSeq::QSeq(std::filesystem::path input_path)
  : input_path(input_path)
  , batch_size(0)
{
  gfile = gzopen(input_path.c_str(), "rb");
  if (gfile == nullptr) {
    std::cerr << "Failed to open gfile at " << input_path << std::endl;
    exit(1);
  }
  kseq = kseq_init(gfile);
}

bool QSeq::read_next_batch()
{
  seq_batch.reserve(BATCH_SIZE);
  name_batch.reserve(BATCH_SIZE);
  bool cont_reading = false;
  uint64_t ix = 0;
  while ((cont_reading = kseq_read(kseq) >= 0) && ix < BATCH_SIZE) {
    seq_batch.emplace_back(kseq->seq.s);
    name_batch.emplace_back(kseq->name.s);
    ix++;
  }
  batch_size = ix;
  return cont_reading;
}

bool QSeq::is_batch_finished()
{
  assert(seq_batch.size() == name_batch.size());
  return seq_batch.empty() && name_batch.empty();
}