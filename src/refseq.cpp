#include "refseq.hpp"
#include "common.hpp"

RefSeq::RefSeq(uint8_t k, uint8_t w, sh_t shash, std::string genomepath, lshf_sptr_t hash_func)
  : k(k)
  , w(w)
  , shash(shash)
  , hash_func(hash_func)
{
  mask_bp = u64m >> (32 - k) * 2;
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  is_url = std::regex_match(genomepath, url_regexp);
  if (is_url) {
    filepath = download_url(genomepath);
  } else {
    filepath = genomepath;
  }
  file = gzopen(filepath.c_str(), "rb");
  if (file == nullptr) {
    std::cerr << "Failed to open file at " << filepath << std::endl;
    exit(1);
  }
  kseq = kseq_init(file);
}

std::string RefSeq::download_url(std::string url)
{
  char tmp_input_path[FILENAME_MAX] = "/tmp/seq";
  const char* sx = std::to_string(gp_hash(url)).c_str();
  strcat(tmp_input_path, sx);
  strcat(tmp_input_path, ".XXXXXX");
  int tmp_fd = mkstemp(tmp_input_path);
  CURL* curl;
  FILE* fp;
  CURLcode resb;
  curl = curl_easy_init();
  if (curl) {
    fp = fopen(tmp_input_path, "wb");
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    resb = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    fclose(fp);
  }
  return tmp_input_path;
}

RefSeq::~RefSeq()
{
  kseq_destroy(kseq);
  gzclose(file);
  if (is_url) {
    std::remove(filepath.c_str());
  }
}

void RefSeq::extract_mers(vvec<mer_t>& table)
{
  uint64_t i, l;
  uint32_t rix, rix_res;
  uint64_t wix = 0, kix = 0;
  uint8_t ldiff = w - k + 1;
  uint32_t cenc32_lr, cenc32_bp;
  uint64_t enc64_bp, enc64_lr, cenc64_bp, cenc64_lr, renc64_bp;
  std::vector<std::pair<uint32_t, enc_t>> lsh_enc_win(ldiff);
  std::pair<uint32_t, enc_t> minimizer;
  uint64_t len_seq = len >= k ? len : 0;
  for (i = l = 0; i < len_seq;) {
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
    cenc64_bp = enc64_bp & mask_bp;
    cenc64_lr = enc64_lr & mask_lr;
    renc64_bp = revcomp_b64(cenc64_bp, k);
    if (cenc64_bp < renc64_bp) {
      cenc64_bp = renc64_bp;
      cenc64_lr = cast_bp64_lr64(renc64_bp);
    }
    hash_func->drop_ppos_enc(cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
    rix = hash_func->compute_hash(cenc64_bp);
    lsh_enc_win[kix % ldiff].first = rix;
    lsh_enc_win[kix % ldiff].second = cenc32_lr;
    if (l >= w || (i == len_seq)) {
      minimizer =
        *std::min_element(lsh_enc_win.begin(),
                          lsh_enc_win.end(),
                          [](std::pair<uint32_t, enc_t> lhs, std::pair<uint32_t, enc_t> rhs) {
                            return xur32_hash(lhs.second) < xur32_hash(rhs.second);
                          });
      rix_res = minimizer.first % hash_func->m;
      if (hash_func->frac ? rix_res <= hash_func->r : rix_res == hash_func->r) {
        minimizer.first = (minimizer.first / hash_func->m);
        assert(minimizer.first < table.size());
        table[minimizer.first].emplace_back(minimizer.second, shash);
        wix++;
      }
    }
    kix++;
  }
}
