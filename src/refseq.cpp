#include "refseq.hpp"

RefSeq::RefSeq(uint8_t k, uint8_t w, sh_t shash, std::string gpath,
               lshf_sptr_t hash_func)
    : k(k), w(w), shash(shash), hash_func(hash_func) {
  mask_bp = u64m >> (32 - k) * 2;
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  is_url = std::regex_match(gpath, url_regexp);
  if (is_url) {
    filepath = download_url(gpath);
  } else {
    filepath = gpath;
  }
  file = gzopen(filepath.c_str(), "rb");
  if (file == nullptr) {
    std::cerr << "Failed to open file at " << filepath << std::endl;
    exit(1);
  }
  kseq = kseq_init(file);
}

std::string RefSeq::download_url(std::string url) {
  char tmp_input_path[FILENAME_MAX] = "/tmp/seq";
  const char *sx = std::to_string(gp_hash(url)).c_str();
  strcat(tmp_input_path, sx);
  strcat(tmp_input_path, ".XXXXXX");
  int tmp_fd = mkstemp(tmp_input_path);
  CURL *curl;
  FILE *fp;
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

RefSeq::~RefSeq() {
  kseq_destroy(kseq);
  gzclose(file);
  if (is_url)
    std::remove(filepath.c_str());
}

void RefSeq::extract_mers(vvec<mer_t> &table) {
  uint64_t i, l;
  uint32_t rix;
  uint64_t wix = 0, kix = 0;
  uint8_t ldiff = w - k + 1;
  uint32_t cenc32_lr, cenc32_bp;
  uint64_t enc64_bp, enc64_lr, cenc64_bp, cenc64_lr;
  mer_t curr_mer;
  curr_mer.shash = shash;
  uint64_t len_seq = len >= k ? len : 0;
  std::pair<uint32_t, uint32_t> minimizer;
  std::vector<std::pair<uint32_t, uint32_t>> lsh_enc_win(ldiff);
  for (i = l = 0; i < len_seq;) {
    if (seq_nt4_table[seq[i]] < 4) {
      l++, i++;
      if (l >= k) {
        if (l == k) {
          compute_encoding(seq + i - k, seq, enc64_lr, enc64_bp);
        } else {
          update_encoding(seq + i - 1, enc64_lr, enc64_bp);
        }
        if (cenc64_bp < bp64revcomp(cenc64_bp, k)) {
          cenc64_bp = bp64revcomp(cenc64_bp, k);
          cenc64_lr = bp64convlr64(cenc64_bp);
        }
        rix = hash_func->compute_hash(cenc64_bp);
        hash_func->drop_ppos_encoding(cenc64_bp, cenc64_lr, cenc32_bp,
                                      cenc32_lr);
        assert(rix <= max_rix);
        if (ldiff > 1) {
          lsh_enc_win[kix % ldiff].first = rix;
          lsh_enc_win[kix % ldiff].second = cenc32_lr;
          kix++;
          if (l >= w || (i == len_seq)) {
            minimizer = *std::min_element(
                lsh_enc_win.begin(), lsh_enc_win.end(), [](auto lhs, auto rhs) {
                  return mur_hash((static_cast<uint64_t>(lhs.first) << 32) |
                                  static_cast<uint64_t>(lhs.second)) <
                         mur_hash((static_cast<uint64_t>(rhs.first) << 32) |
                                  static_cast<uint64_t>(rhs.second));
                });
            curr_mer.encoding = minimizer.second;
            table[minimizer.first].push_back(curr_mer);
            wix++;
          }
        } else {
          curr_mer.encoding = cenc32_lr;
          table[rix].push_back(curr_mer);
          wix++;
        }
      }
    } else {
      l = 0, i++;
    }
  }
}
