#include "reader.hpp"

Reader::Reader(uint8_t k, uint8_t w, sh_t shash, std::string path,
               lshf_sptr_t hash_func)
    : k(k), w(w), shash(shash), hash_func(hash_func) {
  is_url = std::regex_match(path, url_regexp);
  if (is_url) {
    filepath = download_url(path);
  } else {
    filepath = path;
  }
  file = gzopen(filepath.c_str(), "rb");
  if (file == nullptr) {
    std::cerr << "Failed to open file at " << filepath << std::endl;
    exit(1);
  }
  kseq = kseq_init(file);

  mask_bp = u64m >> (32 - k) * 2;
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
}

std::string Reader::download_url(std::string url) {
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

Reader::~Reader() {
  kseq_destroy(kseq);
  gzclose(file);
  if (is_url)
    std::remove(filepath.c_str());
}

bool Reader::make_ready() {
  name = kseq->name.s;
  seq = kseq->seq.s;
  len = kseq->seq.l;
  return len >= k;
}

void Reader::extract_mers(vvec<mer_t> &table) {
  uint64_t i, l;
  uint32_t rix;
  uint64_t wix = 0;
  uint64_t kix = 0;
  uint8_t ldiff = w - k + 1;
  uint32_t cenc32_lr, cenc32_bp;
  uint64_t enc64_bp, enc64_lr, cenc64_bp, cenc64_lr;
  uint64_t minimizer;
  mer_t curr_mer;
  curr_mer.shash = shash;
  uint64_t len_seq = len >= k ? len : 0;
  std::vector<uint64_t> lsh_enc_win(ldiff);
  for (i = l = 0; i < len_seq; ++i) {
    if (seq_nt4_table[seq[i]] < 4) { // not an "N" base
      l++;
      if (l == k) {
        compute_encoding(seq + i + 1 - k, seq + i + 1, enc64_lr, enc64_bp);
      } else if (l > k) {
        update_encoding(seq[i], enc64_lr, enc64_bp);
      }
      if (l >= k) {
        if (cenc64_bp < bp64revcomp(cenc64_bp, k)) {
          cenc64_bp = bp64revcomp(cenc64_bp, k);
          cenc64_lr = bp64convlr64(cenc64_bp);
        }
        rix = hash_func->compute_hash(cenc64_bp);
        hash_func->drop_ppos(cenc64_bp, cenc64_lr, cenc32_bp, cenc32_lr);
        assert(rix <= max_rix);
        if (ldiff > 1) {
          lsh_enc_win[kix % ldiff] =
              static_cast<uint64_t>(rix) |
              ((static_cast<uint64_t>(cenc32_lr) | 0) << 32);
          kix++;
        } else {
          curr_mer.encoding = cenc32_lr;
          table[rix].push_back(curr_mer);
          wix++;
        }
      }
      if ((l >= w || ((i == len_seq - 1) && l >= k)) && ldiff > 1) {
        minimizer = *std::min_element(lsh_enc_win.begin(), lsh_enc_win.end(),
                                      [](uint64_t lhs, uint64_t rhs) {
                                        return mmh3f64(lhs) < mmh3f64(rhs);
                                        /* return lhs < rhs; */
                                      });

        curr_mer.encoding = (minimizer >> 32) & 0xFFFFFFFF;
        table[minimizer & 0xFFFFFFFF].push_back(curr_mer);
        wix++;
      }
    } else {
      l = 0;
    }
  }
}
