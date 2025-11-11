#ifndef _RQSEQ_H
#define _RQSEQ_H

#include "common.hpp"
#include "lshf.hpp"
#include "table.hpp"

/* #define CANONICAL */
#define RBATCH_SIZE 512
#define DSEQ_LEN 150

class HandlerURL
{
protected:
  const std::regex url_regexp = std::regex(
    R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
#if defined _WLCURL && _WLCURL == 1
  static size_t write_data(void* ptr, size_t s, size_t nmb, FILE* fst)
  {
    size_t written = fwrite(ptr, s, nmb, fst);
    return written;
  }

  std::string download_url(std::string url)
  {
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();
    if (!std::filesystem::exists(tmp_dir) || !std::filesystem::is_directory(tmp_dir)) {
      error_exit(std::string("Failed to get temp directory: ") + tmp_dir.string());
    }
    std::string hash_str = std::to_string(gp_hash(url));
    std::string tmp_filename = "rseq_" + hash_str + ".tmp";
    std::filesystem::path tmp_path = tmp_dir / tmp_filename;

    FILE* fp = fopen(tmp_path.string().c_str(), "wb");
    if (!fp) {
      error_exit(std::string("Failed to open temp file for writing: ") + tmp_path.string());
    }
    CURL* curl = curl_easy_init();
    if (!curl) {
      error_exit("Failed to initialize CURL.");
    }
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    CURLcode resb = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    if (resb != CURLE_OK) {
      error_exit(std::string("CURL download failed: ") + curl_easy_strerror(resb));
    }
    fclose(fp);

    return tmp_path.string();
  }
#endif
};

extern "C"
{
#include "kseq.h"
}

KSEQ_INIT(gzFile, gzread)

class RSeq : public HandlerURL
{
  friend class DynHT;

public:
  RSeq(std::string input, lshf_sptr_t lshf, uint8_t w, uint32_t r, bool frac, int sdust_t, int sdust_w);
  ~RSeq();
  bool read_next_seq() { return kseq_read(kseq) >= 0; }
  double get_rho() { return rho; }
  void compute_rho() { rho = static_cast<double>(wcix) / static_cast<double>(wnix); }
  bool set_curr_seq()
  {
    name = kseq->name.s;
    seq = kseq->seq.s;
    len = kseq->seq.l;
    return len >= w;
  }
  template<typename T>
  void extract_mers(vvec<T>& table, sh_t sh = 0);

private:
  uint8_t k;
  uint8_t w;
  uint32_t m;
  uint32_t r;
  char* seq;
  char* name;
  uint64_t len;
  bool frac;
  bool is_url;
  gzFile gfile;
  kseq_t* kseq;
  lshf_sptr_t lshf;
  uint64_t mask_bp = 0;
  uint64_t mask_lr = 0;
  double rho = 0;
  uint64_t wcix = 0;
  uint64_t wnix = 0;
  std::filesystem::path input_path = "";
  int sdust_t;
  int sdust_w;
};

class QSeq : public HandlerURL
{
  friend class IBatch;
  friend class SBatch;

public:
  QSeq(std::string input);
  ~QSeq();
  bool read_next_batch();
  bool is_batch_finished();
  void clear_curr_batch();
  uint64_t get_cbatch_size() { return cbatch_size; }

private:
  bool is_url;
  gzFile gfile;
  kseq_t* kseq;
  vec<std::string> seq_batch;
  vec<std::string> identifer_batch;
  std::filesystem::path input_path;
  uint64_t bpc_limit = RBATCH_SIZE * DSEQ_LEN;
  uint64_t rbatch_size = RBATCH_SIZE;
  uint64_t cbatch_size = 0;
};

#endif
