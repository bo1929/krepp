#ifndef _RQSEQ_H
#define _RQSEQ_H

#include "common.hpp"
#include "lshf.hpp"
#include "table.hpp"

#define BATCH_SIZE 262144

extern "C"
{
#include "kseq.h"
}

KSEQ_INIT(gzFile, gzread)

class RSeq
{
  friend class DynHT;

public:
  RSeq(uint8_t w, uint32_t r, bool frac, sh_t shash, lshf_sptr_t lshashf, std::string input);
  ~RSeq();
  bool read_next_seq() { return kseq_read(kseq) >= 0; }
  bool set_curr_seq()
  {
    name = kseq->name.s;
    seq = kseq->seq.s;
    len = kseq->seq.l;
    return len >= k;
  }
  static size_t write_data(void* ptr, size_t s, size_t nmb, FILE* fst)
  {
    size_t written = fwrite(ptr, s, nmb, fst);
    return written;
  }
  void extract_mers(vvec<mer_t>& table);
  std::string download_url(std::string url);

private:
  uint8_t k;
  uint8_t w;
  uint32_t m;
  uint32_t r;
  bool frac;
  lshf_sptr_t lshashf;
  std::filesystem::path input_path;
  gzFile gfile;
  bool is_url;
  kseq_t* kseq;
  uint64_t len;
  char* seq;
  char* name;
  sh_t shash;
  uint64_t mask_bp;
  uint64_t mask_lr;
  const std::regex url_regexp = std::regex(
    R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
};

class QSeq
{
  friend class QBatch;

public:
  QSeq(std::filesystem::path input_path);
  ~QSeq();
  bool read_next_batch();
  bool is_batch_finished();
  void clear_curr_batch();

private:
  gzFile gfile;
  kseq_t* kseq;
  uint64_t batch_size;
  vec<std::string> seq_batch;
  vec<std::string> name_batch;
  std::filesystem::path input_path;
};

#endif
