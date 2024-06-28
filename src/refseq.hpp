#ifndef _READER_H
#define _READER_H

#include "common.hpp"
#include "lshf.hpp"
#include "table.hpp"

extern "C" {
#include "kseq.h"
}

KSEQ_INIT(gzFile, gzread)

class RefSeq {
  friend class DynTable;

public:
  RefSeq(uint8_t k, uint8_t w, sh_t shash, std::string gpath,
         lshf_sptr_t hash_func);
  ~RefSeq();
  bool read_next_seq() { return kseq_read(kseq) >= 0; }
  bool set_curr_seq() {
    name = kseq->name.s;
    seq = kseq->seq.s;
    len = kseq->seq.l;
    return len >= k;
  }
  static void compute_encoding(char *s1, char *s2, uint64_t &enc_lr,
                               uint64_t &enc_bp) {
    enc_lr = 0;
    enc_bp = 0;
    for (; s1 < s2; s1++) {
      enc_lr <<= 1;
      enc_bp <<= 2;
      enc_bp += nt4_bp_table[seq_nt4_table[*s1]];
      enc_lr += nt4_lr_table[seq_nt4_table[*s1]];
    }
  }
  static void update_encoding(char *s1, uint64_t &enc_lr, uint64_t &enc_bp) {
    enc_lr <<= 1;
    enc_bp <<= 2;
    enc_lr &= 0xFFFFFFFEFFFFFFFE;
    enc_bp += nt4_bp_table[seq_nt4_table[*s1]];
    enc_lr += nt4_lr_table[seq_nt4_table[*s1]];
  }
  static size_t write_data(void *ptr, size_t size, size_t nmemb, FILE *fst) {
    size_t written = fwrite(ptr, size, nmemb, fst);
    return written;
  }
  void extract_mers(vvec<mer_t> &table);
  std::string download_url(std::string url);

private:
  uint8_t k;
  uint8_t w;
  lshf_sptr_t hash_func;
  std::string filepath;
  gzFile file;
  bool is_url;
  kseq_t *kseq;
  uint64_t len;
  char *seq;
  char *name;
  sh_t shash;
  uint64_t mask_bp;
  uint64_t mask_lr;
  const uint32_t max_rix = std::numeric_limits<uint32_t>::max();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const std::regex url_regexp = std::regex(
      R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");
};
#endif
