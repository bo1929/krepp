#include "sketch.hpp"

void Sketch::load()
{
  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  sflatht = std::make_shared<SFlatHT>();
  sflatht->load(sketch_stream);
  sketch_stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  sketch_stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  sketch_stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  sketch_stream.read(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  sketch_stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  sketch_stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
  sketch_stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  vec<uint8_t> ppos_v(h), npos_v(k - h);
  sketch_stream.read(reinterpret_cast<char*>(ppos_v.data()), ppos_v.size() * sizeof(uint8_t));
  sketch_stream.read(reinterpret_cast<char*>(npos_v.data()), npos_v.size() * sizeof(uint8_t));
  sketch_stream.read(reinterpret_cast<char*>(&rho), sizeof(double));

  if (!sketch_stream.good()) {
    std::cerr << "Reading the sketch file has failed!" << std::endl;
    exit(EXIT_FAILURE);
  }
  sketch_stream.close();

  lshf = std::make_shared<LSHF>(m, ppos_v, npos_v);
}

void Sketch::make_rho_partial()
{
  if (frac) {
    rho *= (static_cast<double>(r) + 1.0) / static_cast<double>(m);
  } else {
    rho *= 1.0 / static_cast<double>(m);
  }
}

std::vector<enc_t>::const_iterator Sketch::bucket_start() { return sflatht->bucket_start(offset); }

std::vector<enc_t>::const_iterator Sketch::bucket_next() { return sflatht->bucket_next(offset); }
