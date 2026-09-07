#include "common.hpp"

uint32_t num_threads = 1;
std::string leave_out_ref = "";
std::string invocation = "";
thread_local std::random_device rd;
thread_local std::mt19937 gen;
uint32_t seed = 0;

const unsigned char seq_nt4_table[128] = { // Table to change "ACGTN" to 01234
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

const uint64_t nt4_lr_table[4] = {0, 1, 4294967296, 4294967297};

const uint64_t nt4_bp_table[4] = {0, 1, 2, 3};

// Function-local static: avoids depending on initialization order between
// translation units, since error_exit may run before namespace-scope objects
// in other units are constructed.
static error_handler_t& error_handler()
{
  static error_handler_t handler;
  return handler;
}

void set_error_handler(error_handler_t handler) { error_handler() = std::move(handler); }

[[noreturn]] void error_exit(const std::string& msg, int code)
{
  // Deliberately a reference, not a copy: see the note on set_error_handler
  // about not replacing the handler once krepp is running.
  const error_handler_t& handler = error_handler();
  if (handler) {
    handler(msg, code);
    // The handler returned, which its contract forbids. Report that rather than
    // falling back silently, since the fallback is indistinguishable from
    // having installed no handler at all.
    std::cerr << "[ERROR] error handler returned instead of throwing or exiting; original error (" << code << "): " << msg
              << std::endl;
    std::abort();
  }
  std::cerr << "[ERROR] " << msg << std::endl;
  std::exit(code);
}

inline void warn_msg(const std::string& msg) { std::cerr << "[WARNING] " << msg << std::endl; }

template<typename StreamT>
inline void check_fstream_or_exit(const StreamT& stream, const std::string& msg, const std::string& path)
{
  if (!stream.good()) {
    if (!path.empty()) {
      error_exit(msg + ": " + path);
    } else {
      error_exit(msg);
    }
  }
}

template void check_fstream_or_exit<std::ifstream>(const std::ifstream&, const std::string&, const std::string&);
template void check_fstream_or_exit<std::ofstream>(const std::ofstream&, const std::string&, const std::string&);
