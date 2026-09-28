#ifndef KREPP_FILEMAP_HPP
#define KREPP_FILEMAP_HPP

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string_view>
#include <utility>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace krepp {

  enum class MmapAdvice : uint8_t
  {
    none,
    random,
    willneed
  };

  inline MmapAdvice mmap_advice_from_env()
  {
    const char* value = std::getenv("KREPP_MMAP_ADVICE");
    if (value == nullptr) {
      return MmapAdvice::none;
    }
    const std::string_view name(value);
    if (name == "random") {
      return MmapAdvice::random;
    }
    if (name == "willneed") {
      return MmapAdvice::willneed;
    }
    return MmapAdvice::none;
  }

  class FileMap
  {
  public:
    FileMap() = default;

    explicit FileMap(const std::filesystem::path& path) { open(path); }
    ~FileMap() { close(); }

    FileMap(const FileMap&) = delete;
    FileMap& operator=(const FileMap&) = delete;
    FileMap(FileMap&& other) noexcept { swap(other); }
    FileMap& operator=(FileMap&& other) noexcept
    {
      swap(other);
      return *this;
    }

    void swap(FileMap& other) noexcept
    {
      std::swap(base_, other.base_);
      std::swap(size_, other.size_);
    }

    bool is_open() const { return base_ != nullptr; }
    const char* data() const { return base_; }
    size_t size() const { return size_; }

  private:
    void open(const std::filesystem::path& path)
    {
      close();
      const int fd = ::open(path.c_str(), O_RDONLY);
      if (fd < 0) return;
      struct stat st
      {};
      if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
        ::close(fd);
        return;
      }
      size_ = static_cast<size_t>(st.st_size);
      void* base = ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0);
      // The mapping outlives the descriptor, so it can be closed right away.
      ::close(fd);
      if (base == MAP_FAILED) {
        size_ = 0;
        return;
      }
      base_ = static_cast<const char*>(base);
      advise(base_, size_);
    }

    void advise(const void* base, size_t size)
    {
      switch (mmap_advice_from_env()) {
#if defined(MADV_RANDOM)
        case MmapAdvice::random:
          ::madvise(const_cast<void*>(base), size, MADV_RANDOM);
          break;
#endif
#if defined(MADV_WILLNEED)
        case MmapAdvice::willneed:
          ::madvise(const_cast<void*>(base), size, MADV_WILLNEED);
          break;
#endif
        default:
          break;
      }
    }

    void close()
    {
      if (base_ != nullptr) {
        ::munmap(const_cast<char*>(base_), size_);
        base_ = nullptr;
      }
      size_ = 0;
    }

    const char* base_ = nullptr;
    size_t size_ = 0;
  };

} // namespace krepp

#endif
