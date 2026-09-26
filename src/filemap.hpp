#ifndef KREPP_FILEMAP_HPP
#define KREPP_FILEMAP_HPP

#include <cstddef>
#include <cstring>
#include <filesystem>
#include <utility>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace krepp {

  /* A read-only mapping of a whole file, used to view index arrays in place
 * instead of reading them into private buffers. The mapping keeps the data valid
 * until the FileMap dies; when the file cannot be opened or mapped the object
 * stays empty, which makes the loaders fall back to reading it. */
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
