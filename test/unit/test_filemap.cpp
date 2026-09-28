/* Unit tests for the read-only file mapping used by the index loaders. */

#include "test_helpers.hpp"

using namespace ktest;

TEST_SUITE_BEGIN("filemap");

TEST_CASE("the mapping advice is read from the environment")
{
  unsetenv("KREPP_MMAP_ADVICE");
  CHECK(krepp::mmap_advice_from_env() == krepp::MmapAdvice::none);
  setenv("KREPP_MMAP_ADVICE", "random", 1);
  CHECK(krepp::mmap_advice_from_env() == krepp::MmapAdvice::random);
  setenv("KREPP_MMAP_ADVICE", "willneed", 1);
  CHECK(krepp::mmap_advice_from_env() == krepp::MmapAdvice::willneed);
  setenv("KREPP_MMAP_ADVICE", "nonsense", 1);
  CHECK(krepp::mmap_advice_from_env() == krepp::MmapAdvice::none);
  unsetenv("KREPP_MMAP_ADVICE");
  CHECK(krepp::mmap_advice_from_env() == krepp::MmapAdvice::none);
}

TEST_CASE("a mapping views the whole file and survives being moved")
{
  TempDir dir("filemap-basic");
  const std::string data = "0123456789abcdef";
  const std::filesystem::path path = dir / "data.bin";
  spit(path, data);

  krepp::FileMap map(path);
  REQUIRE(map.is_open());
  CHECK(map.size() == data.size());
  CHECK(std::string(map.data(), map.size()) == data);

  krepp::FileMap moved(std::move(map));
  CHECK(moved.is_open());
  CHECK_FALSE(map.is_open());
  CHECK(std::string(moved.data(), moved.size()) == data);

  CHECK_FALSE(krepp::FileMap(dir / "does-not-exist.bin").is_open());
}

TEST_SUITE_END();
