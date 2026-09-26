/* Unit tests for DynHT/SDynHT and their flattened FlatHT/SFlatHT forms. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

BuildOptions kConfig()
{
  BuildOptions opt;
  opt.k = 15;
  opt.w = 20;
  opt.h = 9;
  opt.m = 4;
  opt.r = 1;
  opt.frac = true;
  opt.seed = 5;
  return opt;
}

/* A tree with a known shape plus a matching Record; returned by value so the
 * tests can share them. */
struct Fixture
{
  tree_sptr_t tree;
  record_sptr_t record;
  uint32_t nrows;
};

Fixture make_fixture(const std::vector<std::string>& names, const BuildOptions& opt = kConfig())
{
  TestLSH lsh;
  lsh.configure(opt.k, opt.w, opt.h, opt.m, opt.r, opt.frac);
  gen.seed(opt.seed);
  lsh.build();

  auto tree = std::make_shared<Tree>();
  std::vector<std::string> mutable_names = names;
  tree->generate_tree(mutable_names);
  tree->reset_traversal();
  return {tree, std::make_shared<Record>(tree), lsh.rows()};
}

} // namespace

TEST_SUITE_BEGIN("table");

TEST_CASE("DynHT fill_table extracts the minimizers of a sequence")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("dynht");
  const std::string seq = rand_dna(3000, 11);
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", seq);

  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();

  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto dynht = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  node_sptr_t leaf_a = find_leaf(fx.tree, "A");
  REQUIRE(leaf_a != nullptr);
  dynht->fill_table(leaf_a->get_sh(), rs);

  // One entry per distinct minimizer of the sequence.
  CHECK(dynht->get_nkmers() > 0);
  CHECK(dynht->get_nkmers() < seq.size());

  // Converting to a flat table preserves the count and keeps the rows sorted
  // and unique.
  auto flatht = std::make_shared<FlatHT>(dynht);
  CHECK(flatht->get_nkmers() == dynht->get_nkmers());
  uint64_t total = 0;
  for (uint32_t rix = 0; rix < fx.nrows; ++rix) {
    auto begin = flatht->bucket_start(rix);
    auto end = flatht->bucket_next(rix);
    CHECK(begin <= end);
    total += static_cast<uint64_t>(std::distance(begin, end));
    // Rows are ordered by LSH hash, so only the entries inside one row are
    // sorted by encoding.
    CHECK(std::is_sorted(begin, end, [](const cmer_t& a, const cmer_t& b) { return a.first < b.first; }));
    for (auto it = begin; it != end; ++it) {
      // Every colour in a single-reference table is that reference.
      CHECK(it->second == leaf_a->get_se());
    }
  }
  CHECK(total == flatht->get_nkmers());
}

TEST_CASE("FlatHT save and load reproduce the same table")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("flatht");
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", rand_dna(2000, 21));

  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();

  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto dynht = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  dynht->fill_table(find_leaf(fx.tree, "A")->get_sh(), rs);
  auto flatht = std::make_shared<FlatHT>(dynht);

  const std::filesystem::path mer_path = dir / "cmer";
  const std::filesystem::path inc_path = dir / "inc";
  {
    std::ofstream mer_out(mer_path, std::ofstream::binary);
    std::ofstream inc_out(inc_path, std::ofstream::binary);
    flatht->save(mer_out, inc_out);
    CHECK(mer_out.good());
    CHECK(inc_out.good());
  }
  auto loaded = std::make_shared<FlatHT>(fx.tree, std::make_shared<CRecord>(fx.record));
  {
    std::ifstream mer_in(mer_path, std::ifstream::binary);
    std::ifstream inc_in(inc_path, std::ifstream::binary);
    loaded->load(mer_in, inc_in);
    CHECK(mer_in.good());
  }
  CHECK(loaded->get_nkmers() == flatht->get_nkmers());
  for (uint32_t rix = 0; rix < fx.nrows; ++rix) {
    auto a = flatht->bucket_start(rix);
    auto b = flatht->bucket_next(rix);
    auto c = loaded->bucket_start(rix);
    auto d = loaded->bucket_next(rix);
    REQUIRE(std::distance(a, b) == std::distance(c, d));
    CHECK(std::equal(a, b, c));
  }
}

TEST_CASE("SFlatHT save and load reproduce the same table")
{
  TempDir dir("sflatht");
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", rand_dna(2000, 33));

  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();

  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto sdynht = std::make_shared<SDynHT>();
  sdynht->fill_table(lsh.rows(), rs);
  auto sflatht = std::make_shared<SFlatHT>(sdynht);
  CHECK(sflatht->bucket_next(0) >= sflatht->bucket_start(0));

  const std::filesystem::path path = dir / "sketch";
  {
    std::ofstream out(path, std::ofstream::binary);
    sflatht->save(out);
    CHECK(out.good());
  }
  auto loaded = std::make_shared<SFlatHT>();
  {
    std::ifstream in(path, std::ifstream::binary);
    loaded->load(in);
    CHECK(in.good());
  }
  for (uint32_t rix = 0; rix < lsh.rows(); ++rix) {
    auto a = sflatht->bucket_start(rix);
    auto b = sflatht->bucket_next(rix);
    auto c = loaded->bucket_start(rix);
    auto d = loaded->bucket_next(rix);
    REQUIRE(std::distance(a, b) == std::distance(c, d));
    CHECK(std::equal(a, b, c));
  }
}

TEST_CASE("union_row merges two sorted rows and unions equal encodings")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  node_sptr_t a = find_leaf(fx.tree, "A");
  node_sptr_t b = find_leaf(fx.tree, "B");
  REQUIRE(a != nullptr);
  REQUIRE(b != nullptr);

  auto dynht = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  vec<mer_t> dest = {{10, a->get_sh()}, {30, a->get_sh()}, {50, a->get_sh()}};
  vec<mer_t> source = {{20, b->get_sh()}, {30, b->get_sh()}, {40, b->get_sh()}};
  dynht->union_row(dest, source);

  REQUIRE(dest.size() == 5);
  CHECK(dest[0].encoding == 10);
  CHECK(dest[0].sh == a->get_sh());
  CHECK(dest[1].encoding == 20);
  CHECK(dest[1].sh == b->get_sh());
  CHECK(dest[2].encoding == 30);
  CHECK(dest[3].encoding == 40);
  CHECK(dest[4].encoding == 50);
  CHECK(dest[4].sh == a->get_sh());

  // The shared encoding now carries the union of the two clades.
  CHECK(color_leaves_sh(fx.record, dest[2].sh) == std::vector<std::string>{"A", "B"});
}

TEST_CASE("union_row handles empty inputs in both directions")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  node_sptr_t a = find_leaf(fx.tree, "A");
  auto dynht = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);

  vec<mer_t> empty;
  vec<mer_t> row = {{1, a->get_sh()}, {2, a->get_sh()}};
  dynht->union_row(empty, row);
  CHECK(empty.size() == 2);

  vec<mer_t> row2 = {{3, a->get_sh()}};
  vec<mer_t> empty2;
  dynht->union_row(row2, empty2);
  CHECK(row2.size() == 1);
  CHECK(row2[0].encoding == 3);
}

TEST_CASE("union_table merges the k-mers of two references")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("union");
  const std::string shared = rand_dna(1200, 5);
  const std::string seq_a = shared + rand_dna(800, 6);
  const std::string seq_b = rand_dna(800, 7) + shared;
  const std::filesystem::path fa_a = dir / "a.fna";
  const std::filesystem::path fa_b = dir / "b.fna";
  write_fasta(fa_a, "A", seq_a);
  write_fasta(fa_b, "B", seq_b);

  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();

  node_sptr_t leaf_a = find_leaf(fx.tree, "A");
  node_sptr_t leaf_b = find_leaf(fx.tree, "B");
  auto table_a = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  auto table_b = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  {
    auto rs = std::make_shared<RSeq>(fa_a.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
    table_a->fill_table(leaf_a->get_sh(), rs);
  }
  {
    auto rs = std::make_shared<RSeq>(fa_b.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
    table_b->fill_table(leaf_b->get_sh(), rs);
  }
  const uint64_t nkmers_a = table_a->get_nkmers();
  const uint64_t nkmers_b = table_b->get_nkmers();
  REQUIRE(nkmers_a > 0);
  REQUIRE(nkmers_b > 0);

  table_a->union_table(table_b);
  CHECK(table_a->get_nkmers() <= nkmers_a + nkmers_b);
  CHECK(table_a->get_nkmers() >= std::max(nkmers_a, nkmers_b));

  // The shared region must produce colours covering both leaves.
  auto flatht = std::make_shared<FlatHT>(table_a);
  auto crecord = std::make_shared<CRecord>(fx.record);
  uint32_t two_leaf_colors = 0;
  for (uint32_t rix = 0; rix < fx.nrows; ++rix) {
    for (auto it = flatht->bucket_start(rix); it != flatht->bucket_next(rix); ++it) {
      const std::vector<std::string> leaves = color_leaves(crecord, it->second);
      if (leaves.size() == 2) {
        CHECK(leaves == std::vector<std::string>{"A", "B"});
        ++two_leaf_colors;
      } else {
        CHECK(leaves.size() == 1);
      }
    }
  }
  CHECK(two_leaf_colors > 0);
}

TEST_CASE("union_table with an empty source is a no-op")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("union-empty");
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", rand_dna(1500, 17));
  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();

  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto table = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  table->fill_table(find_leaf(fx.tree, "A")->get_sh(), rs);
  const uint64_t before = table->get_nkmers();

  auto empty = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  table->union_table(empty);
  CHECK(table->get_nkmers() == before);

  // Moving an empty table into a fresh one adopts the source wholesale.
  auto moved = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  moved->union_table(table);
  CHECK(moved->get_nkmers() == before);
  // The moved-from table reports nothing: its rows belong to the destination.
  CHECK(table->get_nkmers() == 0);
}

TEST_CASE("prune_columns caps the entries per row")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("prune");
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", rand_dna(4000, 23));
  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();
  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto table = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  table->fill_table(find_leaf(fx.tree, "A")->get_sh(), rs);
  const uint64_t before = table->get_nkmers();
  REQUIRE(before > 0);
  table->prune_columns(1);
  CHECK(table->get_nkmers() <= before);
  CHECK(table->get_nkmers() > 0);
  // Pruning to zero empties the table.
  table->prune_columns(0);
  CHECK(table->get_nkmers() == 0);
}

TEST_CASE("clear_rows and the size histogram agree on the table size")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  TempDir dir("hist");
  const std::filesystem::path fa = dir / "a.fna";
  write_fasta(fa, "A", rand_dna(2500, 29));
  const BuildOptions cfg = kConfig();
  TestLSH lsh;
  lsh.configure(cfg.k, cfg.w, cfg.h, cfg.m, cfg.r, cfg.frac);
  gen.seed(kConfig().seed);
  lsh.build();
  auto rs = std::make_shared<RSeq>(fa.string(), lsh.lsh(), lsh.win(), kConfig().r, kConfig().frac, 0, 0);
  auto table = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  table->fill_table(find_leaf(fx.tree, "A")->get_sh(), rs);
  const uint64_t nkmers = table->get_nkmers();
  REQUIRE(nkmers > 0);
  {
    CaptureStream quiet(std::cout);
    table->print_info();
    CHECK(quiet.str().find("size: " + std::to_string(nkmers)) != std::string::npos);
    CHECK(quiet.str().find("H(") != std::string::npos);
  }
  table->clear_rows();
  CHECK(table->get_nkmers() == 0);
}

TEST_CASE("conv_mer_cmer maps a colour handle onto its compact id")
{
  const Fixture fx = make_fixture({"A", "B", "C", "D"});
  // make_compact() runs in the CRecord constructor, so the ids are assigned
  // before anything asks for them.
  auto crecord = std::make_shared<CRecord>(fx.record);
  auto dynht = std::make_shared<DynHT>(fx.nrows, fx.tree, fx.record);
  node_sptr_t a = find_leaf(fx.tree, "A");
  const cmer_t cm = dynht->conv_mer_cmer(mer_t(42, a->get_sh()));
  CHECK(cm.first == 42);
  CHECK(cm.second == a->get_se());
  CHECK(color_leaves(crecord, cm.second) == std::vector<std::string>{"A"});
}

TEST_SUITE_END();
