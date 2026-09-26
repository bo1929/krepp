/* Unit tests for Record (build time clades) and CRecord (query time colours). */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

const std::vector<std::string> kNames = {"A", "B", "C", "D", "E"};

} // namespace

TEST_SUITE_BEGIN("record");

TEST_CASE("Record indexes every node of the tree")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  const std::vector<node_sptr_t> nodes = post_order(tree);
  CHECK(nodes.size() == 2 * kNames.size() - 1);
  CHECK(record->get_size() >= nodes.size());

  // decode_sh() of a node hash yields that node's clade.
  for (const node_sptr_t& nd : nodes) {
    if (!nd->get_sh()) continue;
    CHECK(color_leaves_sh(record, nd->get_sh()) == clade_leaves(nd));
  }
}

TEST_CASE("add_subset returns a clade covering both operands")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t b = find_leaf(tree, "B");
  node_sptr_t c = find_leaf(tree, "C");
  node_sptr_t d = find_leaf(tree, "D");
  REQUIRE(a != nullptr);
  REQUIRE(b != nullptr);
  REQUIRE(c != nullptr);
  REQUIRE(d != nullptr);

  const std::vector<std::string> ab = color_leaves_sh(record, record->add_subset(a->get_sh(), b->get_sh()));
  CHECK(ab == std::vector<std::string>{"A", "B"});

  const std::vector<std::string> cd = color_leaves_sh(record, record->add_subset(c->get_sh(), d->get_sh()));
  CHECK(cd == std::vector<std::string>{"C", "D"});

  // The union of the two unions is the whole tree.
  const std::vector<std::string> abcd = color_leaves_sh(record, record->add_subset(a->get_sh(), c->get_sh()));
  CHECK(abcd == std::vector<std::string>{"A", "C"});

  // Adding the same pair twice must be stable.
  const sh_t first = record->add_subset(a->get_sh(), b->get_sh());
  const sh_t second = record->add_subset(a->get_sh(), b->get_sh());
  CHECK(first == second);
  CHECK(color_leaves_sh(record, second) == ab);
}

TEST_CASE("add_subset rejects unknown handles")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  ThrowingErrorHandler handler;
  const std::string msg = ThrowingErrorHandler::catches([&] { record->add_subset(0xdeadbeefull, 0xfeedfaceull); });
  CHECK(msg.find("Failed for partition") != std::string::npos);
}

TEST_CASE("make_compact assigns a unique id to every node and colour")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);
  CHECK(crecord->get_nsubsets() > 0);
  CHECK(crecord->get_nsubsets() <= record->get_size() + 1);
}

TEST_CASE("CRecord decodes a tree node into its clade")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);
  const std::vector<node_sptr_t> nodes = post_order(tree);
  for (const node_sptr_t& nd : nodes) {
    if (!nd->get_se()) continue;
    CHECK(color_leaves(crecord, nd->get_se()) == clade_leaves(nd));
  }
}

TEST_CASE("every colour decomposes into two smaller colours")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  // Add a merged colour so the subset lattice is more than just the nodes.
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t c = find_leaf(tree, "C");
  const sh_t merged_sh = record->add_subset(a->get_sh(), c->get_sh());
  auto crecord = std::make_shared<CRecord>(record);

  const se_t nsubsets = crecord->get_nsubsets();
  se_t internal = 0, leaves = 0, merged = 0, empty = 0;
  for (se_t se = 1; se < nsubsets; ++se) {
    const pse_t pse = crecord->get_pse(se);
    REQUIRE(pse.first < nsubsets);
    REQUIRE(pse.second < nsubsets);
    node_sptr_t nd = tree->check_node(se) ? tree->get_node(se) : nullptr;
    if (pse.first == 0 && pse.second == 0) {
      // Either a leaf or the reserved empty colour (Record stores a Subset for
      // sh == 0, which make_compact() also gives an id).
      if (nd && nd->check_leaf()) {
        ++leaves;
      } else {
        CHECK(nd == nullptr);
        ++empty;
      }
      continue;
    }
    // Every other colour is the disjoint union of its two parts.
    CHECK(pse.first != se);
    CHECK(pse.second != se);
    std::vector<std::string> parts = color_leaves(crecord, pse.first);
    const std::vector<std::string> other = color_leaves(crecord, pse.second);
    parts.insert(parts.end(), other.begin(), other.end());
    std::sort(parts.begin(), parts.end());
    parts.erase(std::unique(parts.begin(), parts.end()), parts.end());
    CHECK(parts == color_leaves(crecord, se));
    if (nd) {
      ++internal;
      CHECK_FALSE(nd->check_leaf());
    } else {
      ++merged;
    }
  }
  CHECK(leaves == kNames.size());
  CHECK(internal > 0);
  CHECK(merged > 0);
  CHECK(empty == 1);

  // The merged colour covers exactly the two leaves it was built from.
  vec<node_sptr_t> decoded;
  record->decode_sh(merged_sh, decoded);
  CHECK(color_leaves_sh(record, merged_sh) == std::vector<std::string>{"A", "C"});
}

TEST_CASE("CRecord round-trips through save and load")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  record->add_subset(find_leaf(tree, "A")->get_sh(), find_leaf(tree, "C")->get_sh());
  record->insert_rho(find_leaf(tree, "A")->get_sh(), 0.25);
  record->insert_rho(find_leaf(tree, "B")->get_sh(), 0.75);
  auto crecord = std::make_shared<CRecord>(record);

  TempDir dir("crecord");
  const std::filesystem::path path = dir / "crecord";
  {
    std::ofstream out(path, std::ofstream::binary);
    crecord->save(out);
    CHECK(out.good());
  }
  auto loaded = std::make_shared<CRecord>(tree);
  {
    std::ifstream in(path, std::ifstream::binary);
    loaded->load(in);
    CHECK(in.good());
  }
  CHECK(loaded->get_nsubsets() == crecord->get_nsubsets());
  for (se_t se = 0; se < crecord->get_nsubsets(); ++se) {
    CHECK(loaded->get_pse(se) == crecord->get_pse(se));
  }
  // The rho values travel per node, in node id order.
  for (const node_sptr_t& nd : post_order(tree)) {
    CHECK(loaded->get_rho(nd->get_se()) == doctest::Approx(crecord->get_rho(nd->get_se())));
  }
  CHECK(loaded->get_rho(find_leaf(tree, "A")->get_se()) == doctest::Approx(0.25));
  CHECK(loaded->get_rho(find_leaf(tree, "B")->get_se()) == doctest::Approx(0.75));
}

TEST_CASE("apply_rho_coef scales every stored rho")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  for (const node_sptr_t& nd : post_order(tree)) {
    record->insert_rho(nd->get_sh(), 0.5);
  }
  auto crecord = std::make_shared<CRecord>(record);
  const double before = crecord->get_rho(1);
  crecord->apply_rho_coef(0.5);
  CHECK(crecord->get_rho(1) == doctest::Approx(before * 0.5));
}

TEST_CASE("decoding the reserved empty colour terminates")
{
  // 0 is the reserved empty colour: it has no decomposition, so a walk that
  // reaches it has to stop instead of pushing it again forever.
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);

  vec<node_sptr_t> from_sh;
  record->decode_sh(0, from_sh);
  CHECK(from_sh.empty());

  vec<node_sptr_t> from_se;
  crecord->decode_se(0, from_se);
  CHECK(from_se.empty());

  // A leaf decomposes into the empty colour, which is the path that used to
  // loop: walking a leaf's own id must still terminate.
  node_sptr_t a = find_leaf(tree, "A");
  vec<node_sptr_t> leaves;
  crecord->decode_se(a->get_se(), leaves);
  REQUIRE(leaves.size() == 1);
  CHECK(leaves.front() == a);
}

TEST_CASE("a CRecord built from a plain tree can decode every node")
{
  tree_sptr_t tree = make_tree(kNames);
  // This constructor leaves the colours empty on purpose: every node is its own
  // colour and there is no decomposition to follow.
  auto crecord = std::make_shared<CRecord>(tree);
  for (const node_sptr_t& nd : post_order(tree)) {
    if (!nd->get_se()) continue;
    CHECK(color_leaves(crecord, nd->get_se()) == clade_leaves(nd));
  }
  CHECK(crecord->get_pse(1) == std::make_pair(se_t{0}, se_t{0}));
}

TEST_CASE("decoding an out-of-range colour is reported")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);
  ThrowingErrorHandler handler;
  vec<node_sptr_t> out;
  const std::string msg = ThrowingErrorHandler::catches([&] { crecord->decode_se(crecord->get_nsubsets() + 100, out); });
  CHECK(msg.find("Invalid colour ID") != std::string::npos);
}

TEST_CASE("CRecord accepts a plain tree without a Record")
{
  tree_sptr_t tree = make_tree(kNames);
  auto crecord = std::make_shared<CRecord>(tree);
  CHECK(crecord->get_nsubsets() == tree->get_nnodes() + 1);
  for (se_t se = 1; se <= tree->get_nnodes(); ++se) {
    CHECK(crecord->get_rho(se) == doctest::Approx(0.0));
  }
}

TEST_CASE("Record::print_info reports the number of colours and nodes")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);
  CaptureStream quiet(std::cout);
  crecord->print_info();
  CHECK(quiet.str().find("Total number of subsets") != std::string::npos);
  CHECK(quiet.str().find("Number of nodes") != std::string::npos);
}

TEST_CASE("CRecord::display_info reports colour statistics")
{
  tree_sptr_t tree = make_tree(kNames);
  auto record = std::make_shared<Record>(tree);
  auto crecord = std::make_shared<CRecord>(record);
  std::stringstream out;
  vec<uint64_t> se_to_count(crecord->get_nsubsets(), 0);
  se_to_count[1] = 3;
  crecord->display_info(&out, 0, se_to_count);
  CHECK(out.str().find("NUM_COLORS") != std::string::npos);
  CHECK(out.str().find("MER_COUNT") != std::string::npos);
  CHECK(out.str().find("OUTDEGREE_COUNT") != std::string::npos);
}

TEST_CASE("a tree with a single-child node is rejected")
{
  ThrowingErrorHandler handler;
  for (const std::string& nwk_str : {std::string("(A:1);"), std::string("((A:1):1,B:1);"), std::string("(A:1,(B:1):1);")}) {
    auto tree = std::make_shared<Tree>();
    std::stringstream nwk(nwk_str);
    const std::string msg = ThrowingErrorHandler::catches([&] { tree->load(nwk); });
    CHECK(msg.find("single child") != std::string::npos);
  }
}

TEST_SUITE_END();
