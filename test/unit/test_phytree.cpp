/* Unit tests for Newick parsing, traversal, jplace decorations and lineages. */

#include "test_helpers.hpp"

using namespace ktest;

namespace {

tree_sptr_t load_tree(const std::string& nwk)
{
  auto tree = std::make_shared<Tree>();
  std::stringstream ss(nwk);
  tree->load(ss);
  tree->reset_traversal();
  return tree;
}

std::string print_basic(tree_sptr_t tree)
{
  strstream out;
  tree->stream_nwk_basic(out, tree->get_root());
  return out.str();
}

std::string print_jplace(tree_sptr_t tree)
{
  strstream out;
  tree->stream_nwk_jplace(out, tree->get_root());
  return out.str();
}

} // namespace

TEST_SUITE_BEGIN("phytree");

TEST_CASE("generate_tree builds a balanced tree over the given names")
{
  const std::vector<std::string> names = {"A", "B", "C", "D", "E"};
  tree_sptr_t tree = make_tree(names);
  CHECK(tree->get_nnodes() == 2 * names.size() - 1);
  CHECK(tree->get_root()->get_card() == names.size());

  std::vector<std::string> leaves;
  for (const node_sptr_t& nd : post_order(tree)) {
    if (nd->check_leaf()) {
      leaves.push_back(nd->get_name());
      CHECK(nd->get_card() == 1);
      CHECK(nd->get_parent() != nullptr);
    } else {
      CHECK(nd->get_nchildren() >= 2);
      CHECK(nd->get_card() > 1);
    }
    REQUIRE(nd->get_sh() != 0);
    CHECK(nd->get_se() >= 1);
    CHECK(nd->get_se() <= tree->get_nnodes());
  }
  std::sort(leaves.begin(), leaves.end());
  CHECK(leaves == std::vector<std::string>{"A", "B", "C", "D", "E"});
}

TEST_CASE("post-order visits every child before its parent")
{
  tree_sptr_t tree = make_tree({"A", "B", "C", "D", "E", "F", "G"});
  std::vector<node_sptr_t> order = post_order(tree);
  flat_phmap<node_sptr_t, size_t> position;
  for (size_t i = 0; i < order.size(); ++i) position[order[i]] = i;
  CHECK(order.size() == tree->get_nnodes());
  for (const node_sptr_t& nd : order) {
    for (tuint_t i = 0; i < nd->get_nchildren(); ++i) {
      node_sptr_t child = *std::next(nd->get_children(), i);
      CHECK(position[child] < position[nd]);
      CHECK(child->get_parent() == nd);
    }
  }
  // The root comes last.
  CHECK(order.back() == tree->get_root());
}

TEST_CASE("a Newick tree round-trips through load and stream_nwk_basic")
{
  const std::string nwk = "((A:1,B:2)N1:3,C:4)root;";
  tree_sptr_t tree = load_tree(nwk);
  CHECK(tree->get_nnodes() == 5);
  CHECK(tree->get_root()->get_name() == "root");
  const std::string printed = print_basic(tree);
  CHECK(printed.find("A:1") != std::string::npos);
  CHECK(printed.find("B:2") != std::string::npos);
  CHECK(printed.find("C:4") != std::string::npos);
  CHECK(printed.back() == ';');
  // Reloading the printed tree gives the same topology and labels.
  tree_sptr_t again = load_tree(printed);
  CHECK(again->get_nnodes() == tree->get_nnodes());
  CHECK(print_basic(again) == printed);

  // save() writes back the original text.
  std::stringstream saved;
  tree->save(saved);
  CHECK(saved.str() == nwk);
}

TEST_CASE("branch lengths are parsed correctly")
{
  tree_sptr_t tree = load_tree("((A:1,B:2)N1:3,C:4)root;");
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t b = find_leaf(tree, "B");
  node_sptr_t n1 = a->get_parent();
  REQUIRE(a != nullptr);
  REQUIRE(b != nullptr);
  CHECK(a->get_blen() == doctest::Approx(1.0));
  CHECK(b->get_blen() == doctest::Approx(2.0));
  CHECK(n1->get_blen() == doctest::Approx(3.0));
  CHECK(std::isnan(tree->get_root()->get_blen())); // labelled root without a length
  CHECK(a->get_midpoint_pendant() == doctest::Approx(0.5));
  CHECK(n1->get_tblen() > 0);
}

TEST_CASE("node depths describe the tree")
{
  tree_sptr_t tree = load_tree("((A:1,B:2)N1:3,C:4)root;");
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t b = find_leaf(tree, "B");
  node_sptr_t c = find_leaf(tree, "C");
  node_sptr_t n1 = a->get_parent();
  CHECK(tree->get_root()->get_ldepth() == 0);
  CHECK(n1->get_ldepth() == 1);
  CHECK(a->get_ldepth() == 2);
  CHECK(b->get_ldepth() == 2);
  CHECK(c->get_ldepth() == 1);
  // bdepth is the sum of the branch lengths from the root.
  CHECK(tree->get_root()->get_bdepth() == doctest::Approx(0.0));
  CHECK(n1->get_bdepth() == doctest::Approx(3.0));
  CHECK(a->get_bdepth() == doctest::Approx(4.0)); // root -> N1 (3) -> A (1)
  CHECK(b->get_bdepth() == doctest::Approx(5.0)); // root -> N1 (3) -> B (2)
  CHECK(c->get_bdepth() == doctest::Approx(4.0)); // root -> C (4)

  // A labelled root has no branch length of its own; a missing length counts as
  // zero rather than poisoning the depths below it with NaN.
  tree_sptr_t plain = load_tree("((A:1,B:1):2,C:3);");
  node_sptr_t inner = find_leaf(plain, "A")->get_parent();
  CHECK(std::isnan(plain->get_root()->get_blen()));
  CHECK(plain->get_root()->get_bdepth() == doctest::Approx(0.0));
  CHECK(inner->get_bdepth() == doctest::Approx(2.0));
  CHECK(find_leaf(plain, "A")->get_bdepth() == doctest::Approx(3.0));
}

TEST_CASE("generate_tree also fills in the depths")
{
  tree_sptr_t tree = make_tree({"A", "B", "C", "D"});
  CHECK(tree->get_root()->get_ldepth() == 0);
  for (const node_sptr_t& nd : post_order(tree)) {
    if (nd == tree->get_root()) continue;
    CHECK(nd->get_ldepth() == nd->get_parent()->get_ldepth() + 1);
    CHECK(nd->get_bdepth() == doctest::Approx(nd->get_parent()->get_bdepth() + nd->get_blen()));
  }
}

TEST_CASE("compute_lca and compute_distance walk a known tree")
{
  tree_sptr_t tree = load_tree("((A:1,B:2)N1:3,C:4)root;");
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t b = find_leaf(tree, "B");
  node_sptr_t c = find_leaf(tree, "C");
  node_sptr_t n1 = a->get_parent();
  // Reflexive and null cases.
  CHECK(Tree::compute_lca(a, a) == a);
  CHECK(Tree::compute_lca(nullptr, a) == a);
  CHECK(Tree::compute_lca(a, nullptr) == a);
  CHECK(Tree::compute_distance(a, a) == doctest::Approx(0.0));
  CHECK(Tree::compute_distance(nullptr, a) == std::numeric_limits<double>::max());
  // Siblings meet at their parent.
  CHECK(Tree::compute_lca(a, b) == n1);
  CHECK(Tree::compute_distance(a, b) == doctest::Approx(3.0));
  // A leaf and its aunt meet at the root.
  CHECK(Tree::compute_lca(a, c) == tree->get_root());
  CHECK(Tree::compute_distance(a, c) == doctest::Approx(8.0));
}

TEST_CASE("decorated trees report their edge numbers")
{
  tree_sptr_t tree = load_tree("((A:1{1},B:2{2})N1:3{3},C:4{4})root{0};");
  CHECK(tree->get_root()->get_en() == 0);
  CHECK(find_leaf(tree, "A")->get_en() == 1);
  CHECK(find_leaf(tree, "B")->get_en() == 2);
  CHECK(find_leaf(tree, "A")->get_parent()->get_en() == 3);
  // The jplace rendering keeps the decorations.
  const std::string jplace = print_jplace(tree);
  CHECK(jplace.find("{1}") != std::string::npos);
  CHECK(jplace.find("{0}") != std::string::npos);
  CHECK(jplace.back() == ';');

  // Without decorations the edge number is the post-order index minus one.
  tree_sptr_t plain = load_tree("((A:1,B:2)N1:3,C:4)root;");
  CHECK(plain->get_root()->get_en() == plain->get_nnodes() - 1);
  CHECK(plain->get_root()->check_decorated() == false);
}

TEST_CASE("partially decorated and duplicate-edge trees are rejected")
{
  ThrowingErrorHandler handler;
  {
    const std::string msg = ThrowingErrorHandler::catches([] { load_tree("((A:1{1},B:2)N1:3,C:4)root;"); });
    CHECK(msg.find("partially decorated") != std::string::npos);
  }
  {
    const std::string msg = ThrowingErrorHandler::catches([] { load_tree("((A:1{1},B:2{1})N1:3{3},C:4{4})root{0};"); });
    CHECK(msg.find("Duplicate decorated edge number") != std::string::npos);
  }
  {
    const std::string msg = ThrowingErrorHandler::catches([] { load_tree("((A:1{},B:2{2})N1:3{3},C:4{4})root{0};"); });
    CHECK(msg.find("must contain an edge number") != std::string::npos);
  }
  {
    const std::string msg = ThrowingErrorHandler::catches([] { load_tree("((A:1{1,B:2{2})N1:3{3},C:4{4})root{0};"); });
    CHECK(msg.find("missing the closing brace") != std::string::npos);
  }
}

TEST_CASE("malformed Newick input is reported")
{
  ThrowingErrorHandler handler;
  CHECK(ThrowingErrorHandler::catches([] { load_tree(""); }).find("empty") != std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A:1,B:1)"); }).find("ends with a character other than ';'") !=
        std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A:1,B:1);\n(C:1,D:1);"); }).find("multiple trees") !=
        std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A:1,B:1);(C:1,D:1);"); }).find("';'") != std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A:1,B:1)[comment;"); }).find("unquoted label") !=
        std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A B:1,C:1);"); }).find("' ' or newline") != std::string::npos);
}

TEST_CASE("duplicate labels and unifurcations are rejected")
{
  ThrowingErrorHandler handler;
  CHECK(ThrowingErrorHandler::catches([] { load_tree("(A:1,A:1);"); }).find("Duplicate node name") != std::string::npos);
  CHECK(ThrowingErrorHandler::catches([] { load_tree("((A:1)N1:1,B:1);"); }).find("single child") != std::string::npos);
}

TEST_CASE("check_compatible compares label sequences")
{
  tree_sptr_t a = load_tree("((A:1,B:1)N1:1,C:1)root;");
  tree_sptr_t b = load_tree("((A:1,B:1)N1:1,C:1)root;");
  tree_sptr_t c = load_tree("((A:1,C:1)N1:1,B:1)root;");
  CHECK(a->check_compatible(b));
  CHECK(a->check_compatible(nullptr));
  CHECK_FALSE(a->check_compatible(c));
}

TEST_CASE("subtree traversal is confined to the requested subtree")
{
  tree_sptr_t tree = load_tree("((A:1,B:2)N1:3,C:4)root;");
  node_sptr_t a = find_leaf(tree, "A");
  node_sptr_t b = find_leaf(tree, "B");
  node_sptr_t n1 = a->get_parent();

  auto walk = [&](node_sptr_t from) {
    std::vector<std::string> visited;
    tree->reset_traversal();
    tree->set_subtree(from);
    for (tuint_t i = 0; i < tree->get_nnodes(); ++i) {
      node_sptr_t nd = tree->next_post_order();
      if (!nd) break;
      visited.push_back(nd->get_name());
    }
    return visited;
  };

  CHECK(walk(a) == std::vector<std::string>{"A"});
  CHECK(walk(b) == std::vector<std::string>{"B"});
  CHECK(walk(n1) == std::vector<std::string>{"A", "B", "N1"});
  CHECK(walk(tree->get_root()) == std::vector<std::string>{"A", "B", "N1", "C", "root"});

  tree->reset_traversal();
  CHECK(tree->get_subtree() == tree->get_root());
  CHECK(post_order(tree).size() == tree->get_nnodes());
}

TEST_CASE("map_to_qtree re-points the index tree at a query tree")
{
  const std::vector<std::string> names = {"A", "B", "C", "D"};
  tree_sptr_t index_tree = make_tree(names);
  // The query tree has a different shape but the same leaves.
  tree_sptr_t qtree = load_tree("((A:1,C:2)X:1,(B:1,D:1)Y:1)Q;");
  index_tree->map_to_qtree(qtree);
  CHECK(index_tree->get_root() == qtree->get_root());
  for (const std::string& name : names) {
    node_sptr_t leaf = find_leaf(index_tree, name);
    REQUIRE(leaf != nullptr);
    CHECK(leaf->get_tree() == qtree);
    CHECK(leaf->get_name() == name);
  }
  // Only the leaf slots are remapped; interior slots keep the index tree's own
  // nodes (nothing looks those up after a remap).
  for (se_t se = 1; se <= index_tree->get_nnodes(); ++se) {
    node_sptr_t nd = index_tree->get_node(se);
    if (nd && nd->check_leaf()) CHECK(nd->get_tree() == qtree);
  }
  // eff_nchildren counts the covered children of every node.
  CHECK(qtree->get_root()->get_eff_nchildren() == 2);
}

TEST_CASE("map_to_qtree tolerates leaves that are absent from the query tree")
{
  tree_sptr_t index_tree = make_tree({"A", "B", "C", "D"});
  tree_sptr_t qtree = load_tree("((A:1,B:1)X:1,C:1)Q;");
  index_tree->map_to_qtree(qtree);
  CHECK(index_tree->get_root() == qtree->get_root());
  // D has no counterpart, so its slot stays empty.
  bool any_null = false;
  for (se_t se = 1; se <= index_tree->get_nnodes(); ++se) {
    if (!index_tree->get_node(se)) any_null = true;
  }
  CHECK(any_null);
}

TEST_CASE("parse_lineages builds a taxonomic tree")
{
  TempDir dir("lineages");
  const std::filesystem::path path = dir / "lineages.txt";
  spit(path,
       "A\tk__Bacteria; p__P1; g__G1\n"
       "B\tk__Bacteria; p__P1; g__G2\n"
       "C\tk__Bacteria; p__P2; g__G3\n");
  auto tree = std::make_shared<Tree>();
  std::ifstream in(path);
  tree->parse_lineages(in);
  tree->reset_traversal();
  CHECK(tree->get_root()->get_name() == "root");
  CHECK(tree->get_root()->check_taxon());

  std::vector<std::string> leaves;
  for (const node_sptr_t& nd : post_order(tree)) {
    if (nd->check_leaf()) leaves.push_back(nd->get_name());
  }
  std::sort(leaves.begin(), leaves.end());
  CHECK(leaves == std::vector<std::string>{"A", "B", "C"});

  node_sptr_t g1 = nullptr;
  for (const node_sptr_t& nd : post_order(tree)) {
    if (nd->get_name() == "G1") g1 = nd;
  }
  REQUIRE(g1 != nullptr);
  CHECK(g1->check_taxon());
  CHECK(g1->get_parent()->get_name() == "P1");
  CHECK(g1->get_parent()->get_parent()->get_name() == "Bacteria");
}

TEST_CASE("parse_lineages rejects duplicates and colliding names")
{
  ThrowingErrorHandler handler;
  TempDir dir("lineages-bad");
  {
    const std::filesystem::path path = dir / "dup.txt";
    spit(path, "A\tk__Bacteria; g__G1\nA\tk__Bacteria; g__G2\n");
    auto tree = std::make_shared<Tree>();
    std::ifstream in(path);
    CHECK(ThrowingErrorHandler::catches([&] { tree->parse_lineages(in); }).find("Duplicate reference ID") !=
          std::string::npos);
  }
  {
    const std::filesystem::path path = dir / "collide.txt";
    spit(path, "G1\tk__Bacteria; g__G1\n");
    auto tree = std::make_shared<Tree>();
    std::ifstream in(path);
    CHECK(ThrowingErrorHandler::catches([&] { tree->parse_lineages(in); }).find("collides with a taxon name") !=
          std::string::npos);
  }
  {
    const std::filesystem::path path = dir / "malformed.txt";
    spit(path, "no-tab-separated-lineage\n");
    auto tree = std::make_shared<Tree>();
    std::ifstream in(path);
    CHECK(ThrowingErrorHandler::catches([&] { tree->parse_lineages(in); }).find("Failed to reference to lineage") !=
          std::string::npos);
  }
}

TEST_CASE("print_info reports the tree size")
{
  tree_sptr_t tree = make_tree({"A", "B", "C"});
  CaptureStream quiet(std::cout);
  tree->print_info();
  CHECK(quiet.str().find("Number of nodes") != std::string::npos);
  CHECK(quiet.str().find("Total branch length") != std::string::npos);
}

TEST_CASE("node helpers agree with the tree shape")
{
  tree_sptr_t tree = load_tree("((A:1,B:1)N1:2,C:3)root;");
  node_sptr_t root = tree->get_root();
  node_sptr_t n1 = find_leaf(tree, "A")->get_parent();
  CHECK(root->get_nchildren() == 2);
  CHECK(n1->get_nchildren() == 2);
  CHECK(root->check_leaf() == false);
  CHECK(find_leaf(tree, "A")->check_leaf());
  CHECK(root->check_taxon() == false);
  CHECK(root->get_name(true) == "root");
  CHECK(n1->get_name(true) == "N1");
  // An unnamed node reports its index unless asked for a placeholder.
  tree_sptr_t unnamed = load_tree("((A:1,B:1):2,C:3);");
  node_sptr_t inner = find_leaf(unnamed, "A")->get_parent();
  CHECK(inner->get_name(true) == "NA");
  CHECK(inner->get_name() == std::to_string(inner->get_se() - 1));
  CHECK(inner->get_midpoint_pendant() == doctest::Approx(1.0));
  CHECK(root->sum_children_sh() != 0);
}

TEST_SUITE_END();
