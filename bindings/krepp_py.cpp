// pybind11 bindings for krepp's embedding API (src/api.hpp).
// Results come back as plain Python objects (lists of row objects); no
// streams, no temp files. Heavy calls release the GIL while OpenMP works.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "api.hpp"

namespace py = pybind11;
using namespace api;

namespace {
  // Decode C++ strings as UTF-8 with surrogateescape (like Python does for
  // filenames): identifiers from hostile FASTA/FASTQ files may contain raw
  // bytes, and a strict decode would throw UnicodeDecodeError and discard the
  // whole result batch. Surrogates round-trip and print safely.
  py::str pystr(const std::string& s)
  {
    return py::reinterpret_steal<py::str>(PyUnicode_DecodeUTF8(s.data(), static_cast<Py_ssize_t>(s.size()), "surrogateescape"));
  }
} // namespace

PYBIND11_MODULE(_krepp, m)
{
  m.doc() = "krepp: k-mer-based search, distance estimation and phylogenetic placement";
  m.attr("__version__") = VERSION;

  py::class_<Options>(m, "Options", "Query options; field names and defaults match the CLI flags of `krepp dist`/`krepp place`.")
    .def(py::init<>())
    .def_readwrite("hdist_th", &Options::hdist_th, "Maximum Hamming distance for a k-mer to match (--hdist-th, [0,31]). [4]")
    .def_readwrite("chisq_value", &Options::chisq_value, "Chi-square value for the distinguishability test (--chisq). [2.706]")
    .def_readwrite("dist_max", &Options::dist_max, "Maximum distance to report (--dist-max, [1e-8, 0.33]). [NaN]")
    .def_readwrite("tau", &Options::tau, "Highest Hamming distance for the placement threshold (--tau). [2]")
    .def_readwrite("multi", &Options::multi, "Report all candidates, not just the best (--multi). [True]")
    .def_readwrite("p_value", &Options::p_value, "dist: compute the p-value column (--p-value). [False]")
    .def_readwrite("dist_filter", &Options::dist_filter, "dist: filter by statistical significance (dist --filter). [False]")
    .def_readwrite("place_filter", &Options::place_filter, "place: filter by k-mer support (place --filter). [True]")
    .def_readwrite("num_threads", &Options::num_threads, "OpenMP threads (--num-threads); >1 makes result order nondeterministic, as in the CLI. [1]")
    .def_readwrite("verbose", &Options::verbose, "Per-phase load timing on stderr (--verbose). [False]")
    .def_readwrite("use_mmap", &Options::use_mmap, "View index arrays through a file mapping (--mmap/--no-mmap). [True]")
    .def_readwrite("lineage_path", &Options::lineage_path, "place: Greengenes/GTDB lineage file (-l).")
    .def_readwrite("nwk_path", &Options::nwk_path, "place: Newick backbone tree (-t); overrides the index's tree.")
    .def("__repr__", [](const Options& o) {
      return "<krepp.Options hdist_th=" + std::to_string(o.hdist_th) + " tau=" + std::to_string(o.tau) +
             " multi=" + (o.multi ? "True" : "False") + " num_threads=" + std::to_string(o.num_threads) + ">";
    });

  py::class_<dist_row_t>(m, "DistRow", "One `krepp dist` output row.")
    .def_property_readonly("ref_name", [](const dist_row_t& r) { return pystr(r.ref_name); }, "REFERENCE_NAME")
    .def_readonly("dist", &dist_row_t::d_llh, "DIST (maximum-likelihood k-mer distance)")
    .def_readonly("p_value", &dist_row_t::p_value, "P_VALUE; NaN unless Options.p_value is set")
    .def("__repr__", [](const dist_row_t& r) {
      return "<DistRow " + r.ref_name + " dist=" + std::to_string(r.d_llh) + ">";
    });

  py::class_<query_dist_t>(m, "QueryDist", "All distance rows for one query sequence.")
    .def_property_readonly("query", [](const query_dist_t& q) { return pystr(q.query); }, "SEQ_ID")
    .def_readonly("na", &query_dist_t::na, "True when the CLI would print the NA/NaN row")
    .def_readonly("rows", &query_dist_t::rows, "list[DistRow], sorted like CLI output");

  py::class_<placement_t>(m, "Placement", "One placement; fields are the jplace \"p\" array plus the distal node name.")
    .def_readonly("edge_num", &placement_t::edge_num)
    .def_readonly("pendant_length", &placement_t::pendant_length)
    .def_readonly("distal_length", &placement_t::distal_length)
    .def_readonly("likelihood", &placement_t::likelihood)
    .def_readonly("like_weight_ratio", &placement_t::like_weight_ratio, "LWR")
    .def_readonly("distance", &placement_t::distance)
    .def_property_readonly("distal_node", [](const placement_t& p) { return pystr(p.distal_node); })
    .def("__repr__", [](const placement_t& p) {
      return "<Placement edge=" + std::to_string(p.edge_num) + " lwr=" + std::to_string(p.like_weight_ratio) + ">";
    });

  py::class_<query_place_t>(m, "QueryPlace", "All placements for one query sequence.")
    .def_property_readonly("query", [](const query_place_t& q) { return pystr(q.query); })
    .def_readonly("placed", &query_place_t::placed, "False when every candidate was filtered out (query absent from jplace output)")
    .def_readonly("placements", &query_place_t::placements, "list[Placement]");

  py::class_<match_row_t>(m, "MatchRow", "Per-reference k-mer match summary for one query.")
    .def_property_readonly("ref_name", [](const match_row_t& r) { return pystr(r.ref_name); })
    .def_readonly("match_count", &match_row_t::match_count, "query k-mers with at least one match within hdist_th")
    .def_readonly("mismatch_count", &match_row_t::mismatch_count, "onmers - match_count")
    .def_readonly("nmers", &match_row_t::nmers)
    .def_readonly("rho", &match_row_t::rho, "index subsampling rate for this reference")
    .def_readonly("hdist_min", &match_row_t::hdist_min, "smallest observed Hamming distance")
    .def_readonly("hdisthist", &match_row_t::hdisthist, "histogram; index = Hamming distance, value = k-mer count")
    .def_readonly("dist", &match_row_t::d_llh, "ML distance (same value dist() reports for this reference)");

  py::class_<query_match_t>(m, "QueryMatch", "k-mer match summary for one query sequence.")
    .def_property_readonly("query", [](const query_match_t& q) { return pystr(q.query); })
    .def_readonly("onmers", &query_match_t::onmers, "valid k-mers scanned in the query")
    .def_readonly("enmers", &query_match_t::enmers, "expected k-mers: max(0, len - k + 1)")
    .def_readonly("hdist_hist",
                  &query_match_t::hdist_hist,
                  "histogram of match counts by Hamming distance, summed over all matched "
                  "references (a k-mer matching several references counts once per reference); "
                  "index = Hamming distance")
    .def_readonly("rows", &query_match_t::rows, "list[MatchRow], sorted by (dist, ref_name)");

  py::class_<Session>(m, "Session",
                      "A loaded index. Methods mirror krepp's subcommands and return plain data.\n\n"
                      "Corrupt indexes/files raise RuntimeError; invalid paths and option\n"
                      "values raise ValueError. Session calls are serialised process-wide\n"
                      "(the underlying settings are globals); OpenMP parallelism within a\n"
                      "call is unaffected.")
    .def(py::init<const std::filesystem::path&, const Options&>(),
         py::arg("index_dir"),
         py::arg("options") = Options{},
         "Load the index directory (the -i of `krepp dist`/`krepp place`).")
    .def("dist",
         [](Session& s, const std::string& q) {
           py::gil_scoped_release nogil;
           return s.dist(q);
         },
         py::arg("fasta"),
         "Equivalent of `krepp dist`: per-query distance rows (list[QueryDist]).")
    .def("place",
         [](Session& s, const std::string& q) {
           py::gil_scoped_release nogil;
           return s.place(q);
         },
         py::arg("fasta"),
         "Equivalent of `krepp place`: per-query placements (list[QueryPlace]).\n"
         "Needs a tree: the index's backbone, Options.nwk_path, or Options.lineage_path.")
    .def("matches",
         [](Session& s, const std::string& q) {
           py::gil_scoped_release nogil;
           return s.matches(q);
         },
         py::arg("fasta"),
         "Per-query k-mer match summaries per matched reference (list[QueryMatch]).")
    .def_property_readonly("total_qseq", &Session::get_total_qseq, "Number of query sequences processed so far.")
    .def_property_readonly("options", &Session::get_options, py::return_value_policy::reference_internal);
}
