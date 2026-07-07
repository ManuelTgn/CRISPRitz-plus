/**
 * @file bindings.cpp
 * @brief pybind11 surface for the CRISPRitz-plus C++ core.
 *
 * The boundary is deliberately narrow and matches the streaming pipeline: the
 * heavy data (off-target rows) never crosses into Python — it is written to
 * per-partition shard files by the C++ executor, scored in place by Python,
 * then merged back by C++. Only configuration, small profile structs, and file
 * paths cross.
 *
 *   Python owns: partition-level parallelism (ThreadPoolExecutor over
 *   run_search_executor), per-shard CFD scoring (ProcessPoolExecutor), and the
 *   run loop.
 *   C++ owns: load + search + streaming shard writes (run_search_executor),
 *   per-shard sort + k-way merge (merge_sorted_shards), and profile
 *   accumulation + writing (write_merged_profiles).
 *
 * Every long-running entry releases the GIL via py::call_guard so the Python
 * thread pool achieves true C++ parallelism.
 */

#include "nucleotide_encoding.hpp"
#include "profile_data.hpp"         // GuideProfile
#include "profile_merger.hpp"       // write_merged_profiles
#include "result_merger.hpp"        // SortMode, merge_sorted_shards
#include "search_configuration.hpp" // SearchConfiguration, OutputFormat/Mode
#include "search_executor.hpp"      // run_search_executor, PartitionResult
#include "tst.hpp"                  // build_tree
#include "tst_search.hpp"           // BulgeMode, bulge_mode_from_string

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // std::string / std::vector auto-conversion

#include <stdexcept>
#include <string>

namespace py = pybind11;
using namespace crispritz;

PYBIND11_MODULE(_ternary_search_tree, m) {
  m.doc() = "CRISPRitz-plus C++ API bindings (pybind11): TST index "
            "construction and the streaming per-partition search pipeline.";

  // C++ runtime/argument errors -> descriptive Python exceptions.
  py::register_exception<std::runtime_error>(m, "TSTBuildError");
  py::register_exception<std::invalid_argument>(m, "TSTSearchError");

  // =========================================================================
  // TST index construction
  // =========================================================================
  m.def("build_tree", &crispritz::build_tree, py::arg("sequence"),
        py::arg("chr_name"), py::arg("pam_seq"), py::arg("pam_length"),
        py::arg("pam_limit"), py::arg("upstream"), py::arg("outdir"),
        py::arg("max_bulges") = 0, py::arg("num_threads") = 1,
        py::arg("verbosity") = 1,
        "Build a Ternary Search Tree index for one genomic sequence and write "
        "the .bin partition(s) to outdir.");

  // =========================================================================
  // Enums + string helpers
  // =========================================================================
  py::enum_<OutputFormat>(m, "OutputFormat")
      .value("Tsv", OutputFormat::Tsv)
      .value("Targets", OutputFormat::Targets);

  py::enum_<OutputMode>(m, "OutputMode")
      .value("TargetsOnly", OutputMode::TargetsOnly)
      .value("ProfileOnly", OutputMode::ProfileOnly)
      .value("Both", OutputMode::Both);

  py::enum_<SortMode>(m, "SortMode", "Ordering of the final off-target table.")
      .value("EditDistance", SortMode::EditDistance)
      .value("Coordinates", SortMode::Coordinates);

  py::enum_<BulgeMode>(m, "BulgeMode",
                       "Whether one alignment may mix DNA and RNA bulges.")
      .value("MixedBulges", BulgeMode::MixedBulges)
      .value("SingleBulgeType", BulgeMode::SingleBulgeType);

  m.def(
      "output_format_from_string",
      [](const std::string &n) { return output_format_from_string(n); },
      py::arg("name"));
  m.def(
      "output_mode_from_string",
      [](const std::string &n) { return output_mode_from_string(n); },
      py::arg("name"));
  m.def(
      "sort_mode_from_string",
      [](const std::string &n) { return sort_mode_from_string(n); },
      py::arg("name"));
  m.def(
      "bulge_mode_from_string",
      [](const std::string &s) { return bulge_mode_from_string(s); },
      py::arg("s"));

  // =========================================================================
  // SearchConfiguration
  // =========================================================================
  py::class_<SearchConfiguration>(m, "SearchConfiguration")
      .def(py::init<int, int, int, int, OutputFormat, OutputMode>(),
           py::arg("max_mismatches"), py::arg("max_bulges_dna"),
           py::arg("max_bulges_rna"), py::arg("threads"),
           py::arg("output_format") = OutputFormat::Tsv,
           py::arg("output_mode") = OutputMode::Both)
      .def_property_readonly("max_mismatches",
                             &SearchConfiguration::max_mismatches)
      .def_property_readonly("max_bulges_dna",
                             &SearchConfiguration::max_bulges_dna)
      .def_property_readonly("max_bulges_rna",
                             &SearchConfiguration::max_bulges_rna)
      .def_property_readonly("threads", &SearchConfiguration::threads)
      .def_property_readonly("output_format",
                             &SearchConfiguration::output_format)
      .def_property_readonly("output_mode", &SearchConfiguration::output_mode)
      .def_property_readonly("max_total_edits",
                             &SearchConfiguration::max_total_edits)
      .def_property_readonly("write_targets",
                             &SearchConfiguration::write_targets)
      .def_property_readonly("write_profile",
                             &SearchConfiguration::write_profile);

  // =========================================================================
  // GuideProfile (opaque handle)
  // =========================================================================
  py::class_<GuideProfile>(m, "GuideProfile")
      .def_readonly("guide", &GuideProfile::guide)
      .def_readonly("guide_len", &GuideProfile::guide_len)
      .def_readonly("ont_count", &GuideProfile::ont_count)
      .def_readonly("ont_count_complete", &GuideProfile::ont_count_complete)
      .def("__repr__", [](const GuideProfile &p) {
        return "<GuideProfile guide='" + p.guide +
               "' ont=" + std::to_string(p.ont_count_complete) + ">";
      });

  // =========================================================================
  // PartitionResult
  // =========================================================================
  py::class_<PartitionResult>(m, "PartitionResult")
      .def_readonly("source_path", &PartitionResult::source_path)
      .def_readonly("shard_path", &PartitionResult::shard_path)
      .def_readonly("total_hits", &PartitionResult::total_hits)
      .def_readonly("rows_written", &PartitionResult::rows_written)
      .def_readonly("profiles", &PartitionResult::profiles)
      .def("__repr__", [](const PartitionResult &r) {
        return "<PartitionResult src='" + r.source_path +
               "' hits=" + std::to_string(r.total_hits) + ">";
      });

  // =========================================================================
  // run_search_executor — one partition: load + search + stream to shard
  // =========================================================================
  m.def("run_search_executor", &crispritz::run_search_executor,
        py::arg("partition_path"), py::arg("chrom"), py::arg("guides"),
        py::arg("config"), py::arg("pam"), py::arg("pam_at_start"),
        py::arg("shard_path"), py::arg("bulge_mode") = BulgeMode::MixedBulges,
        py::arg("verbosity") = 1, py::call_guard<py::gil_scoped_release>(),
        "Load one .bin partition, search every guide, and stream the hits to a "
        "shard file (targets) and per-guide profiles. Returns a "
        "PartitionResult.");

  // =========================================================================
  // merge_sorted_shards — sort each scored shard, k-way merge to final table
  // =========================================================================
  m.def(
      "merge_sorted_shards", &crispritz::merge_sorted_shards,
      py::arg("shard_paths"), py::arg("final_path"), py::arg("sort_mode"),
      py::arg("write_header") = true, py::arg("remove_inputs") = true,
      py::arg("verbosity") = 1, py::call_guard<py::gil_scoped_release>(),
      "Sort each scored shard by sort_mode, then k-way merge into final_path. "
      "Returns the number of rows written.");

  // =========================================================================
  // write_merged_profiles — merge per-partition profiles and write the files
  // =========================================================================
  m.def("write_merged_profiles", &crispritz::write_merged_profiles,
        py::arg("profiles_by_partition"), py::arg("path_stem"),
        py::arg("verbosity") = 1, py::call_guard<py::gil_scoped_release>(),
        "Sum per-partition, per-guide profiles and write the five .xls profile "
        "files at path_stem. Returns the number of guides written.");
}