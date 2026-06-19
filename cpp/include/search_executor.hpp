/**
 * @file search_executor.hpp
 * @brief Per-partition execution unit: load a partition, search every guide,
 *        and stream the hits to a shard file and/or per-guide profile
 *        accumulators — without ever holding more than one guide's hits in
 *        memory at once.
 *
 * Position in the architecture
 * ----------------------------
 * Python owns the partition-level parallelism: it submits one
 * @c run_search_executor() call per @c .bin partition to a ThreadPoolExecutor,
 * and the pybind11 binding releases the GIL so the calls run as parallel C++.
 * Each call is independent (its own LoadedTST, its own shard file, its own
 * accumulators), so no locking is required — the shard-per-partition layout is
 * what makes the concurrent writes safe without a shared file or mutex.
 *
 * Streaming
 * ---------
 * The executor searches one guide at a time via @c TSTSearcher::search() and
 * sinks each guide's hits immediately into a threshold-flushed
 * @c OutputWriter::Session (targets) and that guide's @c ProfileAccumulator
 * (profile). A guide's hit vector is released before the next guide is
 * searched, so peak memory is bounded by the largest single guide's result
 * set plus the session buffer — not by the whole partition. (Bounding memory
 * below one guide's result set would require a per-hit sink inside the
 * traversal, a separate change to the search hot path.)
 *
 * Output contract
 * ---------------
 * Shard rows use the fixed 10-column scored-TSV schema expected by the Python
 * per-shard scorer (scores/shard_scoring.py), with the trailing @c cfd_score
 * column written as the @c "NA" sentinel for the scorer to fill in place:
 *
 *     chrom  pos  strand  grna  spacer  mismatches  bulge_type
 *     bulge_dna  bulge_rna  cfd_score
 *
 * This schema is independent of SearchConfiguration::output_format; the legacy
 * "targets" layout is a final-presentation transform applied later, at merge
 * time, not at the shard.
 */

#pragma once

#include "offtarget.hpp"            // OffTarget
#include "output_writer.hpp"        // OffTargetFormatter
#include "profile_data.hpp"         // GuideProfile
#include "search_configuration.hpp" // SearchConfiguration

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace crispritz {

// =========================================================================
// ScoredTsvFormatter
// =========================================================================

/**
 * @brief Formatter for the per-shard scored-TSV schema.
 *
 * Emits the 10-column layout shared with the Python scorer. The C++ domain
 * model calls the genomic field @c target(); it is written into the @c spacer
 * column to match the Python @c OffTarget naming. The @c cfd_score column is
 * always the @c "NA" sentinel here — it is filled in place by the Python
 * per-shard scoring pass after the search completes.
 */
class ScoredTsvFormatter final : public OffTargetFormatter {
public:
  /** @brief Sentinel written into the cfd_score column at search time. */
  static constexpr std::string_view SCORE_NA = "NA";

  [[nodiscard]] std::string header() const override;
  [[nodiscard]] std::string format_row(const OffTarget &ot) const override;
  [[nodiscard]] std::string_view name() const noexcept override {
    return "scored_tsv";
  }
};

// =========================================================================
// PartitionResult
// =========================================================================

/**
 * @brief Outcome of executing the search over a single partition.
 *
 * @c profiles holds one @c GuideProfile per input guide (indexed in lockstep
 * with the @c guides argument) when profiling is enabled, or is empty
 * otherwise. Because @c ProfileAccumulator has no cross-partition merge, these
 * per-partition profiles are combined per guide by a separate merge step
 * before @c ProfileWriter::write_all_profiles().
 */
struct PartitionResult {
  /// Originating @c .bin partition path (provenance/logging).
  std::string source_path;
  /// Shard file the targets were written to ("" when targets are disabled).
  std::string shard_path;
  /// Total hits found across all guides in this partition.
  std::size_t total_hits = 0;
  /// Rows actually written to the shard (== total_hits when targets enabled).
  std::size_t rows_written = 0;
  /// One profile per guide (empty when profiling is disabled).
  std::vector<GuideProfile> profiles;
};

// =========================================================================
// Flush threshold
// =========================================================================

/**
 * @brief Records buffered before the shard Session auto-flushes to disk.
 *
 * Set to 100,000 per the streaming-output design (overrides the writer's
 * 1,000,000 default), bounding the session buffer during genome-wide runs.
 */
inline constexpr std::size_t SHARD_FLUSH_THRESHOLD = 100'000;

// =========================================================================
// run_search_executor
// =========================================================================

/**
 * @brief Load one partition, search every guide, and stream the results.
 *
 * Loads @p partition_path, then for each guide searches it and sinks the hits
 * into the shard Session (when @c config.write_targets()) and that guide's
 * @c ProfileAccumulator (when @c config.write_profile()). The shard file is
 * written with a header and auto-flushed every @c SHARD_FLUSH_THRESHOLD rows.
 *
 * @param partition_path Path to a single @c .bin partition file.
 * @param chrom          Chromosome name recorded in every emitted OffTarget
 *                       (resolved Python-side from the partition filename).
 * @param guides         Query guides in canonical orientation.
 * @param config         Validated search parameters (gates targets/profile).
 * @param pam_len        Number of PAM bases (from the PAM model; the .bin
 *                       header does not carry it — see file header).
 * @param pam_at_start   True when the PAM precedes the guide body (Cas12a).
 * @param shard_path     Destination shard file for the targets table. Used
 *                       only when @c config.write_targets(); ignored otherwise.
 *
 * @return A PartitionResult with per-guide profiles and counts.
 *
 * @throws std::runtime_error    if the partition cannot be loaded or a shard
 *                               write fails.
 * @throws std::invalid_argument if the edit budget is incompatible with the
 *                               loaded guide length, or a guide length does
 *                               not match the index.
 */
[[nodiscard]] PartitionResult
run_search_executor(const std::string &partition_path, const std::string &chrom,
                    const std::vector<std::string> &guides,
                    const SearchConfiguration &config, int pam_len,
                    bool pam_at_start, const std::string &shard_path);

} // namespace crispritz