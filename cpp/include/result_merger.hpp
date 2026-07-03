/**
 * @file result_merger.hpp
 * @brief Final assembly of the off-target table: sort each scored shard, then
 *        k-way merge the sorted shards into one globally ordered table.
 *
 * Pipeline position
 * -----------------
 * After the per-partition shards are written (search_executor) and scored
 * in place (Python per-shard scorer), this layer produces the final table.
 * It is memory-bounded: each shard is sorted in memory one at a time and
 * written back, then the sorted shards are merged by streaming a single row
 * from each at a time through a heap — peak memory is one shard for the sort
 * phase and one row per shard for the merge phase, never the whole result set.
 *
 * Sort modes
 * ----------
 *   - EditDistance (default): total edit distance ascending, then mismatches
 *     ascending, then bulges ascending, then CFD score descending (rows whose
 *     score is the "NA" sentinel sort last within a tie group).
 *   - Coordinates: contig lexicographic ascending, then position ascending.
 *
 * Both modes use the same comparator object in the per-shard sort and in the
 * merge heap, so the streamed output is globally ordered.
 *
 * Row schema (shared with the C++ ScoredTsvFormatter and the Python scorer):
 *   0 chrom  1 pos  2 strand  3 grna  4 spacer
 *   5 mismatches  6 bulge_type  7 bulge_dna  8 bulge_rna  9 cfd_score
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace crispritz {

// =============================================================================
// SortMode
// =============================================================================

/**
 * @brief Ordering applied to the final off-target table.
 *
 * API-only (selected via the search_offtargets_tst Python entry point, not a
 * CLI flag). Underlying values are stable.
 */
enum class SortMode : uint8_t {
  EditDistance = 0, ///< edits ↑, mm ↑, bulges ↑, CFD ↓ (NA last). Default.
  Coordinates = 1,  ///< contig (lexicographic) ↑, position ↑.
};

/**
 * @brief Lowercase canonical name of a SortMode
 * ("edit_distance"/"coordinates").
 */
[[nodiscard]] std::string_view to_string(SortMode mode) noexcept;

/**
 * @brief Parse a sort-mode token into the enum.
 * @param name "edit_distance" or "coordinates".
 * @throws std::invalid_argument on an unrecognised token.
 */
[[nodiscard]] SortMode sort_mode_from_string(std::string_view name);

// =============================================================================
// merge_sorted_shards
// =============================================================================

/**
 * @brief Sort each shard by @p mode, then k-way merge into @p final_path.
 *
 * Each shard is sorted in place (header preserved), then the sorted shards are
 * merged by a streaming heap into @p final_path with a single header line.
 *
 * @param shard_paths  Scored shard files (each: header + data rows).
 * @param final_path   Destination for the merged, globally-sorted table.
 * @param mode         Ordering to apply.
 * @param write_header Emit one header line at the top of @p final_path.
 * @param remove_inputs Delete each shard file after a successful merge.
 * @param verbosity   Output verbosity (0=Silent,1=Normal,2=Verbose,3=Debug).
 * @return Number of data rows written to @p final_path.
 *
 * @throws std::runtime_error if a shard or the output cannot be opened/written.
 */
std::size_t merge_sorted_shards(const std::vector<std::string> &shard_paths,
                                const std::string &final_path, SortMode mode,
                                bool write_header = true,
                                bool remove_inputs = true, int verbosity = 1);

} // namespace crispritz