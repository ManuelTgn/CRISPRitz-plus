/**
 * @file profile_merge.hpp
 * @brief Combine per-partition profiles into per-guide totals and write them.
 *
 * The search executor produces one GuideProfile per guide *per partition*
 * (ProfileAccumulator has no cross-partition merge). This layer sums those
 * profiles element-wise per guide and hands the result to ProfileWriter.
 *
 * Input layout: @c by_partition[p][g] is the profile for guide @p g from
 * partition @p p. Every partition must carry the same guides in the same
 * order and with identical geometry (guide length, edit budget, PAM), which
 * holds by construction (all partitions are searched with the same guides and
 * SearchConfiguration).
 */

#pragma once

#include "profile_data.hpp" // GuideProfile

#include <string>
#include <vector>

namespace crispritz {

/**
 * @brief Sum per-partition profiles into one GuideProfile per guide.
 *
 * @param by_partition  Outer index = partition, inner index = guide.
 * @return One merged GuideProfile per guide (empty if @p by_partition is
 *         empty). Geometry fields are taken from the first partition; all
 *         counter fields are summed.
 *
 * @throws std::invalid_argument if partitions disagree on guide count or
 *         per-guide geometry (a sign the inputs were not produced from the
 *         same search).
 */
[[nodiscard]] std::vector<GuideProfile> merge_guide_profiles(
    const std::vector<std::vector<GuideProfile>> &by_partition);

/**
 * @brief Merge per-partition profiles and write the five profile files.
 *
 * Convenience wrapper: @c merge_guide_profiles() followed by
 * @c ProfileWriter::write_all_profiles(merged, path_stem).
 *
 * @param by_partition  Per-partition, per-guide profiles.
 * @param path_stem     Shared path prefix for the five @c .xls files.
 * @return Number of guides written.
 *
 * @throws std::invalid_argument on inconsistent inputs (see above).
 * @throws std::runtime_error if any profile file cannot be written.
 */
std::size_t write_merged_profiles(
    const std::vector<std::vector<GuideProfile>> &by_partition,
    const std::string &path_stem);

} // namespace crispritz