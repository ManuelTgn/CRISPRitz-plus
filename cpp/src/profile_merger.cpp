/**
 * @file profile_merge.cpp
 * @brief Implementation of cross-partition profile merging.
 */

#include "profile_merger.hpp"
#include "verbosity.hpp"

#include "output_writer.hpp" // ProfileWriter

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace crispritz {

namespace {

void add1(std::vector<int> &acc, const std::vector<int> &other,
          const char *field) {
  if (acc.size() != other.size())
    throw std::invalid_argument(
        std::string("merge_guide_profiles: size mismatch in ") + field);
  for (std::size_t i = 0; i < acc.size(); ++i)
    acc[i] += other[i];
}

void add2(std::vector<std::vector<int>> &acc,
          const std::vector<std::vector<int>> &other, const char *field) {
  if (acc.size() != other.size())
    throw std::invalid_argument(
        std::string("merge_guide_profiles: outer size mismatch in ") + field);
  for (std::size_t i = 0; i < acc.size(); ++i)
    add1(acc[i], other[i], field);
}

void add3(std::vector<std::vector<std::vector<int>>> &acc,
          const std::vector<std::vector<std::vector<int>>> &other,
          const char *field) {
  if (acc.size() != other.size())
    throw std::invalid_argument(
        std::string("merge_guide_profiles: outer size mismatch in ") + field);
  for (std::size_t i = 0; i < acc.size(); ++i)
    add2(acc[i], other[i], field);
}

/** @brief Verify two profiles describe the same guide and search geometry. */
void check_geometry(const GuideProfile &a, const GuideProfile &b) {
  if (a.guide != b.guide || a.guide_len != b.guide_len ||
      a.pam_len != b.pam_len || a.pam_at_start != b.pam_at_start ||
      a.max_mm != b.max_mm || a.max_bulge_dna != b.max_bulge_dna ||
      a.max_bulge_rna != b.max_bulge_rna)
    throw std::invalid_argument(
        "merge_guide_profiles: profiles for guide '" + a.guide +
        "' have inconsistent geometry across partitions");
}

/** @brief Accumulate every counter field of @p other into @p acc. */
void add_into(GuideProfile &acc, const GuideProfile &other) {
  check_geometry(acc, other);

  // File 1
  add1(acc.pos_mm_count, other.pos_mm_count, "pos_mm_count");
  acc.ont_count += other.ont_count;
  add1(acc.offt_by_mm, other.offt_by_mm, "offt_by_mm");

  // File 2
  add3(acc.ext_mm_nuc_pos, other.ext_mm_nuc_pos, "ext_mm_nuc_pos");
  add2(acc.ext_total_by_mm, other.ext_total_by_mm, "ext_total_by_mm");
  add2(acc.ext_dna_by_mm_pos, other.ext_dna_by_mm_pos, "ext_dna_by_mm_pos");
  add2(acc.ext_rna_by_mm_pos, other.ext_rna_by_mm_pos, "ext_rna_by_mm_pos");

  // File 3
  add1(acc.pos_bulge_dna, other.pos_bulge_dna, "pos_bulge_dna");
  add1(acc.pos_mm_in_dna, other.pos_mm_in_dna, "pos_mm_in_dna");
  add2(acc.offt_dna, other.offt_dna, "offt_dna");
  acc.ont_count_dna += other.ont_count_dna;

  // File 4
  add1(acc.pos_bulge_rna, other.pos_bulge_rna, "pos_bulge_rna");
  add1(acc.pos_mm_in_rna, other.pos_mm_in_rna, "pos_mm_in_rna");
  add2(acc.offt_rna, other.offt_rna, "offt_rna");
  acc.ont_count_rna += other.ont_count_rna;

  // File 5
  add1(acc.offt_complete_by_mm, other.offt_complete_by_mm,
       "offt_complete_by_mm");
  acc.ont_count_complete += other.ont_count_complete;
  add1(acc.pos_mm_complete, other.pos_mm_complete, "pos_mm_complete");
}

} // namespace

std::vector<GuideProfile> merge_guide_profiles(
    const std::vector<std::vector<GuideProfile>> &by_partition) {
  if (by_partition.empty())
    return {};

  // Seed with the first partition's profiles (geometry + counters).
  std::vector<GuideProfile> merged = by_partition.front();
  const std::size_t guide_count = merged.size();

  for (std::size_t p = 1; p < by_partition.size(); ++p) {
    const std::vector<GuideProfile> &part = by_partition[p];
    if (part.size() != guide_count)
      throw std::invalid_argument(
          "merge_guide_profiles: partition " + std::to_string(p) + " has " +
          std::to_string(part.size()) + " guides, expected " +
          std::to_string(guide_count));
    for (std::size_t g = 0; g < guide_count; ++g)
      add_into(merged[g], part[g]);
  }
  return merged;
}

std::size_t write_merged_profiles(
    const std::vector<std::vector<GuideProfile>> &by_partition,
    const std::string &path_stem, int verbosity) {
  const std::vector<GuideProfile> merged = merge_guide_profiles(by_partition);
  ProfileWriter writer;
  writer.write_all_profiles(merged, path_stem);
  // Single-threaded finalization step; level-2 output is safe here.
  print_verbosity("Wrote profiles for " + std::to_string(merged.size()) +
                      " guide(s) to " + path_stem + ".profile*.xls",
                  verbosity, VERBOSITY_VERBOSE);
  return merged.size();
}

} // namespace crispritz