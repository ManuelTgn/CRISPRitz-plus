/**
 * @file test_profile_merge.cpp
 * @brief Unit tests for cross-partition profile merging
 *        (profile_merge.hpp / profile_merge.cpp).
 *
 * Verifies that every GuideProfile counter rank (scalar / 1D / 2D / 3D) is
 * summed across partitions, that geometry is preserved and validated, and that
 * write_merged_profiles routes the merged result to ProfileWriter (a .xls file
 * is produced).
 */

#include "profile_data.hpp"
#include "profile_merger.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace crispritz;

static int g_total = 0, g_passed = 0, g_failed = 0;
static void record(const std::string &n, bool ok) {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << n << '\n';
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << n << '\n';
  }
}

// Build a GuideProfile (guide_len=2, max_mm=1, bulge_dna=1, bulge_rna=1) with
// every counter set to a constant `v`, so sums are easy to predict.
static GuideProfile mk(const std::string &guide, int v) {
  GuideProfile p;
  p.guide = guide;
  p.guide_len = 2;
  p.pam_len = 3;
  p.pam_at_start = false;
  p.max_mm = 1;
  p.max_bulge_dna = 1;
  p.max_bulge_rna = 1;
  const int mm = p.max_mm + 1, L = p.guide_len;
  p.pos_mm_count.assign(L, v);
  p.ont_count = v;
  p.offt_by_mm.assign(mm, v);
  p.ext_mm_nuc_pos.assign(
      mm, std::vector<std::vector<int>>(4, std::vector<int>(L, v)));
  p.ext_total_by_mm.assign(mm, std::vector<int>(L, v));
  p.ext_dna_by_mm_pos.assign(mm, std::vector<int>(L, v));
  p.ext_rna_by_mm_pos.assign(mm, std::vector<int>(L, v));
  p.pos_bulge_dna.assign(L, v);
  p.pos_mm_in_dna.assign(L, v);
  p.offt_dna.assign(mm, std::vector<int>(p.max_bulge_dna, v));
  p.ont_count_dna = v;
  p.pos_bulge_rna.assign(L, v);
  p.pos_mm_in_rna.assign(L, v);
  p.offt_rna.assign(mm, std::vector<int>(p.max_bulge_rna, v));
  p.ont_count_rna = v;
  p.offt_complete_by_mm.assign(mm, v);
  p.ont_count_complete = v;
  p.pos_mm_complete.assign(L, v);
  return p;
}

int main() {
  std::cout << "=== test_profile_merge ===\n\n";

  // 3 partitions x 2 guides; partition contributions 1, 2, 4 -> every
  // counter 7.
  std::vector<std::vector<GuideProfile>> by_part = {{mk("AA", 1), mk("CC", 1)},
                                                    {mk("AA", 2), mk("CC", 2)},
                                                    {mk("AA", 4), mk("CC", 4)}};

  std::cout << "-- merge_guide_profiles --\n";
  auto merged = merge_guide_profiles(by_part);
  record("2 guides returned", merged.size() == 2u);
  const GuideProfile &g0 = merged[0];
  record("geometry preserved (guide, max_mm)",
         g0.guide == "AA" && g0.max_mm == 1);
  record("scalar ont_count summed", g0.ont_count == 7);
  record("scalar ont_count_dna summed", g0.ont_count_dna == 7);
  record("scalar ont_count_rna summed", g0.ont_count_rna == 7);
  record("scalar ont_count_complete summed", g0.ont_count_complete == 7);
  record("1D pos_mm_count summed",
         g0.pos_mm_count[0] == 7 && g0.pos_mm_count[1] == 7);
  record("1D offt_by_mm summed", g0.offt_by_mm[1] == 7);
  record("1D pos_bulge_dna summed", g0.pos_bulge_dna[0] == 7);
  record("1D pos_mm_in_rna summed", g0.pos_mm_in_rna[1] == 7);
  record("1D offt_complete_by_mm summed", g0.offt_complete_by_mm[1] == 7);
  record("1D pos_mm_complete summed", g0.pos_mm_complete[0] == 7);
  record("2D ext_total_by_mm summed", g0.ext_total_by_mm[1][1] == 7);
  record("2D ext_dna_by_mm_pos summed", g0.ext_dna_by_mm_pos[0][1] == 7);
  record("2D ext_rna_by_mm_pos summed", g0.ext_rna_by_mm_pos[1][0] == 7);
  record("2D offt_dna summed", g0.offt_dna[1][0] == 7);
  record("2D offt_rna summed", g0.offt_rna[0][0] == 7);
  record("3D ext_mm_nuc_pos summed", g0.ext_mm_nuc_pos[1][3][1] == 7);
  record("second guide summed",
         merged[1].guide == "CC" && merged[1].ont_count == 7);
  record("inputs not mutated", by_part[0][0].ont_count == 1);
  record("empty input -> empty output", merge_guide_profiles({}).empty());

  std::cout << "\n-- guards --\n";
  bool threw = false;
  try {
    (void)merge_guide_profiles({{mk("AA", 1)}, {mk("AG", 1)}});
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("geometry mismatch throws invalid_argument", threw);
  threw = false;
  try {
    (void)merge_guide_profiles({{mk("AA", 1), mk("CC", 1)}, {mk("AA", 1)}});
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("guide-count mismatch throws invalid_argument", threw);

  std::cout << "\n-- write_merged_profiles --\n";
  const std::string stem = "/tmp/cz_pm_stem";
  std::size_t n = write_merged_profiles(by_part, stem);
  record("returns guide count (2)", n == 2u);
  record(".profile.xls file produced",
         std::ifstream(stem + ".profile.xls").good());
  for (const char *suf :
       {".profile.xls", ".extended_profile.xls", ".profile_dna.xls",
        ".profile_rna.xls", ".profile_complete.xls"})
    std::remove((stem + suf).c_str());

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";
  return g_failed == 0 ? 0 : 1;
}