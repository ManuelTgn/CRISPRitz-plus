/**
 * @file test_search_executor.cpp
 * @brief Unit tests for the search-executor output contract.
 *
 * Covers the deterministic, filesystem-free part: ScoredTsvFormatter, which
 * defines the 10-column shard schema shared with the Python per-shard scorer
 * (scores/shard_scoring.py). The end-to-end run_search_executor() test
 * (partition -> shard rows + per-guide profiles) needs the .bin builder
 * fixture currently private to test_tst_search.cpp; see the note at the bottom.
 *
 * Uses the same minimal record()/g_* harness as the other CRISPRitz C++ tests
 * so it registers with add_crispritz_test without an external framework.
 */

#include "offtarget.hpp"
#include "search_executor.hpp"

#include <iostream>
#include <string>
#include <vector>

using crispritz::OffTarget;
using crispritz::ScoredTsvFormatter;
using crispritz::Strand;

// =============================================================================
// Minimal harness
// =============================================================================

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

static void record(const std::string &name, bool ok,
                   const std::string &detail = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << name << '\n';
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << name;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << '\n';
  }
}

static std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> parts;
  std::string part;
  for (char c : s) {
    if (c == delim) {
      parts.push_back(part);
      part.clear();
    } else {
      part += c;
    }
  }
  parts.push_back(part);
  return parts;
}

// =============================================================================
// Header
// =============================================================================

static void test_header_columns() {
  ScoredTsvFormatter f;
  const std::string hdr = f.header();
  auto cols = split(hdr, '\t');
  const std::vector<std::string> expected = {
      "chrom",      "pos",        "strand",    "grna",      "spacer",
      "mismatches", "bulge_type", "bulge_dna", "bulge_rna", "cfd_score"};
  record("header has 10 columns", cols.size() == 10u,
         "got " + std::to_string(cols.size()));
  record("header column names/order match the scorer contract",
         cols == expected);
}

// =============================================================================
// Row layout
// =============================================================================

static void test_row_field_count() {
  ScoredTsvFormatter f;
  OffTarget ot{"chr1", 101, Strand::Forward, "ACGTANGG", "ACGTAaGG", 1, 0, 0};
  auto fields = split(f.format_row(ot), '\t');
  record("row has 10 tab-separated fields", fields.size() == 10u,
         "got " + std::to_string(fields.size()));
}

static void test_row_field_values_no_bulge() {
  ScoredTsvFormatter f;
  // 1 mismatch (lowercase 'a' in spacer), no bulges -> bulge_type "X".
  OffTarget ot{"chr1", 101, Strand::Forward, "ACGTANGG", "ACGTAaGG", 1, 0, 0};
  auto v = split(f.format_row(ot), '\t');
  record("field 0 chrom", v[0] == "chr1");
  record("field 1 pos", v[1] == "101");
  record("field 2 strand '+'", v[2] == "+");
  record("field 3 grna", v[3] == "ACGTANGG");
  record("field 4 spacer == OffTarget.target()", v[4] == "ACGTAaGG");
  record("field 5 mismatches", v[5] == "1");
  record("field 6 bulge_type 'X'", v[6] == "X");
  record("field 7 bulge_dna", v[7] == "0");
  record("field 8 bulge_rna", v[8] == "0");
  record("field 9 cfd_score is the NA sentinel", v[9] == "NA");
}

static void test_row_reverse_strand() {
  ScoredTsvFormatter f;
  OffTarget ot{"chrM", 500, Strand::Reverse, "ACGTANGG", "ACGTAGGG", 0, 0, 0};
  auto v = split(f.format_row(ot), '\t');
  record("reverse strand serialises as '-'", v[2] == "-");
}

static void test_row_dna_bulge_type() {
  ScoredTsvFormatter f;
  // DNA bulge: gap in grna, bulge_dna=1 -> bulge_type "DNA".
  OffTarget ot{"chr2", 50, Strand::Forward, "AC-GTANGG", "ACTGTAGGG", 0, 1, 0};
  auto v = split(f.format_row(ot), '\t');
  record("DNA-bulge row: bulge_type 'DNA'", v[6] == "DNA");
  record("DNA-bulge row: bulge_dna 1", v[7] == "1");
  record("DNA-bulge row: bulge_rna 0", v[8] == "0");
  record("DNA-bulge row: cfd_score NA", v[9] == "NA");
}

static void test_row_rna_bulge_type() {
  ScoredTsvFormatter f;
  // RNA bulge: gap in target, bulge_rna=1 -> bulge_type "RNA".
  OffTarget ot{"chr3", 9, Strand::Forward, "ACGGTANGG", "AC-GTAGGG", 0, 0, 1};
  auto v = split(f.format_row(ot), '\t');
  record("RNA-bulge row: bulge_type 'RNA'", v[6] == "RNA");
  record("RNA-bulge row: bulge_rna 1", v[8] == "1");
}

static void test_row_no_trailing_newline() {
  ScoredTsvFormatter f;
  OffTarget ot{"chr1", 1, Strand::Forward, "ACGTANGG", "ACGTAGGG", 0, 0, 0};
  const std::string row = f.format_row(ot);
  record("row has no trailing newline", !row.empty() && row.back() != '\n');
}

static void test_header_row_width_match() {
  ScoredTsvFormatter f;
  OffTarget ot{"chr1", 1, Strand::Forward, "ACGTANGG", "ACGTAGGG", 0, 0, 0};
  record("header column count == row field count",
         split(f.header(), '\t').size() ==
             split(f.format_row(ot), '\t').size());
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_search_executor ===\n\n";

  std::cout << "-- ScoredTsvFormatter header --\n";
  test_header_columns();

  std::cout << "\n-- ScoredTsvFormatter rows --\n";
  test_row_field_count();
  test_row_field_values_no_bulge();
  test_row_reverse_strand();
  test_row_dna_bulge_type();
  test_row_rna_bulge_type();
  test_row_no_trailing_newline();
  test_header_row_width_match();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  // NOTE: the end-to-end run_search_executor() test (build a tiny .bin, run the
  // executor, assert shard rows + per-guide profiles + targets/profile gating)
  // needs the fixture::make_bin / Builder helper currently private to
  // test_tst_search.cpp. The clean step is to promote that helper into a shared
  // cpp/tests/fixtures.hpp and include it here; that fixture refactor is the
  // companion to this test, not something to duplicate.

  return g_failed == 0 ? 0 : 1;
}