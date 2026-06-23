/**
 * @file test_tst_search.cpp
 * @brief End-to-end search tests over a REAL on-disk index.
 *
 * test_tst_search.cpp validates the traversal algorithm in isolation, but it
 * builds its .bin files with a hand-rolled re-implementation of the writer
 * (fixture::Builder + make_bin). That deliberately keeps the algorithm test
 * self-contained, but it means the test writer and the production writer are
 * two separate code paths: a serializer bug (nibble phase drift, sentinel/null
 * alignment, header geometry) can break production while that test stays green.
 *
 * This file closes that gap. It drives the production pipeline end to end:
 *
 *     build_tree()      (real TernarySearchTree builder + serializer)
 *        -> .bin file    (real on-disk partition)
 *        -> load_partition()   (real deserializer)
 *        -> TSTSearcher / search_partition()   (real traversal)
 *
 * so the writer and reader are tested against each other on a real file, which
 * is where production index bugs actually live.
 *
 * ## Determinism strategy
 * Hit counts are only predictable if the genome's PAM-site set is known, so the
 * fixture genome is engineered to contain EXACTLY ONE PAM site: a single "GG"
 * (one forward NGG site), no "CC" (no reverse-strand site), and no "N". That
 * yields exactly one indexed leaf.
 *
 * ## Orientation independence
 * For a PAM-at-end index the builder reverses the guide before storage. To keep
 * these tests independent of that convention, the fixture guide is a
 * reverse-palindrome ("ACGCA"), so the canonical stored form equals the query
 * regardless of whether the builder reverses. The mismatch query ("ACTCA") is
 * likewise a palindrome differing in exactly one position. Do not "simplify"
 * these to non-palindromic strings — the palindrome is load-bearing.
 *
 * ## Scope
 * Exact match, mismatch budget, absent guide, multi-guide grouping, provenance,
 * and the deferred edit-budget / guide-length validation. Bulge ALIGNMENT
 * semantics are covered by test_tst_search.cpp against controlled trees;
 * predicting exact bulge emissions end to end is fragile and out of scope here.
 *
 * Uses the shared lightweight record()/g_* harness.
 */

#include "tst_search.hpp"

#include "offtarget.hpp"
#include "search_configuration.hpp"
#include "tst.hpp" // build_tree

#include <chrono>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace crispritz;

// =============================================================================
// Harness
// =============================================================================

static int g_total = 0, g_passed = 0, g_failed = 0;

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

// =============================================================================
// Fixture
// =============================================================================

namespace {

// PAM-at-end "NGG" index over a 5-bp guide.
const std::string kPamSeq = "NGG";
constexpr int kPamLimit = 3;
constexpr int kGuideLen = 5;
constexpr int kPamLength = kGuideLen + kPamLimit; // 8
constexpr bool kPamAtEnd = false;
const std::string kChr = "chrE2E";

// Exactly one PAM site: one "GG", no "CC", no "N". The core "ACGCATGG" places
// the palindromic guide "ACGCA" immediately 5' of the PAM "TGG"; the "ATATAT"
// padding introduces no further GG/CC.
const std::string kGenome = "ATATATACGCATGGATATAT";

// Reverse-palindrome guide => canonical stored form == query, regardless of
// the builder's PAM-at-end guide reversal.
const std::string kGuide = "ACGCA";       // present, exact
const std::string kGuide1mm = "ACTCA";    // differs from kGuide at index 2
const std::string kGuideAbsent = "TTTTT"; // not in the index

std::string bin_name() { return kPamSeq + "_" + kChr + "_1.bin"; }

/** @brief RAII unique temp directory (created now, removed with contents). */
struct TempDir {
  fs::path path;
  bool ok = false;
  TempDir() {
    std::error_code ec;
    static int counter = 0;
    const auto stamp =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    path = fs::temp_directory_path(ec) /
           ("crispritz_search_e2e_" + std::to_string(stamp) + "_" +
            std::to_string(counter++));
    if (ec)
      return;
    ok = fs::create_directories(path, ec) && !ec;
  }
  ~TempDir() {
    std::error_code ec;
    if (!path.empty())
      fs::remove_all(path, ec);
  }
  std::string str() const { return path.string(); }
};

/** @brief Build the single-site index into @p tmp; return the .bin path. */
std::string build_single_site_index(const TempDir &tmp) {
  build_tree(kGenome, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
             tmp.str(), /*max_bulges=*/0, /*num_threads=*/1);
  return (tmp.path / bin_name()).string();
}

/** @brief A config with the given mismatch budget and no bulges. */
SearchConfiguration cfg(int mm) {
  return SearchConfiguration(mm, /*bdna=*/0, /*brna=*/0, /*threads=*/1);
}

} // namespace

// =============================================================================
// Round-trip into a searchable index
// =============================================================================

/**
 * @brief The real builder + loader produce a single-leaf, correctly-shaped
 *        index that the searcher can consume.
 */
static void test_real_index_roundtrips() {
  TempDir tmp;
  if (!tmp.ok) {
    record("real index round-trips", false, "temp dir setup failed");
    return;
  }
  std::string path;
  try {
    path = build_single_site_index(tmp);
  } catch (const std::exception &e) {
    record("real index round-trips", false,
           std::string("build threw: ") + e.what());
    return;
  }

  if (!fs::exists(path)) {
    record("real index round-trips", false, "no .bin written: " + path);
    return;
  }

  try {
    LoadedTST tst = load_partition(path);
    record("loaded index has exactly one leaf", tst.leaf_count() == 1u,
           "got " + std::to_string(tst.leaf_count()));
    record("loaded guide_length == 5", tst.guide_length() == kGuideLen);
    record("loaded pam_limit == 3", tst.pam_limit() == kPamLimit);
    record("loaded node pool non-empty", tst.node_count() >= 1u);
  } catch (const std::exception &e) {
    record("real index round-trips", false,
           std::string("load_partition threw: ") + e.what());
  }
}

// =============================================================================
// Exact match
// =============================================================================

/**
 * @brief The indexed guide is found with a zero-edit budget, and the hit is a
 *        clean forward-strand 0-mismatch, 0-bulge alignment.
 */
static void test_exact_match_found() {
  TempDir tmp;
  if (!tmp.ok) {
    record("exact match found", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);
  TSTSearcher searcher(cfg(0));

  std::vector<OffTarget> hits =
      searcher.search(tst, kGuide, kChr, kPamSeq, kPamAtEnd);

  record("exact match: exactly one hit", hits.size() == 1u,
         "got " + std::to_string(hits.size()));
  if (hits.size() != 1u)
    return;

  const OffTarget &h = hits.front();
  record("exact match: 0 mismatches", h.mismatches() == 0);
  record("exact match: no bulge", !h.has_bulge());
  record("exact match: forward strand", h.strand() == Strand::Forward);
  record("exact match: chrom propagated", h.chrom() == kChr);
  record("exact match: 1-based position is positive", h.pos() > 0);
}

// =============================================================================
// Mismatch budget
// =============================================================================

/**
 * @brief A 1-mismatch query is found iff the budget allows one mismatch.
 */
static void test_one_mismatch_budget() {
  TempDir tmp;
  if (!tmp.ok) {
    record("one-mismatch budget", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);

  // Budget 0: the 1-mismatch query must not match.
  {
    TSTSearcher searcher(cfg(0));
    std::vector<OffTarget> hits =
        searcher.search(tst, kGuide1mm, kChr, kPamSeq, kPamAtEnd);
    record("1-mm query with budget 0 yields no hit", hits.empty(),
           "got " + std::to_string(hits.size()));
  }

  // Budget 1: the 1-mismatch query matches with exactly one mismatch.
  {
    TSTSearcher searcher(cfg(1));
    std::vector<OffTarget> hits =
        searcher.search(tst, kGuide1mm, kChr, kPamSeq, kPamAtEnd);
    record("1-mm query with budget 1 finds the site", hits.size() == 1u,
           "got " + std::to_string(hits.size()));
    if (hits.size() == 1u)
      record("1-mm query: hit reports exactly one mismatch",
             hits.front().mismatches() == 1,
             "got " + std::to_string(hits.front().mismatches()));
  }

  // Budget 1: the exact guide is still found (0 mismatches), unaffected by the
  // larger budget.
  {
    TSTSearcher searcher(cfg(1));
    std::vector<OffTarget> hits =
        searcher.search(tst, kGuide, kChr, kPamSeq, kPamAtEnd);
    record("exact guide still found under a larger budget",
           hits.size() == 1u && hits.front().mismatches() == 0);
  }
}

// =============================================================================
// Absent guide
// =============================================================================

/**
 * @brief A guide that is not in the index (and not within budget of anything)
 *        returns no hits.
 */
static void test_absent_guide() {
  TempDir tmp;
  if (!tmp.ok) {
    record("absent guide yields no hit", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);
  TSTSearcher searcher(cfg(0));

  std::vector<OffTarget> hits =
      searcher.search(tst, kGuideAbsent, kChr, kPamSeq, kPamAtEnd);
  record("absent guide yields no hit", hits.empty(),
         "got " + std::to_string(hits.size()));
}

// =============================================================================
// Multi-guide grouping + provenance (search_all / search_partition)
// =============================================================================

/**
 * @brief search_all groups results in lockstep with the input guide order and
 *        records the source partition.
 */
static void test_search_all_grouping() {
  TempDir tmp;
  if (!tmp.ok) {
    record("search_all grouping", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);
  TSTSearcher searcher(cfg(0));

  const std::vector<std::string> guides = {kGuide, kGuideAbsent};
  SearchResult r = searcher.search_all(tst, guides, kChr, kPamSeq, kPamAtEnd);

  record("search_all: one result list per guide", r.guide_count() == 2u);
  record("search_all: present guide (index 0) has a hit",
         r.guide_count() == 2u && r.hits_by_guide[0].size() == 1u);
  record("search_all: absent guide (index 1) has no hits",
         r.guide_count() == 2u && r.hits_by_guide[1].empty());
  record("search_all: total_hits counts across guides", r.total_hits() == 1u);
  record("search_all: source_path is the loaded partition",
         r.source_path == path);
}

/**
 * @brief search_partition() loads and searches in one call, the primitive the
 *        Python layer drives per .bin file.
 */
static void test_search_partition_one_call() {
  TempDir tmp;
  if (!tmp.ok) {
    record("search_partition one-call", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);

  const std::vector<std::string> guides = {kGuide, kGuideAbsent};
  SearchResult r =
      search_partition(path, kChr, guides, cfg(0), kPamSeq, kPamAtEnd);

  record("search_partition: one list per guide", r.guide_count() == 2u);
  record("search_partition: present guide has a hit",
         r.guide_count() == 2u && !r.hits_by_guide[0].empty());
  record("search_partition: absent guide has none",
         r.guide_count() == 2u && r.hits_by_guide[1].empty());
  record("search_partition: source_path preserved", r.source_path == path);
}

// =============================================================================
// Deferred validation (edit budget vs index guide length)
// =============================================================================

/**
 * @brief An edit budget exceeding the index guide length is rejected at search
 *        time (the check the standalone SearchConfiguration cannot perform).
 */
static void test_budget_exceeds_guide_length_rejected() {
  TempDir tmp;
  if (!tmp.ok) {
    record("over-budget search rejected", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);

  // guide_length is 5; a 6-mismatch budget cannot fit.
  TSTSearcher searcher(cfg(kGuideLen + 1));
  bool threw = false;
  try {
    searcher.search(tst, kGuide, kChr, kPamSeq, kPamAtEnd);
  } catch (const std::invalid_argument &) {
    threw = true;
  } catch (...) {
    record("over-budget search throws invalid_argument", false,
           "wrong exception type");
    return;
  }
  record("over-budget search throws invalid_argument", threw);
}

/**
 * @brief A query whose length disagrees with the index guide length is
 *        rejected.
 */
static void test_query_length_mismatch_rejected() {
  TempDir tmp;
  if (!tmp.ok) {
    record("query-length mismatch rejected", false, "temp dir setup failed");
    return;
  }
  std::string path = build_single_site_index(tmp);
  LoadedTST tst = load_partition(path);
  TSTSearcher searcher(cfg(0));

  bool threw = false;
  try {
    searcher.search(tst, "ACG", kChr, kPamSeq, kPamAtEnd); // length 3 != 5
  } catch (const std::invalid_argument &) {
    threw = true;
  } catch (...) {
    record("query-length mismatch throws invalid_argument", false,
           "wrong exception type");
    return;
  }
  record("query-length mismatch throws invalid_argument", threw);
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_tst_search_e2e ===\n\n";

  std::cout << "-- build -> load round-trip --\n";
  test_real_index_roundtrips();

  std::cout << "\n-- exact match --\n";
  test_exact_match_found();

  std::cout << "\n-- mismatch budget --\n";
  test_one_mismatch_budget();

  std::cout << "\n-- absent guide --\n";
  test_absent_guide();

  std::cout << "\n-- multi-guide grouping / provenance --\n";
  test_search_all_grouping();
  test_search_partition_one_call();

  std::cout << "\n-- deferred validation --\n";
  test_budget_exceeds_guide_length_rejected();
  test_query_length_mismatch_rejected();

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}