/**
 * @file test_tst.cpp
 * @brief Unit tests for the TST index builder (tst.hpp / tst.cpp).
 *
 * Unlike test_tst_search.cpp — which drives a *mocked* tree builder to exercise
 * the search path in isolation — this file tests the REAL
 * @c crispritz::TernarySearchTree end to end: it builds an index from a
 * synthetic genomic sequence, lets the production serializer write a real
 * @c .bin partition to a temp directory, then reads it back with the
 * production @c crispritz::load_partition. That round-trip is the point: it is
 * the only test that exercises the actual writer and reader against each other
 * on a real file, which is where the serialization bug class (nibble phase
 * drift, sentinel/null alignment, header format) actually lives.
 *
 * Scope:
 *   - Constructor argument validation (guide_length / pam_limit / outdir).
 *   - build() is a silent no-op when the sequence has no PAM site.
 *   - N-contaminated sites are never indexed (throw or zero leaves).
 *   - build() -> .bin file is written with the expected name.
 *   - The serialized header (magic / version / guide_length / pam_limit) is
 *     byte-correct.
 *   - build() -> load_partition() round-trips header geometry, leaf count,
 *     and per-leaf PAM byte width.
 *   - The free-function pybind entry point build_tree() behaves like the class.
 *
 * Assertions deliberately avoid hard-coding the number of PAM sites found:
 * that depends on the (both-strand) PAM search and is out of scope here. The
 * tests assert invariants that must hold for any valid index instead.
 *
 * Uses the same minimal record()/g_* harness as the other CRISPRitz C++ tests
 * so it registers with add_crispritz_test without an external framework.
 */

#include "tst.hpp"        // TernarySearchTree, TSTNode, TSTLeaf, build_tree
#include "tst_search.hpp" // LoadedTST, load_partition
#include "tst_utils.hpp"  // TST_BIN_MAGIC, TST_BIN_VERSION

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using crispritz::build_tree;
using crispritz::load_partition;
using crispritz::LoadedTST;
using crispritz::TernarySearchTree;

// =============================================================================
// Minimal test harness
// =============================================================================

static int g_total = 0;
static int g_passed = 0;
static int g_failed = 0;

static void record(const std::string &name, bool ok,
                   const std::string &detail = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << name << "\n";
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << name;
    if (!detail.empty())
      std::cout << " -- " << detail;
    std::cout << "\n";
  }
}

/**
 * @brief Record a pass iff @p fn throws std::runtime_error.
 *
 * The constructor contract (tst.hpp) is specifically std::runtime_error, so a
 * different exception type is a contract violation, not a pass.
 */
template <typename Fn>
static void expect_runtime_error(const std::string &name, Fn &&fn) {
  try {
    fn();
    record(name, false, "no exception thrown");
  } catch (const std::runtime_error &) {
    record(name, true);
  } catch (const std::exception &e) {
    record(name, false, std::string("wrong exception type: ") + e.what());
  } catch (...) {
    record(name, false, "non-std::exception thrown");
  }
}

// =============================================================================
// Test fixture parameters (project convention: 5-bp guide, 3-bp "NGG" PAM,
// PAM-at-end orientation).
// =============================================================================

namespace {
constexpr int kPamLimit = 3;                      // "NGG"
constexpr int kGuideLen = 5;                      // ACGTA
constexpr int kPamLength = kGuideLen + kPamLimit; // 8
const std::string kPamSeq = "NGG";
const std::string kChr = "testchr";
constexpr bool kPamAtEnd = false; // SpCas9-style: PAM follows the guide

/** @brief ceil(pam_limit / 2): packed PAM byte width per leaf. */
constexpr int kPamBytes = (kPamLimit + 1) / 2; // 2

/** @brief Expected partition filename for part 1: "<pam>_<chr>_1.bin". */
std::string expected_bin_name() { return kPamSeq + "_" + kChr + "_1.bin"; }

/**
 * @brief RAII unique temp directory. Created on construction, removed (with
 *        contents) on destruction so a failing test cannot leak files.
 */
struct TempDir {
  fs::path path;
  bool ok = false;

  TempDir() {
    std::error_code ec;
    static int counter = 0;
    const auto stamp =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    path = fs::temp_directory_path(ec) /
           ("crispritz_tst_test_" + std::to_string(stamp) + "_" +
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

/** @brief Count regular files with a ".bin" extension in @p dir. */
int count_bin_files(const fs::path &dir) {
  int n = 0;
  std::error_code ec;
  for (auto it = fs::directory_iterator(dir, ec);
       !ec && it != fs::directory_iterator(); it.increment(ec)) {
    if (it->is_regular_file(ec) && it->path().extension() == ".bin")
      ++n;
  }
  return n;
}
} // namespace

// =============================================================================
// Constructor validation
// =============================================================================

/**
 * @brief guide_length <= 0 (pam_length <= pam_limit) must be rejected.
 */
static void test_ctor_rejects_nonpositive_guide_length() {
  expect_runtime_error("ctor throws when guide_length == 0", [] {
    // pam_length == pam_limit -> guide_length 0
    TernarySearchTree(/*sequence=*/"ACGT", kChr, kPamSeq,
                      /*pam_length=*/kPamLimit, /*pam_limit=*/kPamLimit,
                      kPamAtEnd, /*outdir=*/".");
  });
  expect_runtime_error("ctor throws when guide_length < 0", [] {
    TernarySearchTree("ACGT", kChr, kPamSeq, /*pam_length=*/2,
                      /*pam_limit=*/5, kPamAtEnd, ".");
  });
}

/**
 * @brief pam_limit <= 0 must be rejected (guide_length check passes first).
 */
static void test_ctor_rejects_nonpositive_pam_limit() {
  expect_runtime_error("ctor throws when pam_limit == 0", [] {
    TernarySearchTree("ACGT", kChr, kPamSeq, /*pam_length=*/8,
                      /*pam_limit=*/0, kPamAtEnd, ".");
  });
}

/**
 * @brief An empty outdir must be rejected.
 */
static void test_ctor_rejects_empty_outdir() {
  expect_runtime_error("ctor throws when outdir is empty", [] {
    TernarySearchTree("ACGT", kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                      /*outdir=*/"");
  });
}

/**
 * @brief A well-formed configuration must construct without throwing. No build
 *        is run, so this only checks the constructor's validation path.
 */
static void test_ctor_accepts_valid_config() {
  TempDir tmp;
  if (!tmp.ok) {
    record("ctor accepts a valid configuration", false,
           "temp dir setup failed");
    return;
  }
  bool threw = false;
  try {
    TernarySearchTree tst("ACGTACGT", kChr, kPamSeq, kPamLength, kPamLimit,
                          kPamAtEnd, tmp.str());
    // Nothing built yet: counters must be zero.
    record("freshly constructed tree has 0 leaves", tst.leaf_count() == 0);
    record("freshly constructed tree has 0 nodes", tst.node_count() == 0);
  } catch (...) {
    threw = true;
  }
  record("ctor accepts a valid configuration", !threw);
}

// =============================================================================
// build(): empty / degenerate inputs
// =============================================================================

/**
 * @brief A sequence with no PAM occurrence (no "GG"/"CC" on either strand) is a
 *        silent no-op: no exception, no file, zero leaves and nodes.
 */
static void test_build_no_pam_site_is_noop() {
  TempDir tmp;
  if (!tmp.ok) {
    record("build() with no PAM site is a no-op", false,
           "temp dir setup failed");
    return;
  }

  // All-A: contains neither "GG" (forward NGG) nor "CC" (reverse-strand NGG).
  const std::string seq(30, 'A');

  bool threw = false;
  TernarySearchTree tst(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                        tmp.str());
  try {
    tst.build();
  } catch (...) {
    threw = true;
  }

  record("build() with no PAM site does not throw", !threw);
  record("build() with no PAM site produces 0 leaves", tst.leaf_count() == 0);
  record("build() with no PAM site writes no .bin file",
         count_bin_files(tmp.path) == 0);
}

/**
 * @brief Every candidate site is N-contaminated, so none can be indexed.
 *
 * The contract (build() doc) is to throw std::runtime_error when all sites are
 * discarded. Some PAM searches may instead pre-filter the N window and report
 * no site at all, making build() a no-op. Both outcomes satisfy the real
 * invariant under test: an N-contaminated window must never become an indexed
 * leaf. The test accepts either, but rejects "succeeded with leaves".
 */
static void test_build_rejects_n_contaminated_sites() {
  TempDir tmp;
  if (!tmp.ok) {
    record("N-contaminated sites are not indexed", false,
           "temp dir setup failed");
    return;
  }

  // Single "AGG" (matches NGG); its 5-bp guide window is all N. No "CC"
  // anywhere, so there is no reverse-strand site to rescue it.
  const std::string seq = "NNNNNAGGNNNNN";

  TernarySearchTree tst(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                        tmp.str());
  bool threw = false;
  try {
    tst.build();
  } catch (const std::runtime_error &) {
    threw = true;
  } catch (...) {
    record("N-contaminated sites are not indexed", false,
           "unexpected exception type");
    return;
  }

  const bool no_leaves = tst.leaf_count() == 0;
  record("N-contaminated sites are not indexed (throw or zero leaves)",
         threw || no_leaves);
  record("N-contaminated build writes no usable .bin",
         count_bin_files(tmp.path) == 0 || tst.leaf_count() == 0);
}

// =============================================================================
// build(): successful index + serialized header
// =============================================================================

/**
 * @brief A clean sequence with at least one PAM site produces exactly one
 *        partition file with the expected name, and a positive leaf/node count.
 */
static void test_build_writes_expected_partition() {
  TempDir tmp;
  if (!tmp.ok) {
    record("build() writes expected partition", false, "temp dir setup failed");
    return;
  }

  // Clean sequence containing PAM sites ("AGG", "CGG") and no N.
  const std::string seq = "TTTTACGTACGTAGGAAAAACGGTTTTACGTA";

  TernarySearchTree tst(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                        tmp.str());
  bool threw = false;
  try {
    tst.build();
  } catch (const std::exception &e) {
    threw = true;
    record("build() on a clean sequence does not throw", false, e.what());
  }
  if (threw)
    return;

  record("build() on a clean sequence does not throw", true);
  record("build() produced at least one leaf", tst.leaf_count() >= 1);
  record("build() allocated at least one node", tst.node_count() >= 1);

  const fs::path bin = tmp.path / expected_bin_name();
  record("expected partition file exists", fs::exists(bin),
         "looked for " + bin.string());
  record("exactly one .bin partition was written",
         count_bin_files(tmp.path) == 1);
}

/**
 * @brief The serialized header is byte-correct: magic, version, and the
 *        guide_length / pam_limit geometry that travels with the index.
 */
static void test_serialized_header_is_correct() {
  TempDir tmp;
  if (!tmp.ok) {
    record("serialized header is byte-correct", false, "temp dir setup failed");
    return;
  }

  const std::string seq = "TTTTACGTACGTAGGAAAAACGGTTTTACGTA";
  TernarySearchTree tst(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                        tmp.str());
  try {
    tst.build();
  } catch (const std::exception &e) {
    record("serialized header is byte-correct", false,
           std::string("build threw: ") + e.what());
    return;
  }

  const fs::path bin = tmp.path / expected_bin_name();
  std::ifstream in(bin, std::ios::in | std::ios::binary);
  if (!in.is_open()) {
    record("serialized header is byte-correct", false,
           "cannot open " + bin.string());
    return;
  }

  std::uint32_t magic = 0, version = 0;
  std::int32_t chunk_size = 0, guide_length = 0, pam_limit = 0;
  in.read(reinterpret_cast<char *>(&magic), sizeof(magic));
  in.read(reinterpret_cast<char *>(&version), sizeof(version));
  in.read(reinterpret_cast<char *>(&chunk_size), sizeof(chunk_size));
  in.read(reinterpret_cast<char *>(&guide_length), sizeof(guide_length));
  in.read(reinterpret_cast<char *>(&pam_limit), sizeof(pam_limit));
  const bool read_ok = static_cast<bool>(in);

  record("header read succeeded", read_ok);
  record("header magic == TST_BIN_MAGIC", magic == crispritz::TST_BIN_MAGIC,
         "got " + std::to_string(magic));
  record("header version == TST_BIN_VERSION",
         version == crispritz::TST_BIN_VERSION);
  record("header guide_length == 5", guide_length == kGuideLen);
  record("header pam_limit == 3", pam_limit == kPamLimit);
  record("header chunk_size matches builder leaf_count",
         chunk_size == tst.leaf_count());
}

// =============================================================================
// build() -> load_partition() round-trip
// =============================================================================

/**
 * @brief The production reader reconstructs a consistent index from what the
 *        production writer emitted. This is the core serializer<->deserializer
 *        contract test on a real .bin file.
 */
static void test_build_load_roundtrip() {
  TempDir tmp;
  if (!tmp.ok) {
    record("build/load round-trip", false, "temp dir setup failed");
    return;
  }

  const std::string seq = "TTTTACGTACGTAGGAAAAACGGTTTTACGTA";
  TernarySearchTree tst(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd,
                        tmp.str());
  try {
    tst.build();
  } catch (const std::exception &e) {
    record("build/load round-trip", false,
           std::string("build threw: ") + e.what());
    return;
  }

  const fs::path bin = tmp.path / expected_bin_name();
  if (!fs::exists(bin)) {
    record("build/load round-trip", false,
           "partition not written: " + bin.string());
    return;
  }

  LoadedTST loaded = [&] {
    try {
      return load_partition(bin.string());
    } catch (const std::exception &e) {
      record("load_partition succeeds on a freshly built index", false,
             e.what());
      throw;
    }
  }();
  record("load_partition succeeds on a freshly built index", true);

  // Geometry travels with the index, not the search config.
  record("loaded guide_length == 5", loaded.guide_length() == kGuideLen);
  record("loaded pam_limit == 3", loaded.pam_limit() == kPamLimit);

  // Single small input -> single partition: leaf counts must agree exactly.
  record("loaded leaf_count == builder leaf_count",
         static_cast<int>(loaded.leaf_count()) == tst.leaf_count(),
         "loaded=" + std::to_string(loaded.leaf_count()) +
             " builder=" + std::to_string(tst.leaf_count()));

  // The reader must produce a non-empty node pool (root + structure).
  record("loaded node pool is non-empty", loaded.node_count() >= 1);

  // Source path is retained for diagnostics.
  record("loaded source_path matches input",
         loaded.source_path() == bin.string());

  // Per-leaf invariants that survive serialization (guide_seq is NOT stored in
  // the leaf section — it lives in the node path — so it is intentionally not
  // checked here; PAM bytes and a non-zero strand-encoding index are).
  bool pam_width_ok = true;
  bool index_nonzero = true;
  for (const auto &leaf : loaded.leaves()) {
    if (static_cast<int>(leaf.pam_seq_enc.size()) != kPamBytes)
      pam_width_ok = false;
    if (leaf.guide_index == 0)
      index_nonzero = false;
  }
  record("every loaded leaf has ceil(pam_limit/2) PAM bytes", pam_width_ok);
  record("every loaded leaf has a non-zero (strand-encoded) index",
         index_nonzero);
}

/**
 * @brief The pybind entry-point free function build_tree() must produce the
 *        same on-disk artifact as the class it wraps.
 */
static void test_free_function_build_tree() {
  TempDir tmp;
  if (!tmp.ok) {
    record("build_tree() free function writes a loadable index", false,
           "temp dir setup failed");
    return;
  }

  const std::string seq = "TTTTACGTACGTAGGAAAAACGGTTTTACGTA";
  try {
    build_tree(seq, kChr, kPamSeq, kPamLength, kPamLimit, kPamAtEnd, tmp.str(),
               /*max_bulges=*/0, /*num_threads=*/1);
  } catch (const std::exception &e) {
    record("build_tree() free function writes a loadable index", false,
           std::string("threw: ") + e.what());
    return;
  }

  const fs::path bin = tmp.path / expected_bin_name();
  if (!fs::exists(bin)) {
    record("build_tree() free function writes a loadable index", false,
           "no partition written");
    return;
  }

  try {
    LoadedTST loaded = load_partition(bin.string());
    record("build_tree() free function writes a loadable index",
           loaded.guide_length() == kGuideLen &&
               loaded.pam_limit() == kPamLimit && loaded.leaf_count() >= 1);
  } catch (const std::exception &e) {
    record("build_tree() free function writes a loadable index", false,
           e.what());
  }
}

// =============================================================================
// main
// =============================================================================

int main() {
  std::cout << "=== test_tst ===\n\n";

  std::cout << "-- constructor validation --\n";
  test_ctor_rejects_nonpositive_guide_length();
  test_ctor_rejects_nonpositive_pam_limit();
  test_ctor_rejects_empty_outdir();
  test_ctor_accepts_valid_config();

  std::cout << "\n-- build(): degenerate inputs --\n";
  test_build_no_pam_site_is_noop();
  test_build_rejects_n_contaminated_sites();

  std::cout << "\n-- build(): partition + header --\n";
  test_build_writes_expected_partition();
  test_serialized_header_is_correct();

  std::cout << "\n-- build() -> load_partition() round-trip --\n";
  test_build_load_roundtrip();
  test_free_function_build_tree();

  std::cout << "\n=== Results: " << g_passed << "/" << g_total << " passed";
  if (g_failed > 0)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";

  return g_failed == 0 ? 0 : 1;
}