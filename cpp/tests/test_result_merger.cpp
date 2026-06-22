/**
 * @file test_result_merger.cpp
 * @brief Unit tests for SortMode + per-shard sort and k-way merge
 *        (result_merger.hpp / result_merger.cpp).
 *
 * Self-contained: builds small shard files under /tmp, merges them, and checks
 * the output against an independent oracle — the merged table must be a
 * permutation of the inputs AND globally ordered by the mode's comparator.
 */

#include "result_merger.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace crispritz;

static int g_total = 0, g_passed = 0, g_failed = 0;
static void record(const std::string &n, bool ok, const std::string &d = "") {
  ++g_total;
  if (ok) {
    ++g_passed;
    std::cout << "  [PASS] " << n << '\n';
  } else {
    ++g_failed;
    std::cout << "  [FAIL] " << n;
    if (!d.empty())
      std::cout << " -- " << d;
    std::cout << '\n';
  }
}

static const char *HDR =
    "chrom\tpos\tstrand\tgrna\tspacer\tmismatches\tbulge_type\t"
    "bulge_dna\tbulge_rna\tcfd_score";

static std::string row(const std::string &chrom, long pos, char strand, int mm,
                       int bd, int br, const std::string &cfd) {
  const std::string bt = (bd == 0 && br == 0)  ? "X"
                         : (bd > 0 && br == 0) ? "DNA"
                         : (bd == 0 && br > 0) ? "RNA"
                                               : "DNA,RNA";
  std::ostringstream o;
  o << chrom << '\t' << pos << '\t' << strand << "\tGRNA\tSPACER\t" << mm
    << '\t' << bt << '\t' << bd << '\t' << br << '\t' << cfd;
  return o.str();
}

static void write_shard(const std::string &p,
                        const std::vector<std::string> &rows) {
  std::ofstream o(p);
  o << HDR << '\n';
  for (const auto &r : rows)
    o << r << '\n';
}

static std::vector<std::string> read_data(const std::string &p) {
  std::ifstream in(p);
  std::string line;
  std::vector<std::string> out;
  bool first = true;
  while (std::getline(in, line)) {
    if (first) {
      first = false;
      continue;
    }
    if (!line.empty())
      out.push_back(line);
  }
  return out;
}

struct K {
  int te, mm, bg;
  double sc;
  std::string chrom;
  long pos;
};
static K key(const std::string &line) {
  std::vector<std::string> f;
  std::string cur;
  for (char c : line) {
    if (c == '\t') {
      f.push_back(cur);
      cur.clear();
    } else
      cur += c;
  }
  f.push_back(cur);
  K k;
  k.chrom = f[0];
  k.pos = std::stol(f[1]);
  k.mm = std::stoi(f[5]);
  const int bd = std::stoi(f[7]), br = std::stoi(f[8]);
  k.bg = bd + br;
  k.te = k.mm + k.bg;
  k.sc = (f[9] == "NA") ? -1.0 : std::stod(f[9]);
  return k;
}
static bool before_edit(const K &a, const K &b) {
  if (a.te != b.te)
    return a.te < b.te;
  if (a.mm != b.mm)
    return a.mm < b.mm;
  if (a.bg != b.bg)
    return a.bg < b.bg;
  return a.sc > b.sc;
}
static bool before_coord(const K &a, const K &b) {
  if (a.chrom != b.chrom)
    return a.chrom < b.chrom;
  return a.pos < b.pos;
}

int main() {
  std::cout << "=== test_result_merger ===\n\n";

  const std::string p0 = "/tmp/cz_rm_sh0.tsv", p1 = "/tmp/cz_rm_sh1.tsv",
                    p2 = "/tmp/cz_rm_sh2.tsv";
  const std::vector<std::string> s0 = {row("chr2", 100, '+', 2, 0, 0, "0.10"),
                                       row("chr1", 500, '-', 0, 0, 0, "0.95"),
                                       row("chr10", 50, '+', 1, 0, 0, "NA")};
  const std::vector<std::string> s1 = {row("chr1", 200, '+', 1, 1, 0, "0.40"),
                                       row("chr1", 200, '+', 1, 0, 0, "0.80"),
                                       row("chr3", 9, '-', 3, 0, 0, "0.00")};
  const std::vector<std::string> s2 = {row("chr1", 10, '+', 0, 0, 0, "NA"),
                                       row("chr1", 10, '+', 0, 0, 0, "0.50"),
                                       row("chr2", 100, '+', 1, 0, 0, "0.99")};
  const std::vector<std::string> paths = {p0, p1, p2};
  auto reset = [&] {
    write_shard(p0, s0);
    write_shard(p1, s1);
    write_shard(p2, s2);
  };

  std::multiset<std::string> all;
  for (const auto &v : {s0, s1, s2})
    for (const auto &r : v)
      all.insert(r);

  std::cout << "-- SortMode tokens --\n";
  record("sort_mode_from_string round-trips edit_distance",
         sort_mode_from_string("edit_distance") == SortMode::EditDistance);
  record("sort_mode_from_string round-trips coordinates",
         sort_mode_from_string("coordinates") == SortMode::Coordinates);
  bool threw = false;
  try {
    (void)sort_mode_from_string("bogus");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  record("sort_mode_from_string rejects unknown token", threw);

  std::cout << "\n-- EditDistance merge --\n";
  reset();
  std::size_t n =
      merge_sorted_shards(paths, "/tmp/cz_rm_edit.tsv", SortMode::EditDistance);
  record("edit: 9 rows written", n == 9u);
  auto outE = read_data("/tmp/cz_rm_edit.tsv");
  record("edit: output is permutation of inputs",
         std::multiset<std::string>(outE.begin(), outE.end()) == all);
  bool sortedE = true;
  for (std::size_t i = 1; i < outE.size(); ++i)
    if (before_edit(key(outE[i]), key(outE[i - 1])))
      sortedE = false;
  record("edit: globally sorted by comparator", sortedE);
  record("edit: shards removed (default)", !std::ifstream(p0).good());
  {
    std::vector<double> te0;
    for (const auto &l : outE) {
      auto k = key(l);
      if (k.te == 0)
        te0.push_back(k.sc);
    }
    record("edit: te0 group score-desc with NA last",
           te0.size() == 3 && te0[0] > te0[1] && te0[1] >= 0 && te0[2] < 0);
  }

  std::cout << "\n-- Coordinates merge --\n";
  reset();
  std::size_t n2 =
      merge_sorted_shards(paths, "/tmp/cz_rm_coord.tsv", SortMode::Coordinates);
  record("coord: 9 rows written", n2 == 9u);
  auto outC = read_data("/tmp/cz_rm_coord.tsv");
  record("coord: output is permutation of inputs",
         std::multiset<std::string>(outC.begin(), outC.end()) == all);
  bool sortedC = true;
  for (std::size_t i = 1; i < outC.size(); ++i)
    if (before_coord(key(outC[i]), key(outC[i - 1])))
      sortedC = false;
  record("coord: globally sorted by comparator", sortedC);
  std::size_t i10 = outC.size(), i2 = outC.size();
  for (std::size_t i = 0; i < outC.size(); ++i) {
    if (key(outC[i]).chrom == "chr10" && i10 == outC.size())
      i10 = i;
    if (key(outC[i]).chrom == "chr2" && i2 == outC.size())
      i2 = i;
  }
  record("coord: chr10 before chr2 (lexicographic)", i10 < i2);

  std::cout << "\n-- empty shard --\n";
  write_shard("/tmp/cz_rm_empty.tsv", {});
  write_shard("/tmp/cz_rm_one.tsv", {row("chr1", 1, '+', 0, 0, 0, "0.5")});
  std::size_t n3 =
      merge_sorted_shards({"/tmp/cz_rm_empty.tsv", "/tmp/cz_rm_one.tsv"},
                          "/tmp/cz_rm_final2.tsv", SortMode::EditDistance);
  record("empty-shard: 1 row written", n3 == 1u);

  for (const char *f :
       {"/tmp/cz_rm_edit.tsv", "/tmp/cz_rm_coord.tsv", "/tmp/cz_rm_final2.tsv"})
    std::remove(f);

  std::cout << "\n=== Results: " << g_passed << '/' << g_total << " passed";
  if (g_failed)
    std::cout << " (" << g_failed << " FAILED)";
  std::cout << " ===\n";
  return g_failed == 0 ? 0 : 1;
}