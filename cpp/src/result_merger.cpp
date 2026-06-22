/**
 * @file result_merger.cpp
 * @brief Implementation of per-shard sort + streaming k-way merge.
 */

#include "result_merger.hpp"

#include <algorithm>
#include <cstdio> // std::remove
#include <fstream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace crispritz {

namespace {

// Shard column indices (shared schema; see result_merger.hpp).
constexpr std::size_t COL_CHROM = 0;
constexpr std::size_t COL_POS = 1;
constexpr std::size_t COL_MM = 5;
constexpr std::size_t COL_BDNA = 7;
constexpr std::size_t COL_BRNA = 8;
constexpr std::size_t COL_CFD = 9;
constexpr std::size_t N_COLS = 10;

// Fallback header when no shard is available to copy one from. Must match
// ScoredTsvFormatter::header().
constexpr const char *CANONICAL_HEADER =
    "chrom\tpos\tstrand\tgrna\tspacer\tmismatches\t"
    "bulge_type\tbulge_dna\tbulge_rna\tcfd_score";

/**
 * @brief Parsed sort keys for one row, plus the raw line for output.
 */
struct Row {
  int total_edits = 0;
  int mm = 0;
  int bulges = 0;
  double score = -1.0; // NA -> -1.0 (sorts last under descending CFD)
  std::string chrom;
  long long pos = 0;
  std::string line; // full raw line (no trailing newline)
};

/** @brief Split a tab-separated line into fields (no allocation reuse). */
std::vector<std::string> split_tabs(const std::string &s) {
  std::vector<std::string> out;
  std::string field;
  for (char c : s) {
    if (c == '\t') {
      out.push_back(std::move(field));
      field.clear();
    } else {
      field += c;
    }
  }
  out.push_back(std::move(field));
  return out;
}

/** @brief Parse the sort keys from a raw data line. */
Row parse_row(std::string line) {
  const std::vector<std::string> f = split_tabs(line);
  if (f.size() < N_COLS)
    throw std::runtime_error("result_merger: malformed row (expected " +
                             std::to_string(N_COLS) + " columns, got " +
                             std::to_string(f.size()) + "): " + line);

  Row r;
  r.chrom = f[COL_CHROM];
  try {
    r.pos = std::stoll(f[COL_POS]);
    r.mm = std::stoi(f[COL_MM]);
    const int bd = std::stoi(f[COL_BDNA]);
    const int br = std::stoi(f[COL_BRNA]);
    r.bulges = bd + br;
    r.total_edits = r.mm + r.bulges;
  } catch (const std::exception &) {
    throw std::runtime_error("result_merger: non-numeric key field in row: " +
                             line);
  }
  const std::string &cfd = f[COL_CFD];
  r.score = (cfd == "NA") ? -1.0 : [&]() -> double {
    try {
      return std::stod(cfd);
    } catch (const std::exception &) {
      return -1.0; // unparseable score treated as NA (sorts last)
    }
  }();
  r.line = std::move(line);
  return r;
}

/** @brief True when @p a should come before @p b under @p mode. */
bool comes_before(const Row &a, const Row &b, SortMode mode) {
  if (mode == SortMode::Coordinates) {
    if (a.chrom != b.chrom)
      return a.chrom < b.chrom; // lexicographic
    return a.pos < b.pos;
  }
  // EditDistance
  if (a.total_edits != b.total_edits)
    return a.total_edits < b.total_edits;
  if (a.mm != b.mm)
    return a.mm < b.mm;
  if (a.bulges != b.bulges)
    return a.bulges < b.bulges;
  return a.score > b.score; // descending; NA (-1) sorts last
}

/** @brief Read a file's header line and all non-empty data lines. */
void read_shard(const std::string &path, std::string &header,
                std::vector<std::string> &rows) {
  std::ifstream in(path);
  if (!in.is_open())
    throw std::runtime_error("result_merger: cannot open shard: " + path);
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (first) {
      header = line; // header line (may be empty if file was empty)
      first = false;
      continue;
    }
    if (!line.empty())
      rows.push_back(line);
  }
}

/** @brief Sort one shard file in place (header preserved). */
void sort_shard_file(const std::string &path, SortMode mode) {
  std::string header;
  std::vector<std::string> raw;
  read_shard(path, header, raw);

  std::vector<Row> rows;
  rows.reserve(raw.size());
  for (std::string &l : raw)
    rows.push_back(parse_row(std::move(l)));

  std::stable_sort(
      rows.begin(), rows.end(),
      [mode](const Row &a, const Row &b) { return comes_before(a, b, mode); });

  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out.is_open())
    throw std::runtime_error("result_merger: cannot rewrite shard: " + path);
  if (!header.empty())
    out << header << '\n';
  for (const Row &r : rows)
    out << r.line << '\n';
  if (!out.good())
    throw std::runtime_error("result_merger: write failed for shard: " + path);
}

/** @brief A cursor over a sorted shard file: skips the header, yields rows. */
class ShardCursor {
public:
  explicit ShardCursor(const std::string &path)
      : in_(std::make_unique<std::ifstream>(path)) {
    if (!in_->is_open())
      throw std::runtime_error("result_merger: cannot open sorted shard: " +
                               path);
    std::string header;
    std::getline(*in_, header); // discard header
    advance();
  }

  [[nodiscard]] bool valid() const noexcept { return valid_; }
  [[nodiscard]] const Row &current() const noexcept { return current_; }

  void advance() {
    std::string line;
    while (std::getline(*in_, line)) {
      if (!line.empty() && line.back() == '\r')
        line.pop_back();
      if (line.empty())
        continue;
      current_ = parse_row(std::move(line));
      valid_ = true;
      return;
    }
    valid_ = false;
  }

private:
  std::unique_ptr<std::ifstream> in_;
  Row current_;
  bool valid_ = false;
};

} // namespace

// =============================================================================
// SortMode free functions
// =============================================================================

std::string_view to_string(SortMode mode) noexcept {
  switch (mode) {
  case SortMode::EditDistance:
    return "edit_distance";
  case SortMode::Coordinates:
    return "coordinates";
  }
  return "edit_distance";
}

SortMode sort_mode_from_string(std::string_view name) {
  if (name == "edit_distance")
    return SortMode::EditDistance;
  if (name == "coordinates")
    return SortMode::Coordinates;
  throw std::invalid_argument(
      "sort_mode_from_string: expected \"edit_distance\" or \"coordinates\", "
      "got \"" +
      std::string(name) + '"');
}

// =============================================================================
// merge_sorted_shards
// =============================================================================

std::size_t merge_sorted_shards(const std::vector<std::string> &shard_paths,
                                const std::string &final_path, SortMode mode,
                                bool write_header, bool remove_inputs) {
  // Phase 1: sort each shard in place (bounded to one shard in memory).
  for (const std::string &path : shard_paths)
    sort_shard_file(path, mode);

  // Determine the header to emit (copy from the first shard if possible).
  std::string header = CANONICAL_HEADER;
  if (!shard_paths.empty()) {
    std::ifstream first(shard_paths.front());
    std::string h;
    if (first.is_open() && std::getline(first, h) && !h.empty())
      header = h;
  }

  std::ofstream out(final_path, std::ios::out | std::ios::trunc);
  if (!out.is_open())
    throw std::runtime_error("result_merger: cannot open output: " +
                             final_path);
  if (write_header)
    out << header << '\n';

  // Phase 2: streaming k-way merge (one row per shard resident at a time).
  std::vector<ShardCursor> cursors;
  cursors.reserve(shard_paths.size());
  for (const std::string &path : shard_paths)
    cursors.emplace_back(path);

  // Min-heap of cursor indices ordered by their current row. The priority
  // queue is a max-heap, so the comparator returns true when a should come
  // AFTER b — making the "first" row the top of the heap.
  const auto heap_cmp = [&](std::size_t a, std::size_t b) {
    const Row &ra = cursors[a].current();
    const Row &rb = cursors[b].current();
    if (comes_before(ra, rb, mode))
      return false; // a is not "after" b
    if (comes_before(rb, ra, mode))
      return true; // a is after b
    return a > b;  // tie: stable by shard index
  };
  std::priority_queue<std::size_t, std::vector<std::size_t>, decltype(heap_cmp)>
      heap(heap_cmp);

  for (std::size_t i = 0; i < cursors.size(); ++i)
    if (cursors[i].valid())
      heap.push(i);

  std::size_t written = 0;
  while (!heap.empty()) {
    const std::size_t i = heap.top();
    heap.pop();
    out << cursors[i].current().line << '\n';
    ++written;
    cursors[i].advance();
    if (cursors[i].valid())
      heap.push(i);
  }

  out.flush();
  if (!out.good())
    throw std::runtime_error("result_merger: write failed for output: " +
                             final_path);

  // Close cursors before removing the files they read from.
  cursors.clear();
  if (remove_inputs)
    for (const std::string &path : shard_paths)
      std::remove(path.c_str());

  return written;
}

} // namespace crispritz