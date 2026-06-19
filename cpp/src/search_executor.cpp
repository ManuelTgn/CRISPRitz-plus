/**
 * @file search_executor.cpp
 * @brief Implementation of the per-partition streaming search executor.
 *
 * See search_executor.hpp for the architectural contract. This file performs
 * no concurrency of its own: it is the body of one parallel task, driven by
 * the Python orchestration layer.
 */

#include "search_executor.hpp"

#include "output_writer.hpp" // OutputWriter, OutputWriter::Session
#include "profile_data.hpp"  // ProfileAccumulator
#include "tst_search.hpp"    // LoadedTST, TSTSearcher, load_partition

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace crispritz {

// =========================================================================
// ScoredTsvFormatter
// =========================================================================
//
// Canonical 10-column shard schema (shared with scores/shard_scoring.py):
//   chrom, pos, strand, grna, spacer, mismatches,
//   bulge_type, bulge_dna, bulge_rna, cfd_score
//
// The C++ OffTarget::target() supplies the "spacer" column. cfd_score is the
// "NA" placeholder; the Python scorer fills it in place.
// =========================================================================

std::string ScoredTsvFormatter::header() const {
  return "chrom\tpos\tstrand\tgrna\tspacer\tmismatches\t"
         "bulge_type\tbulge_dna\tbulge_rna\tcfd_score";
}

std::string ScoredTsvFormatter::format_row(const OffTarget &ot) const {
  std::string row;
  row.reserve(ot.chrom().size() + ot.grna().size() + ot.target().size() + 64u);

  row += ot.chrom();
  row += '\t';
  row += std::to_string(ot.pos());
  row += '\t';
  row += to_char(ot.strand());
  row += '\t';
  row += ot.grna();
  row += '\t';
  row += ot.target(); // -> "spacer" column
  row += '\t';
  row += std::to_string(ot.mismatches());
  row += '\t';
  row += ot.bulge_type();
  row += '\t';
  row += std::to_string(ot.bulge_dna());
  row += '\t';
  row += std::to_string(ot.bulge_rna());
  row += '\t';
  row += SCORE_NA; // filled in place by the Python per-shard scorer
  return row;
}

// =========================================================================
// run_search_executor
// =========================================================================

PartitionResult run_search_executor(const std::string &partition_path,
                                    const std::string &chrom,
                                    const std::vector<std::string> &guides,
                                    const SearchConfiguration &config,
                                    int pam_len, bool pam_at_start,
                                    const std::string &shard_path) {
  LoadedTST tst = load_partition(partition_path);
  TSTSearcher searcher(config);

  const bool want_targets = config.write_targets();
  const bool want_profile = config.write_profile();

  PartitionResult result;
  result.source_path = tst.source_path();

  // Targets sink: a threshold-flushed shard Session with the scored-TSV
  // formatter. The OutputWriter owns the formatter and must outlive the
  // Session, so it is declared first (destroyed last).
  std::optional<OutputWriter> writer;
  std::optional<OutputWriter::Session> session;
  if (want_targets) {
    writer.emplace(std::make_unique<ScoredTsvFormatter>());
    session.emplace(
        writer->open_session_to_file(shard_path, SHARD_FLUSH_THRESHOLD));
    result.shard_path = shard_path;
  }

  // Profile sinks: one accumulator per guide. Each query guide IS the guide
  // body (its length equals the index guide length; PAM was applied at index
  // time), so guide.size() is the body length the accumulator expects.
  std::vector<ProfileAccumulator> accumulators;
  if (want_profile) {
    accumulators.reserve(guides.size());
    for (const std::string &guide : guides) {
      accumulators.emplace_back(guide, static_cast<int>(guide.size()), pam_len,
                                config.max_mismatches(),
                                config.max_bulges_dna(),
                                config.max_bulges_rna(), pam_at_start);
    }
  }

  // Search one guide at a time and sink each guide's hits immediately, so the
  // hit vector is released before the next guide is searched.
  for (std::size_t i = 0; i < guides.size(); ++i) {
    const std::vector<OffTarget> hits = searcher.search(tst, guides[i], chrom);
    for (const OffTarget &ot : hits) {
      if (want_targets)
        session->add(ot);
      if (want_profile)
        accumulators[i].push(ot);
    }
    result.total_hits += hits.size();
  }

  // Finalize the shard (explicit close so a final-flush failure surfaces as an
  // exception rather than being swallowed by the destructor).
  if (session)
    result.rows_written = session->close();

  // Build the per-guide profiles for this partition (merged across partitions
  // by a later step before write_all_profiles()).
  if (want_profile) {
    result.profiles.reserve(accumulators.size());
    for (const ProfileAccumulator &acc : accumulators)
      result.profiles.push_back(acc.build());
  }

  return result;
}

} // namespace crispritz