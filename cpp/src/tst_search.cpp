#include "tst_search.hpp"

#include "tst_utils.hpp" // iupac::*, pack/unpack nibbles, sentinels, reverse_complement

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

// using namespace pam; // NucleotideEncoder::*

namespace crispritz {

// =========================================================================
// LoadedTST
// =========================================================================

LoadedTST::LoadedTST(std::vector<TSTNode> nodes, std::vector<TSTLeaf> leaves,
                     int guide_length, int pam_limit, std::string source_path)
    : nodes_(std::move(nodes)), leaves_(std::move(leaves)),
      guide_length_(guide_length), pam_limit_(pam_limit), source_path_(std::move(source_path)) {
  if (nodes_.empty())
    throw std::invalid_argument("LoadedTST: node pool must not be empty");
  if (guide_length_ <= 0)
    throw std::invalid_argument("LoadedTST: guide_length must be > 0, got " +
                                std::to_string(guide_length_));
  if (pam_limit_ <= 0)
    throw std::invalid_argument("LoadedTST: pam_limit must be > 0, got " +
                                std::to_string(pam_limit_));
}

// =========================================================================
// Partition deserialization
// =========================================================================
//
// .bin layout (see TernarySearchTree::write_partition):
//   [4 B]  magic   (TST_BIN_MAGIC)
//   [4 B]  version (TST_BIN_VERSION)
//   [4 B]  chunk_size   (number of leaves)
//   [4 B]  guide_length
//   [4 B]  pam_limit    (PAM length; PAM byte width is ceil(pam_limit/2))
//   for each leaf:
//     [4 B]                 guide_index (signed)
//     [ceil(pam_limit/2) B] bit-packed PAM nibbles
//     [1 B]                 '0'            -> next == 0
//                           '_' + [4 B]    -> next index
//   [4 B]  node_count
//   [var]  nibble-packed node stream (pre-order: split, lo, hi, eq)
//
// pam_limit now travels in the header, so the PAM byte width is read, not
// guessed from the filename. Every read below is bounds-checked against the
// end of the in-memory buffer, and node reconstruction is capped at the
// header's node_count, so a corrupt or truncated partition raises a
// descriptive std::runtime_error instead of walking off the buffer or
// recursing into a stack overflow (which the OS reports as SIGBUS).
// =========================================================================

namespace {
struct NibbleReader; // defined below
void rebuild_children(NibbleReader &r, std::vector<TSTNode> &nodes,
                      int self_idx, int max_nodes);

/** @brief Bounds-checked little-endian uint32 read; advances the cursor. */
std::uint32_t read_u32(const unsigned char *&p, const unsigned char *end,
                       const char *what) {
  if (static_cast<std::size_t>(end - p) < sizeof(std::uint32_t))
    throw std::runtime_error(
        std::string("load_partition: truncated file while reading ") + what);
  std::uint32_t v;
  std::memcpy(&v, p, sizeof(v));
  p += sizeof(v);
  return v;
}

/** @brief Bounds-checked little-endian int32 read; advances the cursor. */
int read_i32(const unsigned char *&p, const unsigned char *end,
             const char *what) {
  if (static_cast<std::size_t>(end - p) < sizeof(int))
    throw std::runtime_error(
        std::string("load_partition: truncated file while reading ") + what);
  int v;
  std::memcpy(&v, p, sizeof(int));
  p += sizeof(int);
  return v;
}

/**
 * @brief Cursor over the nibble-packed node stream.
 *
 * Every access is bounds-checked against @c end, so a desynced or truncated
 * stream throws instead of reading past the buffer.
 */
struct NibbleReader {
  const unsigned char *data;
  const unsigned char *end;
  std::size_t pos;  // byte offset within the buffer
  bool high = true; // next nibble is the high half

  NibbleReader(const unsigned char *d, std::size_t start,
               const unsigned char *e)
      : data(d), end(e), pos(start) {}

  uint8_t next_nibble() {
    if (data + pos >= end)
      throw std::runtime_error("load_partition: node stream overran the file "
                               "(corrupt or truncated partition)");
    const uint8_t byte = data[pos];
    const uint8_t nib = high ? high_nibble(byte) : low_nibble(byte);
    if (high) {
      high = false;
    } else {
      high = true;
      ++pos;
    }
    return nib;
  }

  /** @brief Align to the next whole byte (before reading a 4-byte ptr). */
  void align_byte() {
    if (!high) {
      high = true;
      ++pos;
    }
  }

  int read_leaf_ptr() {
    align_byte();
    if (static_cast<std::size_t>(end - (data + pos)) < sizeof(int))
      throw std::runtime_error("load_partition: leaf pointer overran the file "
                               "(corrupt or truncated partition)");
    int v;
    std::memcpy(&v, data + pos, sizeof(int));
    pos += sizeof(int);
    return v;
  }
};

/** @brief Append a node, but never exceed the count declared in the header. */
int alloc_checked(std::vector<TSTNode> &nodes, int max_nodes) {
  if (static_cast<int>(nodes.size()) >= max_nodes)
    throw std::runtime_error("load_partition: node stream produced more nodes "
                             "than the header declared (corrupt partition)");
  const int idx = static_cast<int>(nodes.size());
  nodes.emplace_back();
  return idx;
}

int rebuild_node(NibbleReader &r, std::vector<TSTNode> &nodes, int max_nodes) {
  const uint8_t split_nib = r.next_nibble();
  const int self_idx = alloc_checked(nodes, max_nodes);
  nodes[self_idx].splitchar_enc = split_nib;
  nodes[self_idx].splitchar = iupac::decode_genome(split_nib);

  { // lokid
    const uint8_t nib = r.next_nibble();
    if (nib == NULL_CHILD_NIBBLE) {
      nodes[self_idx].lokid = 0;
    } else {
      const int child = alloc_checked(nodes, max_nodes);
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = iupac::decode_genome(nib);
      rebuild_children(r, nodes, child, max_nodes);
      nodes[self_idx].lokid = child;
    }
  }
  { // hikid
    const uint8_t nib = r.next_nibble();
    if (nib == NULL_CHILD_NIBBLE) {
      nodes[self_idx].hikid = 0;
    } else {
      const int child = alloc_checked(nodes, max_nodes);
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = iupac::decode_genome(nib);
      rebuild_children(r, nodes, child, max_nodes);
      nodes[self_idx].hikid = child;
    }
  }
  { // eqkid
    const uint8_t nib = r.next_nibble();
    if (nib == SENTINEL_NIBBLE) {
      nodes[self_idx].eqkid = r.read_leaf_ptr();
    } else {
      const int child = alloc_checked(nodes, max_nodes);
      nodes[child].splitchar_enc = nib;
      nodes[child].splitchar = iupac::decode_genome(nib);
      rebuild_children(r, nodes, child, max_nodes);
      nodes[self_idx].eqkid = child;
    }
  }
  return self_idx;
}

void rebuild_children(NibbleReader &r, std::vector<TSTNode> &nodes,
                      int self_idx, int max_nodes) {
  // lokid
  const uint8_t lo = r.next_nibble();
  if (lo == NULL_CHILD_NIBBLE) {
    nodes[self_idx].lokid = 0;
  } else {
    const int child = alloc_checked(nodes, max_nodes);
    nodes[child].splitchar_enc = lo;
    nodes[child].splitchar = iupac::decode_genome(lo);
    rebuild_children(r, nodes, child, max_nodes);
    nodes[self_idx].lokid = child;
  }

  // hikid
  const uint8_t hi = r.next_nibble();
  if (hi == NULL_CHILD_NIBBLE) {
    nodes[self_idx].hikid = 0;
  } else {
    const int child = alloc_checked(nodes, max_nodes);
    nodes[child].splitchar_enc = hi;
    nodes[child].splitchar = iupac::decode_genome(hi);
    rebuild_children(r, nodes, child, max_nodes);
    nodes[self_idx].hikid = child;
  }

  // eqkid
  const uint8_t eq = r.next_nibble();
  if (eq == SENTINEL_NIBBLE) {
    nodes[self_idx].eqkid = r.read_leaf_ptr();
  } else {
    const int child = alloc_checked(nodes, max_nodes);
    nodes[child].splitchar_enc = eq;
    nodes[child].splitchar = iupac::decode_genome(eq);
    rebuild_children(r, nodes, child, max_nodes);
    nodes[self_idx].eqkid = child;
  }
}

} // namespace

LoadedTST load_partition(const std::string &partition_path) {
  std::ifstream in(partition_path, std::ios::in | std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("load_partition: cannot open " + partition_path);

  std::vector<unsigned char> buf((std::istreambuf_iterator<char>(in)),
                                 std::istreambuf_iterator<char>());
  in.close();

  if (buf.size() < 20u) // magic+version+chunk_size+guide_length+pam_limit
    throw std::runtime_error("load_partition: file too small: " +
                             partition_path);

  const unsigned char *p = buf.data();
  const unsigned char *const end = buf.data() + buf.size();

  const std::uint32_t magic = read_u32(p, end, "magic");
  if (magic != TST_BIN_MAGIC)
    throw std::runtime_error(
        "load_partition: " + partition_path +
        " is not a CRISPRitz-plus index (bad magic). Rebuild the index with "
        "the current index-genome.");

  const std::uint32_t version = read_u32(p, end, "version");
  if (version != TST_BIN_VERSION)
    throw std::runtime_error(
        "load_partition: " + partition_path + " uses index format v" +
        std::to_string(version) + " but this build expects v" +
        std::to_string(TST_BIN_VERSION) + ". Rebuild the index.");

  const int chunk_size = read_i32(p, end, "chunk_size");
  const int guide_length = read_i32(p, end, "guide_length");
  const int pam_limit = read_i32(p, end, "pam_limit");

  if (chunk_size < 0 || guide_length <= 0 || pam_limit <= 0)
    throw std::runtime_error("load_partition: invalid header in " +
                             partition_path);

  const int pam_bytes = (pam_limit + 1) / 2; // ceil(pam_limit / 2)

  // --- leaf section -----------------------------------------------------
  std::vector<TSTLeaf> leaves;
  leaves.reserve(static_cast<std::size_t>(chunk_size));

  for (int i = 0; i < chunk_size; ++i) {
    TSTLeaf leaf;
    leaf.guide_index = read_i32(p, end, "leaf.guide_index");

    if (static_cast<std::size_t>(end - p) < static_cast<std::size_t>(pam_bytes))
      throw std::runtime_error("load_partition: truncated PAM bytes in " +
                               partition_path);
    leaf.pam_seq_enc.assign(p, p + pam_bytes);
    p += pam_bytes;

    if (p >= end)
      throw std::runtime_error("load_partition: truncated leaf marker in " +
                               partition_path);
    const unsigned char marker = *p++;
    if (marker == static_cast<unsigned char>('_'))
      leaf.next = read_i32(p, end, "leaf.next");
    else
      leaf.next = 0; // marker == '0'

    leaves.push_back(std::move(leaf));
  }

  // --- node section -----------------------------------------------------
  const int node_count = read_i32(p, end, "node_count");
  if (node_count < 0)
    throw std::runtime_error("load_partition: negative node_count in " +
                             partition_path);

  std::vector<TSTNode> nodes;
  nodes.reserve(static_cast<std::size_t>(node_count));

  if (node_count > 0) {
    const std::size_t node_start = static_cast<std::size_t>(p - buf.data());
    NibbleReader reader(buf.data(), node_start, end);
    rebuild_node(reader, nodes, node_count);
  } else {
    // Degenerate but valid: a single sentinel root so LoadedTST's
    // invariant (non-empty pool) holds.
    nodes.emplace_back();
  }

  return LoadedTST{std::move(nodes), std::move(leaves), guide_length,
                   pam_limit, partition_path}; 
}

// =========================================================================
// TSTSearcher
// =========================================================================

TSTSearcher::TSTSearcher(SearchConfiguration config)
    : config_(std::move(config)) {
  // Note: config_.threads() is intentionally NOT consulted here. The
  // searcher is single-threaded by design; thread count belongs to the
  // Python orchestration layer, which sizes its ThreadPoolExecutor and
  // dispatches one search_partition() call per .bin file. The searcher
  // only reads the edit-budget fields (mismatches / bulges). See the
  // serial-loop rationale in search_all().
}

namespace {
/**
 * @brief Mutable-by-value recursion frame for the bounded near-search.
 *
 * Passing the state by value (not by reference) is deliberate: each
 * branch of the traversal gets its own copy, so a mismatch or bulge
 * spent down one path cannot corrupt the budget of a sibling path. This
 * removes an entire class of backtracking bugs at the cost of copying a
 * few ints and two short strings per recursive call.
 */
struct Frame {
  int mm_left;            // remaining substitution budget
  int bdna_left;          // remaining DNA-bulge budget
  int brna_left;          // remaining RNA-bulge budget
  std::size_t guide_pos;  // next index into the query guide
  std::string aln_guide;  // aligned guide chars accumulated so far
  std::string aln_target; // aligned target chars accumulated so far
  bool dna_bulge_used = false;   // DNA bulge was opened on this path
  bool rna_bulge_used = false;   // RNA bulge was opened on this path
};

/** @brief Bundle of read-only context shared across the recursion. */
struct Context {
  const LoadedTST *tst;
  std::string_view guide; // full query guide
  int guide_len;          // guide.size()
  std::string_view chrom; // chromosome name for emitted hits
  BulgeMode bulge_mode = BulgeMode::MixedBulges;
  Strand strand_for_index(int guide_index) const {
    return guide_index >= 0 ? Strand::Forward : Strand::Reverse;
  }
};

// Forward declaration: the traversal recurses through three node slots.
void traverse(const Context &ctx, int node_idx, Frame frame,
              std::vector<OffTarget> &out, const std::string &pam, bool pam_at_start);

/**
 * @brief Decode a bit-packed PAM (TSTLeaf::pam_seq_enc) into nucleotides.
 *
 * Exact inverse of TernarySearchTree::encode_pam_bytes: two 4-bit IUPAC codes
 * per byte, high nibble first, decoded via iupac::decode_genome. The result
 * is returned in *storage order* (the order the builder packed it); recovering
 * genomic 5'->3' orientation depends on strand and PAM placement and is left
 * to the caller (see emit_leaf_chain).
 *
 * @param enc        Packed PAM bytes; must hold >= ceil(pam_limit/2) bytes.
 * @param pam_limit  Number of PAM nucleotides to decode (> 0).
 * @return           Decoded nucleotide string of length @p pam_limit.
 * @throws std::runtime_error if @p enc is too short for @p pam_limit.
 * @complexity O(pam_limit).
 */
std::string decode_pam(const std::vector<uint8_t> &enc, int pam_limit) {
  if (pam_limit <= 0)
    return {};
  const std::size_t need = static_cast<std::size_t>((pam_limit + 1) / 2);
  if (enc.size() < need)
    throw std::runtime_error(
        "decode_pam: encoded PAM shorter than index pam_limit");
  std::string out;
  out.reserve(static_cast<std::size_t>(pam_limit));
  for (int i = 0; i < pam_limit; ++i) {
    const uint8_t byte = enc[static_cast<std::size_t>(i) >> 1];
    const uint8_t nib = (i & 1) ? low_nibble(byte) : high_nibble(byte);
    out.push_back(iupac::decode_genome(nib));
  }
  return out;
}

/**
 * @brief Emit OffTargets for the leaf chain rooted at @p leaf_ptr.
 *
 * @p leaf_ptr is the negative encoding -(within_chunk_index + 1) stored
 * in a terminal node's eqkid. All leaves chained via TSTLeaf::next share
 * the same guide path, so they share the same aligned strings and edit
 * counts; only their genomic position/strand differ.
 */
void emit_leaf_chain(const Context &ctx, int leaf_ptr, const Frame &frame,
                     std::vector<OffTarget> &out, const std::string &pam, bool pam_at_start) {

  const auto &leaves = ctx.tst->leaves();
  int idx = -leaf_ptr - 1; // decode within-chunk leaf index

  while (idx >= 0 && idx < static_cast<int>(leaves.size())) {
    const TSTLeaf &leaf = leaves[static_cast<std::size_t>(idx)];

    const Strand strand = ctx.strand_for_index(leaf.guide_index);
    // guide_index stores the 0-based leftmost genomic coordinate of
    // the window (same convention for both strands; the sign only
    // encodes strand). OffTarget::pos is 1-based, so add 1.
    const int leftmost0 =
        leaf.guide_index >= 0 ? leaf.guide_index : -leaf.guide_index;
    const int pos = leftmost0 + 1;

    // Edit-distance breakdown = initial budget minus remaining.
    // The caller seeds the frame with the configured maxima, so the
    // spent amount is (max - left). We recover the maxima from the
    // accumulated alignment instead: count lowercase (mismatch) and
    // '-' (bulge) characters, which are written during traversal.
    int mm_count = 0, bdna_count = 0, brna_count = 0;
    for (std::size_t k = 0; k < frame.aln_target.size(); ++k) {
      char tg = frame.aln_target[k];
      char gd = (k < frame.aln_guide.size()) ? frame.aln_guide[k] : '?';
      if (gd == '-')
        ++bdna_count; // gap in guide  = DNA bulge
      else if (tg == '-')
        ++brna_count; // gap in target = RNA bulge
      else if (tg >= 'a' && tg <= 'z')
        ++mm_count; // lowercase = mismatch
    }

    std::string aln_guide = static_cast<std::string>(frame.aln_guide);
    std::string aln_target = static_cast<std::string>(frame.aln_target);

    // Index stores PAM-at-end guide/target in 3'->5' (reversed) order; flip
    // the alignment back to genomic 5'->3' for output.
    if (!pam_at_start) {
      std::reverse(aln_guide.begin(), aln_guide.end());
      std::reverse(aln_target.begin(), aln_target.end());
    }

    // Decode the ACTUAL genomic PAM and re-orient to 5'->3'. encode_pam_bytes
    // stores the PAM reversed in every case EXCEPT reverse-strand + PAM-at-end
    // (see extract_forward / extract_reverse), so undo exactly that reversal.
    std::string pam_target = decode_pam(leaf.pam_seq_enc, ctx.tst->pam_limit());
    const bool pam_stored_reversed =
        (strand == Strand::Forward) || pam_at_start;
    if (pam_stored_reversed)
      std::reverse(pam_target.begin(), pam_target.end());

    // Guide carries the PAM *motif* (e.g. "NGG"); target carries the decoded
    // actual PAM bases at this genomic site.
    if (!pam_at_start) {
      aln_guide  += pam;
      aln_target += pam_target;
    } else {
      aln_guide  = pam + aln_guide;
      aln_target = pam_target + aln_target;
    }

    // Adjust target position according to strand orientation and PAM
    int pos_t = pos;
    if (!pam_at_start && strand == Strand::Forward) {
      pos_t = pos_t - (aln_target.size() - brna_count) + 1;
    } else if (pam_at_start && strand == Strand::Forward)
    {
      pos_t = pos_t - (aln_target.size() - brna_count) + 1;
    }
    

    out.emplace_back(
        /*chrom  */ std::string(ctx.chrom),
        /*pos    */ pos_t,
        /*strand */ strand,
        /*grna   */ aln_guide,
        /*target */ aln_target,
        /*mm     */ mm_count,
        /*bdna   */ bdna_count,
        /*brna   */ brna_count);

    // TSTLeaf::next holds the *encoded* pointer to the next leaf in
    // the collision chain — the same -(within_chunk + 1) convention
    // used by the head pointer (see TernarySearchTree::insert, which
    // assigns next = node.eqkid). 0 terminates the chain.
    if (leaf.next == 0)
      break;
    idx = -leaf.next - 1;
  }
}

// Mirrors legacy saveIndices: collect every leaf below node_idx, emitting
// the alignment accumulated so far. The suffix bases below this node are
// not part of the guide match, so the alignment strings are NOT extended.
void harvest(const Context &ctx, int node_idx, const Frame &frame,
             std::vector<OffTarget> &out, const std::string &pam, bool pam_at_start) {
  if (node_idx <= 0 || node_idx >= (int)ctx.tst->nodes().size()) 
    return;
  const TSTNode &n = ctx.tst->nodes()[(std::size_t)node_idx];
  if (n.lokid > 0) 
    harvest(ctx, n.lokid, frame, out, pam, pam_at_start);
  if (n.hikid > 0) 
    harvest(ctx, n.hikid, frame, out, pam, pam_at_start);
  if (n.eqkid < 0)      
    emit_leaf_chain(ctx, n.eqkid, frame, out, pam, pam_at_start);
  else if (n.eqkid > 0) 
    harvest(ctx, n.eqkid, frame, out, pam, pam_at_start);
}

void traverse(const Context &ctx, int node_idx, Frame frame,
              std::vector<OffTarget> &out, const std::string &pam, bool pam_at_start) {
  if (node_idx < 0 || node_idx >= static_cast<int>(ctx.tst->nodes().size()))
    return;

  const TSTNode &node = ctx.tst->nodes()[static_cast<std::size_t>(node_idx)];

  // A node with no splitchar encoding and no children is the reserved
  // root sentinel (index 0). Descend into its eqkid only.
  const bool is_root_sentinel = (node_idx == 0 && node.splitchar_enc == 0 &&
                                 node.lokid == 0 && node.hikid == 0);

  if (is_root_sentinel) {
    if (node.eqkid < 0)
      emit_leaf_chain(ctx, node.eqkid, frame, out, pam, pam_at_start);
    else if (node.eqkid > 0)
      traverse(ctx, node.eqkid, frame, out, pam, pam_at_start);
    return;
  }

  // alternative splitchars at this depth — no guide consumed
  if (node.lokid > 0) 
    traverse(ctx, node.lokid, frame, out, pam, pam_at_start);
  if (node.hikid > 0) 
    traverse(ctx, node.hikid, frame, out, pam, pam_at_start);

  // GUIDE EXHAUSTED -> harvest leaves in this subtree (was: return)
  if (frame.guide_pos >= static_cast<std::size_t>(ctx.guide_len)) {
    if (node.eqkid < 0)      
      emit_leaf_chain(ctx, node.eqkid, frame, out, pam, pam_at_start);
    else if (node.eqkid > 0) 
      harvest(ctx, node.eqkid, frame, out, pam, pam_at_start);
    return;
  }

  // Not exhausted: read the query base. Index stores guides reversed
  // (PAM-at-end), so consume the query from the 3' end inward
  const int gpos = ctx.guide_len - 1 - static_cast<int>(frame.guide_pos);
  const char q_char  = ctx.guide[static_cast<size_t>(gpos)];
  const uint8_t q_enc = iupac::encode_genome(q_char);
  const bool is_match = iupac::matches(q_enc, node.splitchar_enc);
  const char node_char = node.splitchar;

  // step into the equal child; emit directly only if it's a terminal AND
  // this step just exhausted the guide (handles a max_bulges_==0 index too)
  auto step_eq = [&](Frame nx) {
    if (node.eqkid > 0)
      traverse(ctx, node.eqkid, std::move(nx), out, pam, pam_at_start);
    else if (node.eqkid < 0 && nx.guide_pos >= static_cast<std::size_t>(ctx.guide_len))
      emit_leaf_chain(ctx, node.eqkid, nx, out, pam, pam_at_start);
  };

  if (is_match) {
    Frame nx = frame; 
    nx.aln_guide += node_char; 
    nx.aln_target += node_char; 
    nx.guide_pos++;
    step_eq(std::move(nx));
  } else if (frame.mm_left > 0) {
    Frame nx = frame; 
    nx.mm_left--; 
    nx.aln_guide += q_char;// query base (upper)
    nx.aln_target += static_cast<char>( // genome base (lower)
        node_char >= 'A' && node_char <= 'Z' ? node_char - 'A' + 'a'
                                             : node_char);
    nx.guide_pos++;
    step_eq(std::move(nx));
  }

  // ---- DNA bulge: extra base in the genome/target, gap in the guide.
  //      Consume this node's splitchar into the target, advance the
  //      node (eqkid) but NOT the guide position.
  if (frame.bdna_left > 0 &&
      !(ctx.bulge_mode == BulgeMode::SingleBulgeType && frame.rna_bulge_used)) {
    Frame nx = frame; 
    nx.bdna_left--; 
    nx.dna_bulge_used = true;
    nx.aln_guide += '-';        // gap in guide
    nx.aln_target += node_char; // extra genomic base
    // guide_pos unchanged
    if (node.eqkid > 0) 
      traverse(ctx, node.eqkid, std::move(nx), out, pam, pam_at_start);
    // (a DNA bulge cannot land exactly on a leaf terminal without a
    //  following matched base, so we do not emit on eqkid < 0 here)
  }

  // ---- RNA bulge: extra base in the guide, gap in the genome/target.
  //      Consume a guide character into the alignment as a gap on the
  //      target side, advance the guide but stay on the same node.
  if (frame.brna_left > 0 &&
      !(ctx.bulge_mode == BulgeMode::SingleBulgeType && frame.dna_bulge_used)) {
    Frame nx = frame; 
    nx.brna_left--; 
    nx.rna_bulge_used = true;
    nx.aln_guide += q_char; // extra guide base
    nx.aln_target += '-';   // gap in target
    nx.guide_pos++;
    traverse(ctx, node_idx, std::move(nx), out, pam, pam_at_start); // same node, advanced guide
  }
}

} // namespace

std::vector<OffTarget> TSTSearcher::search(const LoadedTST &tst,
                                           std::string_view guide_seq,
                                           const std::string &chrom,
                                           const std::string &pam,
                                           bool pam_at_start,
                                           BulgeMode bulge_mode) const {  
  // Deferred guide-length / edit-budget validation (see header contract).
  if (config_.max_total_edits() > tst.guide_length())
    throw std::invalid_argument("TSTSearcher::search: max_total_edits (" +
                                std::to_string(config_.max_total_edits()) +
                                ") exceeds index guide_length (" +
                                std::to_string(tst.guide_length()) + ')');

  if (static_cast<int>(guide_seq.size()) != tst.guide_length())
    throw std::invalid_argument("TSTSearcher::search: guide length (" +
                                std::to_string(guide_seq.size()) +
                                ") does not match index guide_length (" +
                                std::to_string(tst.guide_length()) + ')');

  Context ctx;
  ctx.tst = &tst;
  ctx.guide = guide_seq;
  ctx.guide_len = static_cast<int>(guide_seq.size());
  ctx.chrom = chrom;
  ctx.bulge_mode = bulge_mode; 

  Frame root;
  root.mm_left = config_.max_mismatches();
  root.bdna_left = config_.max_bulges_dna();
  root.brna_left = config_.max_bulges_rna();
  root.guide_pos = 0;
   // dna_bulge_used / rna_bulge_used default to false
  root.aln_guide.reserve(
      static_cast<std::size_t>(ctx.guide_len + config_.max_bulges_total()));
  root.aln_target.reserve(
      static_cast<std::size_t>(ctx.guide_len + config_.max_bulges_total()));

  std::vector<OffTarget> out;
  traverse(ctx, /*root index*/ 0, std::move(root), out, pam, pam_at_start);
  return out;
}

SearchResult TSTSearcher::search_all(const LoadedTST &tst,
                                     const std::vector<std::string> &guides,
                                     const std::string &chrom,
                                     const std::string &pam, bool pam_at_start,
                                     BulgeMode bulge_mode) const {
  SearchResult result;
  result.source_path = tst.source_path();
  result.hits_by_guide.reserve(guides.size());

  // The loop over guides is DELIBERATELY serial. Do not parallelise it
  // (e.g. with "#pragma omp parallel for").
  //
  // The parallelism axis in this architecture is the *partition*, not the
  // guide: the Python orchestration layer drives a memory-aware
  // ThreadPoolExecutor that calls search_partition() once per .bin file,
  // and the pybind11 binding releases the GIL so those calls run as true
  // parallel C++. Real workloads have very few guides (1-20) but many
  // partitions, so partition-level parallelism saturates the cores while
  // guide-level parallelism would not.
  //
  // Adding OpenMP here would nest guide-parallelism inside the Python
  // thread pool and oversubscribe the CPU (threads x guides), degrading
  // performance. Keeping this serial is the correct design, not an
  // oversight.
  for (const std::string &g : guides)
    result.hits_by_guide.push_back(search(tst, g, chrom, pam, pam_at_start, bulge_mode));

  return result;
}

// =========================================================================
// search_partition — atomic unit of parallel work
// =========================================================================

SearchResult search_partition(const std::string &partition_path,
                              const std::string &chrom,
                              const std::vector<std::string> &guides,
                              const SearchConfiguration &config,
                              const std::string &pam, bool pam_at_start,
                              BulgeMode bulge_mode) {
  LoadedTST tst = load_partition(partition_path);
  TSTSearcher searcher(config);
  return searcher.search_all(tst, guides, chrom, pam, pam_at_start, bulge_mode);
}

} // namespace crispritz