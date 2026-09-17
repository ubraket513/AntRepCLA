// Construction of the Hamming graph over CDR3 sequences.
//
// Two distinct CDR3s are joined when they have the same length, share at least
// one V gene and one J gene, and differ in at most `threshold_ratio` of their
// length. The first three conditions are equality constraints and so act as
// blocking keys; only the distance is measured inside a block.
//
// Blocking alone is not enough: the key space is bounded (lengths x V genes x
// J genes), so block membership grows linearly with the repertoire and pairs
// within a block grow quadratically. A pigeonhole filter inside each block
// cuts that down without approximating -- two sequences within d mismatches
// cannot differ in all of d + 1 contiguous chunks, so only pairs colliding in
// a chunk index are measured.
#pragma once

#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

#include "types.hpp"

namespace ar {

using Edge = std::pair<NodeId, NodeId>;

// Below this many within-block pairs, a serial scan finishes before threads
// would have paid for themselves.
inline constexpr std::uint64_t kParallelPairThreshold = 250000;

// Number of differing characters, comparing only up to the shorter sequence.
int hamming_distance(std::string_view a, std::string_view b);

// Partition into `count` contiguous chunks whose lengths differ by at most one.
std::vector<std::string_view> split_into_chunks(std::string_view sequence, int count);

// Collapse reads to distinct CDR3s, unioning gene calls across reads.
//
// The union matters: the same CDR3 is reported with different gene calls on
// different reads, and keeping only the first read's calls silently drops
// edges. Insertion order of first appearance is preserved.
std::vector<UniqueCdr3> collapse_to_unique_cdr3(const std::vector<Read>& reads);

// Edges within one block, found without enumerating all of its pairs.
//
// `block_v` and `block_j` identify the block. A pair whose records share
// several gene calls belongs to several blocks; it is emitted only by the one
// with the smallest key, so the edges returned across all blocks are already
// distinct and need no deduplication pass.
std::vector<Edge> scan_block(const std::vector<NodeId>& members,
                             const std::vector<UniqueCdr3>& nodes,
                             int cdr_length,
                             double threshold_ratio,
                             GeneId block_v,
                             GeneId block_j);

// All edges of the Hamming graph. `workers` <= 0 chooses automatically.
std::vector<Edge> build_edges(const std::vector<UniqueCdr3>& nodes,
                              double threshold_ratio,
                              int workers);

// Every pair compared directly, with no blocking, no chunk filter and no
// canonical-block rule.
//
// This is the correctness oracle, and it is **not dead code**. The fast path
// is three optimisations deep -- any one of them could silently drop an edge,
// and a missing edge merges or splits a lineage without any other symptom. The
// golden files catch a regression against a frozen result; this catches the
// algorithm itself being wrong. Quadratic by design, so it is only run over
// small inputs and the reference dataset, from the test suite.
std::vector<Edge> build_edges_reference(const std::vector<UniqueCdr3>& nodes,
                                        double threshold_ratio);

}  // namespace ar
