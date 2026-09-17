// Clonal lineage extraction and statistics.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "hamming.hpp"
#include "types.hpp"

namespace ar {

// A lineage is a connected component: the indices of its member CDR3s, sorted
// by sequence.
using Lineage = std::vector<NodeId>;

// Connected components in canonical order: largest first, ties broken by the
// lexicographically smallest member.
//
// The tie-break is load-bearing. Ordering on size alone leaves equal-sized
// lineages in whatever order the components were enumerated, which depends on
// insertion order and on the graph implementation. Lineage 0 is the one the
// WebLogo query is built from, so that ordering has to be total.
std::vector<Lineage> get_clonal_lineages(const std::vector<UniqueCdr3>& nodes,
                                         const std::vector<Edge>& edges);

// Metric name and value, in report order.
using Stats = std::vector<std::pair<std::string, long long>>;

Stats get_lineage_stats(const std::vector<Read>& reads,
                        const std::vector<UniqueCdr3>& nodes,
                        const std::vector<Lineage>& lineages,
                        int min_lineage_size);

// Reads per V gene. Every V gene seen in the input appears, including those
// with a count of zero.
std::vector<std::pair<std::string, long long>> get_v_gene_usage(
    const std::vector<Read>& reads, const GeneTable& genes);

// WebLogo query built from the largest lineage: its CDR3 amino acid sequences
// in input read order, spaces replaced by '*', padded to the longest.
std::string get_cdr3_aa_from_largest_lineage(const std::vector<Read>& reads,
                                             const std::vector<UniqueCdr3>& nodes,
                                             const std::vector<Lineage>& lineages);

}  // namespace ar
