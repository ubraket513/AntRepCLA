#include "lineage.hpp"

#include <algorithm>

#include "log.hpp"

namespace ar {
namespace {

// Disjoint-set forest with union by size and path halving. This replaces the
// graph library the reference implementation used: the only thing asked of the
// graph is its connected components.
class DisjointSets {
public:
    explicit DisjointSets(std::size_t count) : parent_(count), size_(count, 1) {
        for (std::size_t i = 0; i < count; ++i) parent_[i] = static_cast<NodeId>(i);
    }

    NodeId find(NodeId id) {
        while (parent_[id] != id) {
            parent_[id] = parent_[parent_[id]];
            id = parent_[id];
        }
        return id;
    }

    void unite(NodeId a, NodeId b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (size_[a] < size_[b]) std::swap(a, b);
        parent_[b] = a;
        size_[a] += size_[b];
    }

private:
    std::vector<NodeId> parent_;
    std::vector<std::uint32_t> size_;
};

}  // namespace

std::vector<Lineage> get_clonal_lineages(const std::vector<UniqueCdr3>& nodes,
                                         const std::vector<Edge>& edges) {
    DisjointSets sets(nodes.size());
    for (const Edge& edge : edges) sets.unite(edge.first, edge.second);

    // Group by representative. Isolated CDR3s are lineages of size one, so
    // every node participates.
    ankerl::unordered_dense::map<NodeId, std::size_t> position;
    std::vector<Lineage> lineages;
    for (NodeId id = 0; id < nodes.size(); ++id) {
        const NodeId root = sets.find(id);
        auto found = position.find(root);
        if (found == position.end()) {
            position.emplace(root, lineages.size());
            lineages.push_back(Lineage{id});
        } else {
            lineages[found->second].push_back(id);
        }
    }

    for (Lineage& lineage : lineages) {
        std::sort(lineage.begin(), lineage.end(), [&nodes](NodeId a, NodeId b) {
            return nodes[a].cdr3 < nodes[b].cdr3;
        });
    }

    // Largest first, ties by smallest member. Members are already sorted, so
    // the smallest is at the front.
    std::sort(lineages.begin(), lineages.end(),
              [&nodes](const Lineage& a, const Lineage& b) {
                  if (a.size() != b.size()) return a.size() > b.size();
                  return nodes[a.front()].cdr3 < nodes[b.front()].cdr3;
              });

    log_info("Identified " + std::to_string(lineages.size()) + " clonal lineages.");
    return lineages;
}

Stats get_lineage_stats(const std::vector<Read>& reads,
                        const std::vector<UniqueCdr3>& nodes,
                        const std::vector<Lineage>& lineages,
                        int min_lineage_size) {
    // Map each distinct CDR3 to its lineage index, then tally reads per
    // lineage. An integer label is used rather than the member set because
    // set-derived keys are not safe to group on.
    ankerl::unordered_dense::map<std::string_view, std::size_t> cdr3_to_lineage;
    for (std::size_t index = 0; index < lineages.size(); ++index) {
        for (const NodeId id : lineages[index]) {
            cdr3_to_lineage.emplace(std::string_view(nodes[id].cdr3), index);
        }
    }

    std::vector<long long> reads_per_lineage(lineages.size(), 0);
    for (const Read& read : reads) {
        auto found = cdr3_to_lineage.find(std::string_view(read.cdr3));
        if (found != cdr3_to_lineage.end()) reads_per_lineage[found->second] += 1;
    }

    const long long largest_unique =
        lineages.empty() ? 0 : static_cast<long long>(lineages.front().size());
    const long long largest_reads = lineages.empty() ? 0 : reads_per_lineage.front();

    long long by_unique = 0;
    for (const Lineage& lineage : lineages) {
        if (static_cast<int>(lineage.size()) >= min_lineage_size) ++by_unique;
    }
    long long by_reads = 0;
    for (const long long count : reads_per_lineage) {
        if (count >= min_lineage_size) ++by_reads;
    }

    const std::string min_text = std::to_string(min_lineage_size);
    Stats stats{
        {"number of clonal lineages", static_cast<long long>(lineages.size())},
        {"number of unique CDR3s in the largest clonal lineage", largest_unique},
        {"number of sequences in the largest clonal lineage", largest_reads},
        {"number of lineages represented by at least " + min_text + " unique CDR3s", by_unique},
        {"number of lineages represented by at least " + min_text + " sequences", by_reads},
    };

    log_info("Computed clonal lineage statistics.");
    return stats;
}

std::vector<std::pair<std::string, long long>> get_v_gene_usage(
    const std::vector<Read>& reads, const GeneTable& genes) {
    // Keyed on V genes only: a J gene that never appears as a V call must not
    // show up in the table, even though both share an intern table.
    ankerl::unordered_dense::map<GeneId, long long> counts;
    for (const Read& read : reads) {
        for (const GeneId gene : read.v_call) counts[gene] += 1;
    }

    std::vector<std::pair<std::string, long long>> usage;
    usage.reserve(counts.size());
    for (const auto& [gene, count] : counts) {
        usage.emplace_back(genes.name(gene), count);
    }

    log_info("Computed usage across " + std::to_string(usage.size()) + " V genes.");
    return usage;
}

std::string get_cdr3_aa_from_largest_lineage(const std::vector<Read>& reads,
                                             const std::vector<UniqueCdr3>& nodes,
                                             const std::vector<Lineage>& lineages) {
    if (lineages.empty()) {
        log_warning("No lineages available; returning an empty query.");
        return "";
    }

    ankerl::unordered_dense::set<std::string_view> largest;
    for (const NodeId id : lineages.front()) largest.insert(std::string_view(nodes[id].cdr3));

    // Input read order, which is what an inner join against the read table
    // preserves in the reference implementation.
    std::vector<std::string> sequences;
    for (const Read& read : reads) {
        if (!largest.contains(std::string_view(read.cdr3))) continue;
        std::string amino = read.cdr3_aa;
        std::replace(amino.begin(), amino.end(), ' ', '*');
        sequences.push_back(std::move(amino));
    }

    if (sequences.empty()) {
        log_warning("Largest lineage has no amino acid sequences.");
        return "";
    }

    // Pad to the longest, not the mean: padding to the mean leaves every
    // above-average sequence unpadded and the output ragged.
    std::size_t target = 0;
    for (const std::string& sequence : sequences) target = std::max(target, sequence.size());

    std::string query;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        if (i > 0) query.push_back('\n');
        query.append(sequences[i]);
        query.append(target - sequences[i].size(), '*');
    }

    log_info("Extracted " + std::to_string(sequences.size()) + " CDR3 amino acid sequences.");
    return query;
}

}  // namespace ar
