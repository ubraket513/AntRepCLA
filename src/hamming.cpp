#include "hamming.hpp"

#include <algorithm>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "log.hpp"

namespace ar {
namespace {

// Canonical key for an unordered pair, so a pair found twice is recognised.
inline std::uint64_t pair_key(NodeId a, NodeId b) {
    const NodeId lo = a < b ? a : b;
    const NodeId hi = a < b ? b : a;
    return (static_cast<std::uint64_t>(lo) << 32) | hi;
}

// Blocking key: (cdr_length, v_gene, j_gene) packed into one integer. Lengths
// are small and gene tables hold tens of entries, so the ranges are safe.
inline std::uint64_t block_key(int cdr_length, GeneId v, GeneId j) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cdr_length)) << 48) |
           (static_cast<std::uint64_t>(v) << 24) |
           static_cast<std::uint64_t>(j);
}

// Smallest element common to two sorted, deduplicated lists. Both nodes of a
// pair belong to the block being scanned, so a common element always exists.
inline bool first_common(const std::vector<GeneId>& a, const std::vector<GeneId>& b,
                         GeneId& out) {
    std::size_t i = 0, k = 0;
    while (i < a.size() && k < b.size()) {
        if (a[i] == b[k]) { out = a[i]; return true; }
        if (a[i] < b[k]) ++i; else ++k;
    }
    return false;
}

// Whether this block is the one designated to emit the pair.
//
// The blocks two nodes share are the full cross product of their common V and
// J genes, all at the same CDR3 length. The smallest key in that product is
// therefore (smallest common V, smallest common J), so each pair has exactly
// one canonical block and every other block can skip it.
inline bool is_canonical_block(const UniqueCdr3& a, const UniqueCdr3& b,
                               GeneId block_v, GeneId block_j) {
    GeneId smallest_v = 0, smallest_j = 0;
    if (!first_common(a.v_call, b.v_call, smallest_v)) return false;
    if (!first_common(a.j_call, b.j_call, smallest_j)) return false;
    return block_v == smallest_v && block_j == smallest_j;
}

struct Block {
    int cdr_length = 0;
    GeneId v = 0;
    GeneId j = 0;
    std::vector<NodeId> members;
};

std::vector<Block> build_blocks(const std::vector<UniqueCdr3>& nodes) {
    ankerl::unordered_dense::map<std::uint64_t, std::size_t> index;
    std::vector<Block> blocks;

    for (NodeId id = 0; id < nodes.size(); ++id) {
        const UniqueCdr3& node = nodes[id];
        // A record with several V or J calls is indexed under each
        // combination, which is what keeps the scan exact, not approximate.
        for (const GeneId v : node.v_call) {
            for (const GeneId j : node.j_call) {
                const std::uint64_t key = block_key(node.cdr_length, v, j);
                auto found = index.find(key);
                if (found == index.end()) {
                    index.emplace(key, blocks.size());
                    blocks.push_back(Block{node.cdr_length, v, j, {id}});
                } else {
                    blocks[found->second].members.push_back(id);
                }
            }
        }
    }
    return blocks;
}

int resolve_workers(int workers, const std::vector<Block>& blocks) {
    if (workers > 0) return workers;

    std::uint64_t estimated = 0;
    for (const Block& block : blocks) {
        const std::uint64_t n = block.members.size();
        estimated += n * (n - 1) / 2;
    }
    if (estimated < kParallelPairThreshold) return 1;
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

}  // namespace

int hamming_distance(std::string_view a, std::string_view b) {
    const std::size_t length = std::min(a.size(), b.size());
    int distance = 0;
    for (std::size_t i = 0; i < length; ++i) {
        distance += (a[i] != b[i]);
    }
    return distance;
}

std::vector<std::string_view> split_into_chunks(std::string_view sequence, int count) {
    std::vector<std::string_view> chunks;
    chunks.reserve(static_cast<std::size_t>(count));

    const std::size_t base = sequence.size() / static_cast<std::size_t>(count);
    const std::size_t remainder = sequence.size() % static_cast<std::size_t>(count);

    std::size_t start = 0;
    for (int i = 0; i < count; ++i) {
        const std::size_t size = base + (static_cast<std::size_t>(i) < remainder ? 1 : 0);
        chunks.push_back(sequence.substr(start, size));
        start += size;
    }
    return chunks;
}

std::vector<UniqueCdr3> collapse_to_unique_cdr3(const std::vector<Read>& reads) {
    ankerl::unordered_dense::map<std::string, std::size_t> index;
    std::vector<UniqueCdr3> unique;
    unique.reserve(reads.size() / 4 + 1);

    const auto merge = [](std::vector<GeneId>& into, const std::vector<GeneId>& from) {
        into.insert(into.end(), from.begin(), from.end());
        std::sort(into.begin(), into.end());
        into.erase(std::unique(into.begin(), into.end()), into.end());
    };

    for (const Read& read : reads) {
        auto found = index.find(read.cdr3);
        if (found == index.end()) {
            index.emplace(read.cdr3, unique.size());
            unique.push_back(UniqueCdr3{read.cdr3, read.cdr_length, read.v_call, read.j_call});
        } else {
            UniqueCdr3& entry = unique[found->second];
            merge(entry.v_call, read.v_call);
            merge(entry.j_call, read.j_call);
        }
    }
    return unique;
}

std::vector<Edge> scan_block(const std::vector<NodeId>& members,
                             const std::vector<UniqueCdr3>& nodes,
                             int cdr_length,
                             double threshold_ratio,
                             GeneId block_v,
                             GeneId block_j) {
    std::vector<Edge> edges;
    if (members.size() < 2) return edges;

    const double max_distance = threshold_ratio * static_cast<double>(cdr_length);
    const int num_chunks = static_cast<int>(max_distance) + 1;

    // Chunked once per member, not once per position: recomputing the split
    // inside the position loop makes the scan quadratic in the chunk count.
    // Each member is split by its own length, matching the reference
    // implementation -- which is what keeps the two identical if a block ever
    // holds sequences whose reported length and actual length disagree.
    std::vector<std::vector<std::string_view>> chunked;
    chunked.reserve(members.size());
    for (const NodeId id : members) {
        chunked.push_back(split_into_chunks(nodes[id].cdr3, num_chunks));
    }

    // Only accepted edges are retained; a rejected pair colliding on several
    // chunks is measured more than once, which is the bounded price of never
    // materialising the candidate set.
    ankerl::unordered_dense::set<std::uint64_t> accepted;
    ankerl::unordered_dense::map<std::string_view, std::vector<NodeId>> by_chunk;

    // One chunk position at a time, so the index key stays a plain string_view
    // and the map can be reused rather than keyed on (position, chunk).
    for (int position = 0; position < num_chunks; ++position) {
        by_chunk.clear();
        for (std::size_t m = 0; m < members.size(); ++m) {
            by_chunk[chunked[m][static_cast<std::size_t>(position)]].push_back(members[m]);
        }

        for (const auto& [chunk, group] : by_chunk) {
            if (group.size() < 2) continue;
            for (std::size_t i = 0; i < group.size(); ++i) {
                for (std::size_t k = i + 1; k < group.size(); ++k) {
                    // A pair with several gene calls in common sits in several
                    // blocks at once. Emitting it from only its canonical
                    // block removes the duplicates at the source, rather than
                    // sorting them out of the finished edge list -- which was
                    // a serial pass over every edge in the graph. Checked
                    // before the distance, so the redundant measurements go
                    // away too.
                    if (!is_canonical_block(nodes[group[i]], nodes[group[k]], block_v, block_j)) {
                        continue;
                    }
                    const std::uint64_t key = pair_key(group[i], group[k]);
                    if (accepted.contains(key)) continue;
                    if (hamming_distance(nodes[group[i]].cdr3, nodes[group[k]].cdr3) <= max_distance) {
                        accepted.insert(key);
                        edges.emplace_back(group[i], group[k]);
                    }
                }
            }
        }
    }
    return edges;
}

std::vector<Edge> build_edges_reference(const std::vector<UniqueCdr3>& nodes,
                                        double threshold_ratio) {
    const auto shares_gene = [](const std::vector<GeneId>& a, const std::vector<GeneId>& b) {
        GeneId unused = 0;
        return first_common(a, b, unused);
    };

    std::vector<Edge> edges;
    for (NodeId a = 0; a < nodes.size(); ++a) {
        for (NodeId b = a + 1; b < nodes.size(); ++b) {
            if (nodes[a].cdr_length != nodes[b].cdr_length) continue;
            if (!shares_gene(nodes[a].v_call, nodes[b].v_call)) continue;
            if (!shares_gene(nodes[a].j_call, nodes[b].j_call)) continue;
            const double limit = threshold_ratio * static_cast<double>(nodes[a].cdr_length);
            if (hamming_distance(nodes[a].cdr3, nodes[b].cdr3) > limit) continue;
            edges.emplace_back(a, b);
        }
    }
    return edges;
}

std::vector<Edge> build_edges(const std::vector<UniqueCdr3>& nodes,
                              double threshold_ratio,
                              int workers) {
    std::vector<Block> blocks = build_blocks(nodes);
    const int resolved = resolve_workers(workers, blocks);

    // Largest blocks first, so one straggler does not define the wall time.
    std::sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) {
        return a.members.size() > b.members.size();
    });

    std::vector<Edge> edges;
#ifdef _OPENMP
    if (resolved > 1) {
        std::vector<std::vector<Edge>> per_thread(static_cast<std::size_t>(resolved));
        const long count = static_cast<long>(blocks.size());
#pragma omp parallel for schedule(dynamic) num_threads(resolved)
        for (long i = 0; i < count; ++i) {
            const int thread = omp_get_thread_num();
            std::vector<Edge> found = scan_block(blocks[i].members, nodes,
                                                 blocks[i].cdr_length, threshold_ratio,
                                                 blocks[i].v, blocks[i].j);
            auto& sink = per_thread[static_cast<std::size_t>(thread)];
            sink.insert(sink.end(), found.begin(), found.end());
        }
        for (const auto& chunk : per_thread) {
            edges.insert(edges.end(), chunk.begin(), chunk.end());
        }
    } else
#endif
    {
        for (const Block& block : blocks) {
            std::vector<Edge> found = scan_block(block.members, nodes, block.cdr_length,
                                                 threshold_ratio, block.v, block.j);
            edges.insert(edges.end(), found.begin(), found.end());
        }
    }

    // No global deduplication pass: each pair is emitted by its canonical
    // block alone, so the edge list is already distinct. That pass used to be
    // a serial sort over every edge in the graph and capped the speedup at
    // roughly 2x however many cores were available.
    log_info("Hamming graph: " + std::to_string(nodes.size()) + " nodes, " +
             std::to_string(edges.size()) + " edges, " + std::to_string(blocks.size()) +
             " blocks scanned on " + std::to_string(resolved) + " worker(s).");
    return edges;
}

}  // namespace ar
