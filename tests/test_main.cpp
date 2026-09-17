// Unit tests, carried over from the Python suite this replaced.
//
// These cover the pieces in isolation. Whole-pipeline equivalence is proved
// separately, by diffing this binary's output against the frozen golden files
// in tests/golden/ (see tools/verify_output.sh).
#include "doctest.h"

#include <algorithm>
#include <filesystem>
#include <random>
#include <tuple>
#include <string>
#include <vector>

#include "export.hpp"
#include "hamming.hpp"
#include "igblast.hpp"
#include "lineage.hpp"
#include "log.hpp"
#include "svg.hpp"

using namespace ar;

namespace {

// Build reads from (cdr3, v genes, j genes) triples, mirroring make_df().
std::vector<Read> make_reads(
    const std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>>& rows,
    GeneTable& genes) {
    std::vector<Read> reads;
    for (const auto& [cdr3, v, j] : rows) {
        Read read;
        read.cdr3 = cdr3;
        read.cdr3_aa = std::string(cdr3.size() / 3, 'X');
        read.cdr_length = static_cast<int>(cdr3.size());
        for (const std::string& name : v) read.v_call.push_back(genes.intern(name));
        for (const std::string& name : j) read.j_call.push_back(genes.intern(name));
        std::sort(read.v_call.begin(), read.v_call.end());
        std::sort(read.j_call.begin(), read.j_call.end());
        reads.push_back(std::move(read));
    }
    return reads;
}

std::size_t count_edges(const std::vector<Read>& reads, double ratio = 0.1) {
    const auto nodes = collapse_to_unique_cdr3(reads);
    return build_edges(nodes, ratio, 1).size();
}

// Both edge lists as comparable sets: the fast path emits per block, in
// whatever order the threads finished.
std::vector<Edge> as_set(std::vector<Edge> edges) {
    for (Edge& edge : edges) {
        if (edge.second < edge.first) std::swap(edge.first, edge.second);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

std::string mutate(std::string sequence, const std::vector<std::size_t>& positions) {
    for (const std::size_t pos : positions) {
        sequence[pos] = sequence[pos] != 'A' ? 'A' : 'C';
    }
    return sequence;
}

const std::string kBase32 = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT";
const std::string kBase40 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

}  // namespace

TEST_CASE("hamming distance counts mismatches") {
    CHECK(hamming_distance("AAAA", "AAAA") == 0);
    CHECK(hamming_distance("AAAA", "AAAT") == 1);
    CHECK(hamming_distance("AAAA", "TTTT") == 4);
}

TEST_CASE("chunks partition the sequence") {
    const auto chunks = split_into_chunks("AAAACCCCGGGG", 5);
    CHECK(chunks.size() == 5);
    std::string joined;
    for (const auto& chunk : chunks) joined += std::string(chunk);
    CHECK(joined == "AAAACCCCGGGG");
}

TEST_CASE("chunks are balanced when the length is indivisible") {
    const auto chunks = split_into_chunks("AAAAAAAAAA", 4);
    std::size_t smallest = 99, largest = 0, total = 0;
    for (const auto& chunk : chunks) {
        smallest = std::min(smallest, chunk.size());
        largest = std::max(largest, chunk.size());
        total += chunk.size();
    }
    CHECK(largest - smallest <= 1);
    CHECK(total == 10);
}

TEST_CASE("a chunk count exceeding the length still partitions") {
    const auto chunks = split_into_chunks("ACGT", 7);
    CHECK(chunks.size() == 7);
    std::string joined;
    for (const auto& chunk : chunks) joined += std::string(chunk);
    CHECK(joined == "ACGT");
}

TEST_CASE("identical sequences form one node with no self-loop") {
    GeneTable genes;
    const auto reads = make_reads({{kBase32, {"V1"}, {"J1"}}, {kBase32, {"V1"}, {"J1"}}}, genes);
    CHECK(collapse_to_unique_cdr3(reads).size() == 1);
    CHECK(count_edges(reads) == 0);
}

TEST_CASE("close sequences are joined") {
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase32, {"V1"}, {"J1"}}, {mutate(kBase32, {31}), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads) == 1);
}

TEST_CASE("distant sequences are not joined") {
    GeneTable genes;
    const auto reads = make_reads({{std::string(32, 'A'), {"V1"}, {"J1"}},
                                   {std::string(32, 'T'), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads) == 0);
}

TEST_CASE("the threshold ratio is honoured") {
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase32, {"V1"}, {"J1"}}, {mutate(kBase32, {30, 31}), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads, 0.0) == 0);
    CHECK(count_edges(reads, 0.5) == 1);
}

TEST_CASE("different lengths are never joined") {
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase32, {"V1"}, {"J1"}}, {kBase32.substr(0, 31), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads) == 0);
}

TEST_CASE("a gene mismatch blocks an edge") {
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase32, {"V1"}, {"J1"}}, {mutate(kBase32, {31}), {"V2"}, {"J2"}}}, genes);
    CHECK(count_edges(reads) == 0);
}

TEST_CASE("partial gene overlap allows an edge") {
    GeneTable genes;
    const auto reads = make_reads({{kBase32, {"V1", "V9"}, {"J1", "J9"}},
                                   {mutate(kBase32, {31}), {"V9"}, {"J9"}}}, genes);
    CHECK(count_edges(reads) == 1);
}

TEST_CASE("a pair sharing several blocks yields exactly one edge") {
    // Both records carry two V and two J calls, so the pair sits in four
    // blocking buckets at once. It must still produce a single edge -- which
    // is what lets the scan emit each pair from one designated block rather
    // than deduplicating the whole edge list afterwards.
    GeneTable genes;
    const auto reads = make_reads({{kBase32, {"V1", "V9"}, {"J1", "J9"}},
                                   {mutate(kBase32, {31}), {"V1", "V9"}, {"J1", "J9"}}}, genes);
    CHECK(count_edges(reads) == 1);
}

TEST_CASE("three mutually close records sharing many blocks give three edges") {
    GeneTable genes;
    const auto reads = make_reads({{kBase32, {"V1", "V9"}, {"J1", "J9"}},
                                   {mutate(kBase32, {31}), {"V1", "V9"}, {"J1", "J9"}},
                                   {mutate(kBase32, {30}), {"V1", "V9"}, {"J1", "J9"}}}, genes);
    CHECK(count_edges(reads) == 3);
}

TEST_CASE("gene calls are unioned across duplicate reads") {
    // The shared CDR3 is reported with V1 on one read and V9 on another.
    // Keeping only the first would hide the V9 route and lose the edge.
    GeneTable genes;
    const auto reads = make_reads({{kBase32, {"V1"}, {"J1"}},
                                   {kBase32, {"V9"}, {"J9"}},
                                   {mutate(kBase32, {31}), {"V9"}, {"J9"}}}, genes);
    const auto nodes = collapse_to_unique_cdr3(reads);
    CHECK(nodes.size() == 2);
    CHECK(build_edges(nodes, 0.1, 1).size() == 1);
}

TEST_CASE("mismatches spread across every chunk but one are still joined") {
    // The adversarial case for pigeonhole: 40 bp admits 4 mismatches and splits
    // into 5 chunks. Four mismatches land in four different chunks; the fifth
    // chunk survives intact and must still find the pair.
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase40, {"V1"}, {"J1"}}, {mutate(kBase40, {0, 8, 16, 24}), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads) == 1);
}

TEST_CASE("mismatches in every chunk are not joined") {
    GeneTable genes;
    const auto reads = make_reads(
        {{kBase40, {"V1"}, {"J1"}}, {mutate(kBase40, {0, 8, 16, 24, 32}), {"V1"}, {"J1"}}}, genes);
    CHECK(count_edges(reads) == 0);
}

TEST_CASE("a parallel scan finds the same edges as a serial one") {
    GeneTable genes;
    std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>> rows;
    for (int block = 0; block < 6; ++block) {
        std::string base = kBase40;
        for (std::size_t i = 0; i < base.size(); ++i) {
            base[i] = "ACGT"[(block * 7 + static_cast<int>(i)) % 4];
        }
        for (std::size_t i = 0; i < 12; ++i) {
            rows.push_back({mutate(base, {i}), {"V" + std::to_string(block)},
                            {"J" + std::to_string(block)}});
        }
    }
    const auto reads = make_reads(rows, genes);
    const auto nodes = collapse_to_unique_cdr3(reads);

    // Compared as sets. Edges are no longer emitted in a globally sorted
    // order -- each block emits its own and they are concatenated -- so the
    // order depends on thread scheduling while the set does not. Determinism
    // of the *output* comes from the canonical lineage ordering, which
    // tools/verify_output.sh checks end to end.
    const auto normalise = [](std::vector<Edge> edges) {
        for (Edge& edge : edges) {
            if (edge.second < edge.first) std::swap(edge.first, edge.second);
        }
        std::sort(edges.begin(), edges.end());
        return edges;
    };

    const auto serial = normalise(build_edges(nodes, 0.1, 1));
    for (int workers : {2, 3, 4}) {
        CHECK(normalise(build_edges(nodes, 0.1, workers)) == serial);
    }
}

TEST_CASE("the edge list contains no duplicates without a deduplication pass") {
    // Every record shares two V and two J calls with every other, so each pair
    // sits in four blocks. Only the canonical block may emit it.
    GeneTable genes;
    std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>> rows;
    for (std::size_t i = 0; i < 8; ++i) {
        rows.push_back({mutate(kBase32, {i}), {"V1", "V9"}, {"J1", "J9"}});
    }
    const auto reads = make_reads(rows, genes);
    const auto nodes = collapse_to_unique_cdr3(reads);

    for (int workers : {1, 4}) {
        auto edges = build_edges(nodes, 0.1, workers);
        const std::size_t emitted = edges.size();
        for (Edge& edge : edges) {
            if (edge.second < edge.first) std::swap(edge.first, edge.second);
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        CHECK(edges.size() == emitted);
    }
}

TEST_CASE("allele suffixes are stripped everywhere they appear") {
    CHECK(strip_alleles("IGHV3-30*14") == "IGHV3-30");
    CHECK(strip_alleles("IGHV3-30*14,IGHV3-30-3*01") == "IGHV3-30,IGHV3-30-3");
    CHECK(strip_alleles("IGHV3-30") == "IGHV3-30");
    // A '*' not followed by digits is not an allele suffix and stays put.
    CHECK(strip_alleles("IGHV*X") == "IGHV*X");
}

TEST_CASE("gene calls parse to sorted unique ids") {
    GeneTable genes;
    const auto ids = parse_gene_call("IGHV3-30*14,IGHV3-30*02,IGHV1-2*01", genes);
    // The first two collapse to the same gene once alleles are stripped.
    CHECK(ids.size() == 2);
    CHECK(std::is_sorted(ids.begin(), ids.end()));
    CHECK(parse_gene_call("", genes).empty());
}

TEST_CASE("lineages of equal size are ordered by smallest member") {
    std::vector<UniqueCdr3> nodes{{"CCCC", 4, {}, {}}, {"AAAA", 4, {}, {}}, {"BBBB", 4, {}, {}}};
    const auto lineages = get_clonal_lineages(nodes, {});
    REQUIRE(lineages.size() == 3);
    CHECK(nodes[lineages[0].front()].cdr3 == "AAAA");
    CHECK(nodes[lineages[1].front()].cdr3 == "BBBB");
    CHECK(nodes[lineages[2].front()].cdr3 == "CCCC");
}

TEST_CASE("larger lineages still come first") {
    std::vector<UniqueCdr3> nodes{
        {"ZZZZ", 4, {}, {}}, {"ZZZA", 4, {}, {}}, {"ZZAA", 4, {}, {}}, {"AAAA", 4, {}, {}}};
    const auto lineages = get_clonal_lineages(nodes, {{0, 1}, {0, 2}});
    REQUIRE(lineages.size() == 2);
    CHECK(lineages[0].size() == 3);
    CHECK(nodes[lineages[1].front()].cdr3 == "AAAA");
}

TEST_CASE("digits are grouped like Python's comma format") {
    CHECK(group_digits(0) == "0");
    CHECK(group_digits(623) == "623");
    CHECK(group_digits(1000) == "1,000");
    CHECK(group_digits(12345) == "12,345");
    CHECK(group_digits(1234567) == "1,234,567");
    CHECK(group_digits(-1234) == "-1,234");
}

TEST_CASE("the stats report is right-aligned in eight columns") {
    CHECK(format_stats({{"number of clonal lineages", 623}}) ==
          "     623  number of clonal lineages\n");
}

TEST_CASE("the usage table is sorted by gene name") {
    CHECK(format_usage_table({{"IGHV3", 5}, {"IGHV1", 2}, {"IGHV2", 9}}) ==
          "v_gene\treads\nIGHV1\t2\nIGHV2\t9\nIGHV3\t5\n");
}

TEST_CASE("the lineage table counts reads, not unique CDR3s") {
    GeneTable genes;
    const auto reads = make_reads(
        {{"AAAA", {"V1"}, {"J1"}}, {"AAAA", {"V1"}, {"J1"}}, {"AAAA", {"V1"}, {"J1"}}}, genes);
    const auto nodes = collapse_to_unique_cdr3(reads);
    const auto lineages = get_clonal_lineages(nodes, {});
    CHECK(format_lineage_table(reads, nodes, lineages) ==
          "lineage\tn_unique_cdr3\tn_reads\tcdr3s\n0\t1\t3\tAAAA\n");
}

TEST_CASE("log lines carry timestamp, padded level and message") {
    // The layout matches what the reference pipeline emits, so the two can be
    // read side by side: timestamp, level padded to eight columns, message.
    CHECK(format_log_line("2026-09-17 15:34:47", LogLevel::info, "Analysis complete.") ==
          "2026-09-17 15:34:47 | info     | Analysis complete.");
    CHECK(format_log_line("2026-09-17 15:34:47", LogLevel::warning, "careful") ==
          "2026-09-17 15:34:47 | warning  | careful");
    CHECK(format_log_line("2026-09-17 15:34:47", LogLevel::error, "boom") ==
          "2026-09-17 15:34:47 | error    | boom");
}

TEST_CASE("messages below the configured level are suppressed") {
    CHECK(should_log(LogLevel::info, LogLevel::info));
    CHECK(should_log(LogLevel::info, LogLevel::warning));
    CHECK(should_log(LogLevel::info, LogLevel::error));
    CHECK_FALSE(should_log(LogLevel::info, LogLevel::debug));
    CHECK_FALSE(should_log(LogLevel::error, LogLevel::warning));
    CHECK(should_log(LogLevel::debug, LogLevel::debug));
}

TEST_CASE("log levels parse from their names") {
    CHECK(parse_log_level("DEBUG") == LogLevel::debug);
    CHECK(parse_log_level("INFO") == LogLevel::info);
    CHECK(parse_log_level("WARNING") == LogLevel::warning);
    CHECK(parse_log_level("ERROR") == LogLevel::error);
    CHECK_THROWS_AS(parse_log_level("LOUD"), std::invalid_argument);
}


// --- The blocked scan against brute force -----------------------------------
//
// Blocking, the pigeonhole chunk filter and the canonical-block rule are three
// independent optimisations. Each one can only ever *drop* an edge, and a
// dropped edge silently splits a lineage. These compare against the definition.

TEST_CASE("the fast scan agrees with brute force on generated repertoires") {
    // Families of near-neighbours over a shared gene pool, so blocks overlap,
    // pairs land in several blocks at once, and mutations spread across chunk
    // boundaries.
    std::mt19937 rng(20260917);
    for (int trial = 0; trial < 12; ++trial) {
        GeneTable genes;
        std::vector<std::tuple<std::string, std::vector<std::string>, std::vector<std::string>>> rows;

        const int families = 6;
        for (int f = 0; f < families; ++f) {
            std::string base;
            const std::size_t length = 30 + (rng() % 16);   // varying CDR3 lengths
            for (std::size_t i = 0; i < length; ++i) base += "ACGT"[rng() % 4];

            for (int k = 0; k < 6; ++k) {
                std::string variant = base;
                const int mutations = static_cast<int>(rng() % 6);
                for (int m = 0; m < mutations; ++m) {
                    variant[rng() % variant.size()] = "ACGT"[rng() % 4];
                }
                // Overlapping, sometimes multiple, gene calls.
                std::vector<std::string> v{"V" + std::to_string(rng() % 4)};
                std::vector<std::string> j{"J" + std::to_string(rng() % 3)};
                if (rng() % 3 == 0) v.push_back("V" + std::to_string(rng() % 4));
                if (rng() % 3 == 0) j.push_back("J" + std::to_string(rng() % 3));
                std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
                std::sort(j.begin(), j.end()); j.erase(std::unique(j.begin(), j.end()), j.end());
                rows.push_back({variant, v, j});
            }
        }

        const auto reads = make_reads(rows, genes);
        const auto nodes = collapse_to_unique_cdr3(reads);
        const auto expected = as_set(build_edges_reference(nodes, 0.1));
        CHECK(as_set(build_edges(nodes, 0.1, 1)) == expected);
        CHECK(as_set(build_edges(nodes, 0.1, 4)) == expected);
    }
}

TEST_CASE("the fast scan agrees with brute force on the real dataset") {
    const std::string path = "data/igblast_results.tsv";
    if (!std::filesystem::is_regular_file(path)) {
        MESSAGE("real dataset not present; skipping");
        return;
    }
    const auto repertoire = load_igblast_results(path);
    const auto nodes = collapse_to_unique_cdr3(repertoire.reads);

    const auto expected = as_set(build_edges_reference(nodes, 0.1));
    CHECK(expected.size() == 11294);        // pinned: the graph the pipeline has always produced
    CHECK(as_set(build_edges(nodes, 0.1, 1)) == expected);
    CHECK(as_set(build_edges(nodes, 0.1, 8)) == expected);
}


// --- SVG plotting ------------------------------------------------------------

namespace {
bool contains(const std::string& haystack, const std::string& needle) {
    return haystack.find(needle) != std::string::npos;
}
std::size_t count_of(const std::string& haystack, const std::string& needle) {
    std::size_t n = 0, pos = 0;
    while ((pos = haystack.find(needle, pos)) != std::string::npos) { ++n; pos += needle.size(); }
    return n;
}
}  // namespace

TEST_CASE("a vertical bar chart renders one rect per bar") {
    const std::string svg = render_vertical_bars(
        "V Gene Usage Statistics", "V Genes", "Read Count",
        {{"IGHV1-2", 719}, {"IGHV1-18", 249}, {"IGHV1-24", 208}});

    CHECK(contains(svg, "<svg"));
    CHECK(contains(svg, "</svg>"));
    CHECK(contains(svg, "V Gene Usage Statistics"));
    CHECK(contains(svg, "IGHV1-2"));
    CHECK(count_of(svg, "<rect class=\"bar\"") == 3);
}

TEST_CASE("a horizontal bar chart labels each bar with its value") {
    const std::string svg = render_horizontal_bars(
        "Clonal Lineage Statistics", "Value",
        {{"number of clonal lineages", 623}, {"largest lineage", 103}});

    CHECK(count_of(svg, "<rect class=\"bar\"") == 2);
    CHECK(contains(svg, "623"));
    CHECK(contains(svg, "number of clonal lineages"));
}

TEST_CASE("plot text is XML escaped") {
    // Gene names are tame, but a label reaching the renderer unescaped would
    // produce a file no SVG reader will open.
    const std::string svg = render_vertical_bars("t", "x", "y", {{"a<b&c\"d", 1}});
    CHECK(contains(svg, "a&lt;b&amp;c&quot;d"));
    CHECK_FALSE(contains(svg, "a<b&c"));
}

TEST_CASE("an empty series is rejected rather than rendered blank") {
    CHECK_THROWS_AS(render_vertical_bars("t", "x", "y", {}), std::invalid_argument);
    CHECK_THROWS_AS(render_horizontal_bars("t", "x", {}), std::invalid_argument);
}

TEST_CASE("bars are scaled to the largest value, and zero is handled") {
    // All-zero counts must not divide by zero; every bar simply has no height.
    const std::string svg = render_vertical_bars("t", "x", "y", {{"a", 0}, {"b", 0}});
    CHECK(count_of(svg, "<rect class=\"bar\"") == 2);
    CHECK_FALSE(contains(svg, "nan"));
    CHECK_FALSE(contains(svg, "inf"));
}
