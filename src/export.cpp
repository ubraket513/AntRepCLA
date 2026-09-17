#include "export.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace ar {

std::string group_digits(long long value) {
    const bool negative = value < 0;
    std::string digits = std::to_string(negative ? -value : value);

    std::string grouped;
    grouped.reserve(digits.size() + digits.size() / 3 + 1);
    const std::size_t leading = digits.size() % 3 == 0 ? 3 : digits.size() % 3;
    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i == leading || (i > leading && (i - leading) % 3 == 0)) grouped.push_back(',');
        grouped.push_back(digits[i]);
    }
    return negative ? "-" + grouped : grouped;
}

std::string format_stats(const Stats& stats) {
    std::string out;
    for (const auto& [metric, value] : stats) {
        std::string number = group_digits(value);
        if (number.size() < 8) number.insert(0, 8 - number.size(), ' ');
        out += number + "  " + metric + "\n";
    }
    return out;
}

std::string format_lineage_table(const std::vector<Read>& reads,
                                 const std::vector<UniqueCdr3>& nodes,
                                 const std::vector<Lineage>& lineages) {
    ankerl::unordered_dense::map<std::string_view, long long> reads_per_cdr3;
    for (const Read& read : reads) reads_per_cdr3[std::string_view(read.cdr3)] += 1;

    std::string out = kLineageHeader;
    out.push_back('\n');

    for (std::size_t index = 0; index < lineages.size(); ++index) {
        const Lineage& lineage = lineages[index];
        long long n_reads = 0;
        for (const NodeId id : lineage) {
            n_reads += reads_per_cdr3[std::string_view(nodes[id].cdr3)];
        }

        out += std::to_string(index);
        out.push_back('\t');
        out += std::to_string(lineage.size());
        out.push_back('\t');
        out += std::to_string(n_reads);
        out.push_back('\t');
        // Members arrive sorted by sequence from get_clonal_lineages.
        for (std::size_t i = 0; i < lineage.size(); ++i) {
            if (i > 0) out.push_back(',');
            out += nodes[lineage[i]].cdr3;
        }
        out.push_back('\n');
    }
    return out;
}

std::string format_usage_table(std::vector<std::pair<std::string, long long>> usage) {
    // Sorted by gene name rather than left in discovery order: hash-map
    // ordering is an implementation detail and would not survive a port.
    std::sort(usage.begin(), usage.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    std::string out = kUsageHeader;
    out.push_back('\n');
    for (const auto& [gene, count] : usage) {
        out += gene;
        out.push_back('\t');
        out += std::to_string(count);
        out.push_back('\n');
    }
    return out;
}

void write_file(const std::string& path, const std::string& content) {
    const std::filesystem::path target(path);
    if (target.has_parent_path() && !target.parent_path().empty()) {
        std::filesystem::create_directories(target.parent_path());
    }
    std::ofstream out(target, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write to '" + path + "'.");
    out << content;
    if (!out) throw std::runtime_error("Failed while writing '" + path + "'.");
}

}  // namespace ar
