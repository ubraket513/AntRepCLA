// Portable text dumps of the pipeline's results.
//
// These formats are a contract, inherited from the Python implementation this
// replaced: tab-separated, sorted, ASCII, trailing newline, so the output can
// be diffed byte for byte against the frozen results in tests/golden/. The
// five summary statistics are too coarse to prove equivalence on their own --
// two different lineage partitions can produce identical counts -- so the
// lineage table is what actually does it. Changing any of this changes the
// contract, and there is no longer a second implementation to re-derive the
// golden files from.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "lineage.hpp"
#include "types.hpp"

namespace ar {

inline constexpr const char* kLineageHeader = "lineage\tn_unique_cdr3\tn_reads\tcdr3s";
inline constexpr const char* kUsageHeader = "v_gene\treads";

// Digits grouped in threes, as Python's ',' format specifier renders them.
std::string group_digits(long long value);

// The report written to stdout: value right-aligned in eight columns, two
// spaces, metric name.
std::string format_stats(const Stats& stats);

std::string format_lineage_table(const std::vector<Read>& reads,
                                 const std::vector<UniqueCdr3>& nodes,
                                 const std::vector<Lineage>& lineages);

std::string format_usage_table(std::vector<std::pair<std::string, long long>> usage);

// Write `content` to `path`, creating parent directories. Throws on failure.
void write_file(const std::string& path, const std::string& content);

}  // namespace ar
