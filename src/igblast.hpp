// Loading and normalising NCBI IgBLAST AIRR-format (-outfmt 19) output.
#pragma once

#include <cstddef>
#include <string>

#include "types.hpp"

namespace ar {

// Strip IMGT allele suffixes: "IGHV3-30*14" -> "IGHV3-30".
//
// Matches the pipeline's `str.replace(r"\*\d+", '')`, which removes every
// occurrence, not just a trailing one.
std::string strip_alleles(const std::string& call);

// Split a comma-separated gene call field into interned, sorted, unique ids.
std::vector<GeneId> parse_gene_call(const std::string& field, GeneTable& genes);

// Minimum input size before the parser is allowed to split the file across
// threads. csv-parser defaults this to 50 MB, which never triggers on a
// repertoire of this size; below a few megabytes the split costs more than it
// saves.
inline constexpr std::size_t kParallelParseMinBytes = 2 * 1024 * 1024;

// Read the TSV. Throws std::runtime_error if the file is missing or lacks a
// required column.
//
// `workers` <= 0 lets the parser choose its own thread count.
Repertoire load_igblast_results(const std::string& path, int workers = 0);

}  // namespace ar
