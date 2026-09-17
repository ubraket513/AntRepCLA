#include "igblast.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

#include "csv.hpp"
#include "log.hpp"

namespace ar {

std::string strip_alleles(const std::string& call) {
    std::string out;
    out.reserve(call.size());
    for (std::size_t i = 0; i < call.size(); ++i) {
        if (call[i] == '*' && i + 1 < call.size() &&
            std::isdigit(static_cast<unsigned char>(call[i + 1]))) {
            // Drop the '*' and the digit run that follows it.
            std::size_t j = i + 1;
            while (j < call.size() && std::isdigit(static_cast<unsigned char>(call[j]))) ++j;
            i = j - 1;
            continue;
        }
        out.push_back(call[i]);
    }
    return out;
}

std::vector<GeneId> parse_gene_call(const std::string& field, GeneTable& genes) {
    std::vector<GeneId> ids;
    if (field.empty()) return ids;

    const std::string stripped = strip_alleles(field);
    std::size_t start = 0;
    while (start <= stripped.size()) {
        const std::size_t comma = stripped.find(',', start);
        const std::size_t end = (comma == std::string::npos) ? stripped.size() : comma;
        ids.push_back(genes.intern(stripped.substr(start, end - start)));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }

    // Sorted and unique, so that intersection tests are a linear merge and so
    // that the same call set always compares equal regardless of input order.
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

Repertoire load_igblast_results(const std::string& path, int workers) {
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("Input file '" + path + "' does not exist.");
    }
    log_info("Loading IgBLAST results from " + path);

    csv::CSVFormat format;
    format.delimiter('\t');
    // Parsing dominates the run once the scan is threaded, and csv-parser can
    // split the source across threads -- but only above a size threshold that
    // defaults far higher than these files ever reach.
    format.speculative_parallel_min_bytes(kParallelParseMinBytes);
    if (workers > 0) {
        format.speculative_parallel_threads(static_cast<std::size_t>(workers));
    }
    csv::CSVReader reader(path, format);

    // Resolve column positions once. The AIRR output has 96 columns and
    // looking each one up by name per row is a measurable cost.
    const auto column = [&reader, &path](const char* name) {
        const int index = reader.index_of(name);
        if (index < 0) {
            throw std::runtime_error("Column '" + std::string(name) + "' missing from " + path);
        }
        return static_cast<std::size_t>(index);
    };
    const std::size_t i_v = column("v_call");
    const std::size_t i_j = column("j_call");
    const std::size_t i_cdr3 = column("cdr3");
    const std::size_t i_cdr3_aa = column("cdr3_aa");
    const std::size_t i_start = column("cdr3_start");
    const std::size_t i_end = column("cdr3_end");

    Repertoire repertoire;
    ankerl::unordered_dense::set<std::string> distinct;

    for (csv::CSVRow& row : reader) {
        Read read;
        read.cdr3 = row[i_cdr3].get<std::string>();
        read.cdr3_aa = row[i_cdr3_aa].get<std::string>();

        // Length comes from the reported coordinates, not from the sequence
        // itself, matching the pipeline this is a port of.
        read.cdr_length = row[i_end].get<int>() - row[i_start].get<int>() + 1;

        read.v_call = parse_gene_call(row[i_v].is_null() ? "" : row[i_v].get<std::string>(),
                                      repertoire.genes);
        read.j_call = parse_gene_call(row[i_j].is_null() ? "" : row[i_j].get<std::string>(),
                                      repertoire.genes);

        distinct.insert(read.cdr3);
        repertoire.reads.push_back(std::move(read));
    }

    log_info("Loaded " + std::to_string(repertoire.reads.size()) + " sequences (" +
             std::to_string(distinct.size()) + " unique CDR3s).");
    return repertoire;
}

}  // namespace ar
