// Command-line entry point.
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "export.hpp"
#include "hamming.hpp"
#include "igblast.hpp"
#include "lineage.hpp"
#include "log.hpp"
#include "svg.hpp"

namespace {

struct Options {
    std::string igblast;
    std::string dump_lineages;
    std::string dump_usage;
    std::string weblogo_query;
    std::string usage_plot;
    std::string stats_plot;
    std::string logfile;
    std::string log_level = "INFO";
    double threshold_ratio = 0.1;
    int min_lineage_size = 10;
    int workers = 0;  // 0 means choose automatically
};

void print_usage() {
    std::cout <<
        "Clonal lineage analysis of antibody repertoire sequencing data,\n"
        "from NCBI IgBLAST output.\n\n"
        "Usage: antibody-repertoire --igblast PATH [options]\n\n"
        "  --igblast PATH          NCBI IgBLAST results TSV (-outfmt 19). Required.\n"
        "  --dump-lineages PATH    Write the full lineage partition as TSV.\n"
        "  --dump-usage PATH       Write the V gene usage table as TSV.\n"
        "  --weblogo-query PATH    Write the WebLogo query for the largest lineage.\n"
        "  --usage-plot PATH       Write the V gene usage bar chart as SVG.\n"
        "  --stats-plot PATH       Write the lineage statistics bar chart as SVG.\n"
        "  --workers N             Threads for the block scan. Default: automatic.\n"
        "  --threshold RATIO       Max Hamming distance as a fraction of CDR3\n"
        "                          length. Default: 0.1\n"
        "  --min-lineage-size N    Cutoff for the \"well represented\" tallies.\n"
        "                          Default: 10\n"
        "  --log-file PATH         Also write logs to this file.\n"
        "  --log-level LEVEL       DEBUG, INFO, WARNING or ERROR. Default: INFO\n"
        "  --help                  Show this message.\n";
}

// Returns false if the program should exit after printing help.
bool parse_args(int argc, char** argv, Options& options) {
    const auto value = [argc, argv](int& i, const char* flag) -> std::string {
        if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + flag);
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--help" || flag == "-h") { print_usage(); return false; }
        else if (flag == "--igblast")          options.igblast = value(i, "--igblast");
        else if (flag == "--dump-lineages")    options.dump_lineages = value(i, "--dump-lineages");
        else if (flag == "--dump-usage")       options.dump_usage = value(i, "--dump-usage");
        else if (flag == "--weblogo-query")    options.weblogo_query = value(i, "--weblogo-query");
        else if (flag == "--usage-plot")       options.usage_plot = value(i, "--usage-plot");
        else if (flag == "--stats-plot")       options.stats_plot = value(i, "--stats-plot");
        else if (flag == "--log-file")         options.logfile = value(i, "--log-file");
        else if (flag == "--log-level")        options.log_level = value(i, "--log-level");
        else if (flag == "--workers")          options.workers = std::stoi(value(i, "--workers"));
        else if (flag == "--threshold")        options.threshold_ratio = std::stod(value(i, "--threshold"));
        else if (flag == "--min-lineage-size") options.min_lineage_size = std::stoi(value(i, "--min-lineage-size"));
        else throw std::runtime_error("unrecognised argument: " + flag);
    }

    if (options.igblast.empty()) throw std::runtime_error("--igblast is required");
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        if (!parse_args(argc, argv, options)) return 0;
        ar::configure_logging(ar::parse_log_level(options.log_level), options.logfile);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n\n";
        print_usage();
        return 2;
    }

    try {
        const ar::Repertoire repertoire = ar::load_igblast_results(options.igblast, options.workers);
        const std::vector<ar::UniqueCdr3> nodes = ar::collapse_to_unique_cdr3(repertoire.reads);
        const std::vector<ar::Edge> edges =
            ar::build_edges(nodes, options.threshold_ratio, options.workers);
        const std::vector<ar::Lineage> lineages = ar::get_clonal_lineages(nodes, edges);

        const ar::Stats stats = ar::get_lineage_stats(repertoire.reads, nodes, lineages,
                                                      options.min_lineage_size);
        // Only the report goes to stdout; logs go to stderr, so redirecting
        // stdout yields a file that diffs against the reference pipeline.
        std::cout << ar::format_stats(stats);

        if (!options.dump_lineages.empty()) {
            ar::write_file(options.dump_lineages,
                           ar::format_lineage_table(repertoire.reads, nodes, lineages));
        }
        // Computed once and reused: the table and the chart are the same data.
        std::vector<std::pair<std::string, long long>> usage;
        if (!options.dump_usage.empty() || !options.usage_plot.empty()) {
            usage = ar::get_v_gene_usage(repertoire.reads, repertoire.genes);
        }
        if (!options.dump_usage.empty()) {
            ar::write_file(options.dump_usage, ar::format_usage_table(usage));
        }
        if (!options.usage_plot.empty()) {
            std::vector<ar::Bar> bars;
            bars.reserve(usage.size());
            for (const auto& [gene, count] : usage) {
                bars.push_back(ar::Bar{gene, static_cast<double>(count)});
            }
            ar::write_file(options.usage_plot,
                           ar::render_vertical_bars("V Gene Usage Statistics", "V Genes",
                                                    "Read Count", bars));
            ar::log_info("Plot written to " + options.usage_plot);
        }
        if (!options.stats_plot.empty()) {
            std::vector<ar::Bar> bars;
            bars.reserve(stats.size());
            for (const auto& [metric, value] : stats) {
                bars.push_back(ar::Bar{metric, static_cast<double>(value)});
            }
            ar::write_file(options.stats_plot,
                           ar::render_horizontal_bars("Clonal Lineage Statistics", "Value", bars));
            ar::log_info("Plot written to " + options.stats_plot);
        }
        if (!options.weblogo_query.empty()) {
            ar::write_file(options.weblogo_query,
                           ar::get_cdr3_aa_from_largest_lineage(repertoire.reads, nodes, lineages));
        }

        ar::log_info("Analysis complete.");
        return 0;
    } catch (const std::exception& error) {
        ar::log_error(error.what());
        return 1;
    }
}
