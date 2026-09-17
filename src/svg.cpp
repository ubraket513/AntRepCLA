#include "svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace ar {
namespace {

// Viridis sampled at nine points; linear interpolation between them is close
// enough for a bar chart and avoids carrying the full 256-entry table.
constexpr std::array<std::array<int, 3>, 9> kViridis{{
    {68, 1, 84}, {72, 40, 120}, {62, 74, 137}, {49, 104, 142}, {38, 130, 142},
    {31, 158, 137}, {53, 183, 121}, {109, 205, 89}, {180, 222, 44},
}};

std::string number(double value) {
    std::array<char, 32> buffer{};
    std::snprintf(buffer.data(), buffer.size(), "%.2f", value);
    return buffer.data();
}

std::string integer_text(double value) {
    std::array<char, 32> buffer{};
    std::snprintf(buffer.data(), buffer.size(), "%.0f", value);
    return buffer.data();
}

void require_non_empty(const std::vector<Bar>& bars) {
    if (bars.empty()) throw std::invalid_argument("bars must not be empty.");
}

double largest_value(const std::vector<Bar>& bars) {
    double largest = 0;
    for (const Bar& bar : bars) largest = std::max(largest, bar.value);
    return largest;
}

// Shared chrome: the document header and a style block, so both charts look
// like one family and the CSS is not repeated per element.
std::string open_svg(double width, double height, const std::string& title) {
    std::string out;
    out += "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + number(width) +
           "\" height=\"" + number(height) + "\" viewBox=\"0 0 " + number(width) + " " +
           number(height) + "\" font-family=\"DejaVu Sans, Helvetica, Arial, sans-serif\">\n";
    out += "<style>\n"
           "  .bg { fill: #ffffff; }\n"
           "  .title { font-size: 17px; font-weight: bold; fill: #1a1a1a; }\n"
           "  .axis-label { font-size: 13px; font-weight: bold; fill: #1a1a1a; }\n"
           "  .tick { font-size: 10px; fill: #333333; }\n"
           "  .value { font-size: 11px; fill: #1a1a1a; }\n"
           "  .grid { stroke: #b0b0b0; stroke-width: 1; stroke-dasharray: 4 3; opacity: 0.7; }\n"
           "  .axis { stroke: #1a1a1a; stroke-width: 1.2; }\n"
           "</style>\n";
    out += "<rect class=\"bg\" width=\"100%\" height=\"100%\"/>\n";
    out += "<text class=\"title\" x=\"" + number(width / 2) +
           "\" y=\"28\" text-anchor=\"middle\">" + escape_xml(title) + "</text>\n";
    return out;
}

// Round a maximum up to a readable axis bound (1, 2 or 5 times a power of ten).
double nice_bound(double value) {
    if (value <= 0) return 1;
    const double magnitude = std::pow(10.0, std::floor(std::log10(value)));
    const double normalised = value / magnitude;
    double step = 10;
    if (normalised <= 1) step = 1;
    else if (normalised <= 2) step = 2;
    else if (normalised <= 5) step = 5;
    return step * magnitude;
}

}  // namespace

std::string escape_xml(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (const char character : text) {
        switch (character) {
            case '&':  out += "&amp;";  break;
            case '<':  out += "&lt;";   break;
            case '>':  out += "&gt;";   break;
            case '"':  out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default:   out.push_back(character);
        }
    }
    return out;
}

std::string viridis(double position) {
    position = std::clamp(position, 0.0, 1.0);
    const double scaled = position * static_cast<double>(kViridis.size() - 1);
    const std::size_t low = static_cast<std::size_t>(std::floor(scaled));
    const std::size_t high = std::min(low + 1, kViridis.size() - 1);
    const double fraction = scaled - static_cast<double>(low);

    std::array<char, 16> buffer{};
    const auto channel = [&](int index) {
        return static_cast<int>(std::lround(kViridis[low][index] * (1 - fraction) +
                                            kViridis[high][index] * fraction));
    };
    std::snprintf(buffer.data(), buffer.size(), "#%02x%02x%02x",
                  channel(0), channel(1), channel(2));
    return buffer.data();
}

std::string render_vertical_bars(const std::string& title,
                                 const std::string& x_label,
                                 const std::string& y_label,
                                 const std::vector<Bar>& bars) {
    require_non_empty(bars);

    std::vector<Bar> ordered = bars;
    std::stable_sort(ordered.begin(), ordered.end(),
                     [](const Bar& a, const Bar& b) { return a.value > b.value; });

    // Width follows the series so labels never overlap, however many genes
    // there are.
    const double left = 80, right = 30, top = 50, bottom = 130;
    const double bar_slot = std::max(18.0, 900.0 / static_cast<double>(ordered.size()));
    const double plot_width = bar_slot * static_cast<double>(ordered.size());
    const double plot_height = 420;
    const double width = left + plot_width + right;
    const double height = top + plot_height + bottom;

    const double bound = nice_bound(largest_value(ordered));
    std::string out = open_svg(width, height, title);

    // Horizontal grid and y ticks.
    for (int step = 0; step <= 5; ++step) {
        const double value = bound * step / 5.0;
        const double y = top + plot_height - (value / bound) * plot_height;
        out += "<line class=\"grid\" x1=\"" + number(left) + "\" y1=\"" + number(y) +
               "\" x2=\"" + number(left + plot_width) + "\" y2=\"" + number(y) + "\"/>\n";
        out += "<text class=\"tick\" x=\"" + number(left - 8) + "\" y=\"" + number(y + 4) +
               "\" text-anchor=\"end\">" + integer_text(value) + "</text>\n";
    }

    for (std::size_t i = 0; i < ordered.size(); ++i) {
        const double fraction = ordered.size() == 1
            ? 0.0
            : static_cast<double>(i) / static_cast<double>(ordered.size() - 1);
        const double bar_height = bound > 0 ? (ordered[i].value / bound) * plot_height : 0.0;
        const double x = left + bar_slot * static_cast<double>(i) + bar_slot * 0.15;
        const double y = top + plot_height - bar_height;

        out += "<rect class=\"bar\" x=\"" + number(x) + "\" y=\"" + number(y) +
               "\" width=\"" + number(bar_slot * 0.7) + "\" height=\"" + number(bar_height) +
               "\" fill=\"" + viridis(fraction) + "\"/>\n";

        const double label_x = x + bar_slot * 0.35;
        const double label_y = top + plot_height + 10;
        out += "<text class=\"tick\" x=\"" + number(label_x) + "\" y=\"" + number(label_y) +
               "\" text-anchor=\"end\" transform=\"rotate(-90 " + number(label_x) + " " +
               number(label_y) + ")\">" + escape_xml(ordered[i].label) + "</text>\n";
    }

    out += "<line class=\"axis\" x1=\"" + number(left) + "\" y1=\"" + number(top + plot_height) +
           "\" x2=\"" + number(left + plot_width) + "\" y2=\"" + number(top + plot_height) + "\"/>\n";
    out += "<text class=\"axis-label\" x=\"" + number(left + plot_width / 2) + "\" y=\"" +
           number(height - 12) + "\" text-anchor=\"middle\">" + escape_xml(x_label) + "</text>\n";
    out += "<text class=\"axis-label\" x=\"18\" y=\"" + number(top + plot_height / 2) +
           "\" text-anchor=\"middle\" transform=\"rotate(-90 18 " +
           number(top + plot_height / 2) + ")\">" + escape_xml(y_label) + "</text>\n";
    out += "</svg>\n";
    return out;
}

std::string render_horizontal_bars(const std::string& title,
                                   const std::string& x_label,
                                   const std::vector<Bar>& bars) {
    require_non_empty(bars);

    // Metric names are long, so the left margin is generous and the labels sit
    // outside the plot rather than on top of the bars.
    const double left = 380, right = 90, top = 50, bottom = 60;
    const double row = 44;
    const double plot_width = 520;
    const double plot_height = row * static_cast<double>(bars.size());
    const double width = left + plot_width + right;
    const double height = top + plot_height + bottom;

    const double bound = nice_bound(largest_value(bars));
    std::string out = open_svg(width, height, title);

    for (int step = 0; step <= 5; ++step) {
        const double value = bound * step / 5.0;
        const double x = left + (value / bound) * plot_width;
        out += "<line class=\"grid\" x1=\"" + number(x) + "\" y1=\"" + number(top) +
               "\" x2=\"" + number(x) + "\" y2=\"" + number(top + plot_height) + "\"/>\n";
        out += "<text class=\"tick\" x=\"" + number(x) + "\" y=\"" +
               number(top + plot_height + 18) + "\" text-anchor=\"middle\">" +
               integer_text(value) + "</text>\n";
    }

    for (std::size_t i = 0; i < bars.size(); ++i) {
        const double fraction = bars.size() == 1
            ? 0.0
            : static_cast<double>(i) / static_cast<double>(bars.size() - 1);
        const double bar_width = bound > 0 ? (bars[i].value / bound) * plot_width : 0.0;
        const double y = top + row * static_cast<double>(i) + row * 0.18;

        out += "<rect class=\"bar\" x=\"" + number(left) + "\" y=\"" + number(y) +
               "\" width=\"" + number(bar_width) + "\" height=\"" + number(row * 0.64) +
               "\" fill=\"" + viridis(fraction) + "\"/>\n";
        out += "<text class=\"tick\" x=\"" + number(left - 10) + "\" y=\"" +
               number(y + row * 0.45) + "\" text-anchor=\"end\">" +
               escape_xml(bars[i].label) + "</text>\n";
        out += "<text class=\"value\" x=\"" + number(left + bar_width + 8) + "\" y=\"" +
               number(y + row * 0.45) + "\">" + integer_text(bars[i].value) + "</text>\n";
    }

    out += "<line class=\"axis\" x1=\"" + number(left) + "\" y1=\"" + number(top) +
           "\" x2=\"" + number(left) + "\" y2=\"" + number(top + plot_height) + "\"/>\n";
    out += "<text class=\"axis-label\" x=\"" + number(left + plot_width / 2) + "\" y=\"" +
           number(height - 14) + "\" text-anchor=\"middle\">" + escape_xml(x_label) + "</text>\n";
    out += "</svg>\n";
    return out;
}

}  // namespace ar
