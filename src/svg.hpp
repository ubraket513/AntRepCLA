// Bar charts rendered straight to SVG.
//
// The pipeline's two figures are bar charts. SVG is text, so drawing them
// needs no plotting library, no runtime binary and no build system beyond the
// one already here -- Matplot++, the obvious candidate, requires CMake to
// build and a Gnuplot installation at runtime, which is the kind of dependency
// this port exists to remove. The output is vector, opens in any browser, and
// diffs in version control.
#pragma once

#include <string>
#include <vector>

namespace ar {

struct Bar {
    std::string label;
    double value = 0;
};

// Vertical bars, tallest first, with rotated tick labels. Throws
// std::invalid_argument if `bars` is empty.
std::string render_vertical_bars(const std::string& title,
                                 const std::string& x_label,
                                 const std::string& y_label,
                                 const std::vector<Bar>& bars);

// Horizontal bars in the order given, each annotated with its value. Throws
// std::invalid_argument if `bars` is empty.
std::string render_horizontal_bars(const std::string& title,
                                   const std::string& x_label,
                                   const std::vector<Bar>& bars);

// XML-escape text destined for an SVG text node or attribute.
std::string escape_xml(const std::string& text);

// Viridis at `position` in [0, 1], as "#rrggbb". Matches the palette the
// reference plots used, so the figures stay recognisable.
std::string viridis(double position);

}  // namespace ar
