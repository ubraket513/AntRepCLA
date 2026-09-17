// Logging.
//
// Small enough to own outright. The program emits a handful of progress lines,
// which does not justify a logging library: spdlog was 1.2 MB of vendored
// headers and ~7.7 s of compile time per translation unit that included it.
#pragma once

#include <string>
#include <string_view>

namespace ar {

enum class LogLevel { debug, info, warning, error };

// Configure the global logger. `logfile` may be empty for stderr only.
// Throws std::runtime_error if the log file cannot be opened.
void configure_logging(LogLevel level, const std::string& logfile = "");

void log_debug(const std::string& message);
void log_info(const std::string& message);
void log_warning(const std::string& message);
void log_error(const std::string& message);

LogLevel parse_log_level(const std::string& name);

// Whether a message at `message_level` is emitted when the logger is
// configured at `configured`.
bool should_log(LogLevel configured, LogLevel message_level);

// One rendered line: timestamp, level padded to eight columns, message.
// Pure, so the layout is testable without touching the clock.
std::string format_log_line(std::string_view timestamp, LogLevel level,
                            std::string_view message);

}  // namespace ar
