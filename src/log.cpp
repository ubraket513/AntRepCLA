#include "log.hpp"

#include <array>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>

namespace ar {
namespace {

struct State {
    LogLevel level = LogLevel::info;
    std::ofstream file;
    std::mutex mutex;
};

State& state() {
    static State instance;
    return instance;
}

const char* level_name(LogLevel level) {
    switch (level) {
        case LogLevel::debug:   return "debug";
        case LogLevel::info:    return "info";
        case LogLevel::warning: return "warning";
        case LogLevel::error:   return "error";
    }
    return "info";
}

std::string now_timestamp() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t seconds = std::chrono::system_clock::to_time_t(now);

    std::tm broken{};
    // localtime_r rather than localtime: the scan runs under OpenMP and the
    // non-reentrant form shares a static buffer.
    ::localtime_r(&seconds, &broken);

    std::array<char, 32> buffer{};
    const std::size_t written =
        std::strftime(buffer.data(), buffer.size(), "%Y-%m-%d %H:%M:%S", &broken);
    return std::string(buffer.data(), written);
}

void emit(LogLevel level, const std::string& message) {
    State& current = state();
    if (!should_log(current.level, level)) return;

    const std::string line = format_log_line(now_timestamp(), level, message);

    std::lock_guard<std::mutex> guard(current.mutex);
    // stderr, so that stdout carries only the report and stays diffable.
    std::cerr << line << '\n';
    if (current.file.is_open()) current.file << line << '\n';
}

}  // namespace

bool should_log(LogLevel configured, LogLevel message_level) {
    return static_cast<int>(message_level) >= static_cast<int>(configured);
}

std::string format_log_line(std::string_view timestamp, LogLevel level,
                            std::string_view message) {
    std::string name = level_name(level);
    if (name.size() < 8) name.append(8 - name.size(), ' ');

    std::string line;
    line.reserve(timestamp.size() + name.size() + message.size() + 8);
    line.append(timestamp);
    line.append(" | ");
    line.append(name);
    line.append(" | ");
    line.append(message);
    return line;
}

void configure_logging(LogLevel level, const std::string& logfile) {
    State& current = state();
    std::lock_guard<std::mutex> guard(current.mutex);
    current.level = level;
    if (!logfile.empty()) {
        current.file.open(logfile, std::ios::app);
        if (!current.file) throw std::runtime_error("Cannot open log file '" + logfile + "'.");
    }
}

void log_debug(const std::string& message)   { emit(LogLevel::debug, message); }
void log_info(const std::string& message)    { emit(LogLevel::info, message); }
void log_warning(const std::string& message) { emit(LogLevel::warning, message); }
void log_error(const std::string& message)   { emit(LogLevel::error, message); }

LogLevel parse_log_level(const std::string& name) {
    if (name == "DEBUG")   return LogLevel::debug;
    if (name == "INFO")    return LogLevel::info;
    if (name == "WARNING") return LogLevel::warning;
    if (name == "ERROR")   return LogLevel::error;
    throw std::invalid_argument("unknown log level: " + name);
}

}  // namespace ar
