#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <vector>

namespace treellh { class TreeEnsemble; }

namespace fit {

struct Axis {
    double lower;
    double upper;
    size_t count;
};

struct ProgressEvent {
    enum class Phase { grid, optimizer } phase;
    size_t done = 0;
    size_t total = 0;
    size_t start = 0; // One-based optimizer start; zero for the grid.
    size_t starts = 0;
    bool complete = false;
    double best_log_likelihood = 0;
    std::array<double, 9> best_parameters{};
};

struct Options {
    std::array<double, 9> parameters;
    std::array<Axis, 4> axes; // migration time, probability, Nbase, Nghost
    size_t max_evaluations = 250;
    size_t starts = 3;
    std::function<void(const ProgressEvent&)> progress;
};

struct Result {
    std::array<double, 9> parameters;
    double log_likelihood;
    size_t evaluations;
};

// Each evaluation is written as TSV. The returned point is the best evaluated one.
Result search(treellh::TreeEnsemble& ensemble, const Options& options,
    const char* mode, std::ostream& output);

} // namespace fit
