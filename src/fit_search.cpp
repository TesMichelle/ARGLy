#include "fit_search.h"
#include <likelihood_cache.h>
#include <nlopt.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace fit {
namespace {
constexpr std::array<size_t, 4> indices = {0, 3, 6, 7};

double value(const Axis& axis, size_t dimension, double unit) {
    if (axis.lower == axis.upper) return axis.lower;
    if (dimension >= 2)
        return std::exp(std::log(axis.lower) + unit * (std::log(axis.upper) - std::log(axis.lower)));
    return axis.lower + unit * (axis.upper - axis.lower);
}

double unit(const Axis& axis, size_t dimension, double x) {
    if (axis.lower == axis.upper) return 0;
    if (dimension >= 2)
        return (std::log(x) - std::log(axis.lower)) / (std::log(axis.upper) - std::log(axis.lower));
    return (x - axis.lower) / (axis.upper - axis.lower);
}

struct Point { std::array<double, 4> x; double score; };

struct Context {
    treellh::TreeEnsemble& ensemble;
    const Options& options;
    std::ostream& out;
    Result best;
    size_t calls = 0;
    const char* stage = "";
    std::vector<size_t> active;
    std::array<double, 4> seed{};
    ProgressEvent::Phase phase = ProgressEvent::Phase::grid;
    size_t phase_start_calls = 0;
    size_t phase_total = 0;
    size_t start_index = 0;
    size_t start_count = 0;

    void report(bool complete = false) const {
        if (!options.progress) return;
        options.progress({phase, calls - phase_start_calls, phase_total,
            start_index, start_count, complete, best.log_likelihood, best.parameters});
    }

    double evaluate(const std::array<double, 4>& x) {
        auto p = options.parameters;
        for (size_t i = 0; i < 4; ++i) p[indices[i]] = value(options.axes[i], i, x[i]);
        const double score = ensemble.evaluate(p);
        out << stage << '\t' << calls++ << '\t' << p[0] << '\t' << p[3] << '\t'
            << p[6] << '\t' << p[7] << '\t' << score << '\n';
        if (calls == 1 || score > best.log_likelihood) { best.parameters = p; best.log_likelihood = score; }
        report();
        return score;
    }
};

double objective(const std::vector<double>& x, std::vector<double>&, void* data) {
    auto& context = *static_cast<Context*>(data);
    auto point = context.seed;
    for (size_t i = 0; i < x.size(); ++i) point[context.active[i]] = x[i];
    const double score = context.evaluate(point);
    return std::isfinite(score) ? score : -1e300;
}

void optimize(Context& context, std::array<double, 4> start) {
    context.seed = start;
    if (context.active.empty()) { context.evaluate(start); return; }
    nlopt::opt opt(context.active.size() == 1 ? nlopt::LN_NELDERMEAD : nlopt::LN_BOBYQA,
        context.active.size());
    std::vector<double> x;
    for (size_t i : context.active) x.push_back(std::clamp(start[i], 0., 1.));
    opt.set_lower_bounds(0.);
    opt.set_upper_bounds(1.);
    opt.set_max_objective(objective, &context);
    opt.set_maxeval(static_cast<int>(context.options.max_evaluations));
    opt.set_xtol_abs(1e-5);
    opt.set_initial_step(0.05);
    double optimum;
    opt.optimize(x, optimum);
}
} // namespace

Result search(treellh::TreeEnsemble& ensemble, const Options& options,
    const char* mode, std::ostream& output) {
    const std::string strategy(mode);
    if (strategy != "fit-grid" && strategy != "fit-opt" && strategy != "fit-hybrid")
        throw std::invalid_argument("unknown fit mode");
    if (!options.max_evaluations || !options.starts || options.max_evaluations > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("optimizer evaluations and starts must be positive");
    for (size_t i = 0; i < 4; ++i) {
        const auto& a = options.axes[i];
        if (!std::isfinite(a.lower) || !std::isfinite(a.upper) || a.lower > a.upper || !a.count
            || (i == 0 && (a.lower < 0 || a.upper >= options.parameters[1]))
            || (i == 1 && (a.lower < 0 || a.upper > 1))
            || (i >= 2 && a.lower <= 0))
            throw std::invalid_argument("invalid fit axis or bounds");
    }
    Context context{ensemble, options, output,
        {options.parameters, -std::numeric_limits<double>::infinity(), 0}};
    for (size_t i = 0; i < 4; ++i)
        if (options.axes[i].lower != options.axes[i].upper) context.active.push_back(i);
    output << "stage\tevaluation\tmigration_time\tmigration_prob\tN_base\tN_ghost\tlog_llh\n";
    std::vector<Point> ranked;
    if (strategy != "fit-opt") {
        size_t total = 1;
        for (const auto& a : options.axes) {
            if (a.count > 1000000 / total) throw std::invalid_argument("fit grid exceeds one million points");
            total *= a.count;
        }
        ranked.reserve(total);
        context.stage = "grid";
        context.phase = ProgressEvent::Phase::grid;
        context.phase_start_calls = context.calls;
        context.phase_total = total;
        context.report();
        for (size_t p = 0; p < options.axes[1].count; ++p)
            for (size_t t = 0; t < options.axes[0].count; ++t)
                for (size_t n = 0; n < options.axes[2].count; ++n)
                    for (size_t g = 0; g < options.axes[3].count; ++g) {
                        std::array<size_t, 4> ix{t, p, n, g};
                        std::array<double, 4> x{};
                        for (size_t i = 0; i < 4; ++i)
                            x[i] = options.axes[i].count == 1 ? 0 : double(ix[i]) / (options.axes[i].count - 1);
                        ranked.push_back({x, context.evaluate(x)});
                    }
        context.report(true);
        std::sort(ranked.begin(), ranked.end(), [](const Point& a, const Point& b) { return a.score > b.score; });
    }
    if (strategy != "fit-grid") {
        context.stage = "optimizer";
        std::vector<std::array<double, 4>> seeds;
        if (strategy == "fit-opt") {
            std::array<double, 4> x{};
            for (size_t i = 0; i < 4; ++i)
                x[i] = std::clamp(unit(options.axes[i], i, options.parameters[indices[i]]), 0., 1.);
            seeds.push_back(x);
        } else {
            // Preserve the supplied starting point as one of the local runs.
            // A coarse grid can miss a narrow peak near this point.
            std::array<double, 4> initial{};
            for (size_t i = 0; i < 4; ++i)
                initial[i] = std::clamp(unit(options.axes[i], i, options.parameters[indices[i]]), 0., 1.);
            for (const auto& point : ranked) {
                if (!std::isfinite(point.score)) continue;
                bool separate = true;
                for (const auto& seed : seeds) {
                    double distance = 0;
                    for (size_t i = 0; i < 4; ++i) distance += std::pow(seed[i] - point.x[i], 2);
                    if (distance < 0.04) { separate = false; break; }
                }
                if (separate) seeds.push_back(point.x);
                if (seeds.size() == options.starts) break;
            }
            if (seeds.empty()) seeds.push_back(initial);
            if (options.starts > 1) {
                bool present = false;
                for (const auto& seed : seeds) {
                    double distance = 0;
                    for (size_t i = 0; i < 4; ++i) distance += std::pow(seed[i] - initial[i], 2);
                    if (distance < 1e-12) { present = true; break; }
                }
                if (!present) {
                    if (seeds.size() == options.starts) seeds.pop_back();
                    seeds.push_back(initial);
                }
            }
        }
        for (size_t i = 0; i < seeds.size(); ++i) {
            context.phase = ProgressEvent::Phase::optimizer;
            context.phase_start_calls = context.calls;
            context.phase_total = options.max_evaluations;
            context.start_index = i + 1;
            context.start_count = seeds.size();
            context.report();
            optimize(context, seeds[i]);
            context.report(true);
        }
    }
    context.best.evaluations = context.calls;
    output << "best\t" << context.calls << '\t' << context.best.parameters[0] << '\t'
        << context.best.parameters[3] << '\t' << context.best.parameters[6] << '\t'
        << context.best.parameters[7] << '\t' << context.best.log_likelihood << '\n';
    return context.best;
}
} // namespace fit
