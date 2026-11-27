#pragma once
#include <treellh.h>
#include <array>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace epochfit {
struct Parameters {
    double time = 0, probability = 0;
    std::vector<double> ghost_sizes;
    std::vector<double> base_sizes; // Empty means fixed Nbase from --params.
    std::vector<double> ghost_right_sizes, base_right_sizes;
    std::vector<double> ghost_cuts, base_cuts; // Movable two-step mode only.
};
struct RateEvidence {
    double left_exposure = 0, right_exposure = 0;
    std::vector<std::pair<double,double>> events; // Position in [0,1], posterior mass.
    std::vector<std::pair<double,double>> cut_events; // Absolute time, posterior mass.
    std::vector<std::array<double,3>> spans; // Absolute low/high, posterior pair rate.
};
struct Sufficient {
    std::vector<double> exposure, coalescences;
    std::vector<double> base_exposure, base_coalescences;
    std::vector<RateEvidence> ghost_rates, base_rates;
    double log_likelihood = 0;
};
struct LikelihoodDerivatives {
    double log_likelihood = 0;
    double time = 0; // d log-likelihood / d migration time.
    double probability = 0; // d log-likelihood / d migration probability.
};
struct KernelEstimate {
    double time = 0, coalescences = 0, exposure = 0;
    double size = 0;
};
std::vector<double> balanced_epoch_boundaries(
    const std::vector<std::pair<double,double>>& event_masses,
    size_t epochs, double split_time);
std::vector<std::pair<double,double>> raw_coalescence_counts(
    const tsk_treeseq_t& trees, double split_time);
std::vector<double> coalescences_per_epoch(
    const std::vector<std::pair<double,double>>& events,
    const std::vector<double>& boundaries);
struct RawEpochStatistics {
    double coalescences = 0;
    double lineage_time = 0; // Integral of the number of active lineages.
    double pair_time = 0;    // Integral of the number of unordered lineage pairs.
};
std::array<std::vector<RawEpochStatistics>,2> raw_epoch_statistics(
    const tsk_treeseq_t& trees, const std::vector<double>& base_boundaries,
    const std::vector<double>& ghost_boundaries);
struct Step {
    std::array<int, 2> lineages{};
    std::array<unsigned, 2> coalescences{};
};
struct TreeCache {
    std::vector<double> times; // Descending, shared by all colourings, includes epoch boundaries.
    std::vector<std::array<size_t, 2>> epochs; // Base and ghost epoch for each interval.
    std::vector<std::array<unsigned, 2>> totals; // Scenario-independent lineages and events.
    std::vector<unsigned char> ghost_steps; // Boundary-major ghost counts, two values per scenario.
    unsigned count_bytes = 1;
    std::vector<double> prefix; // Only rows reachable by the migration-time bounds.
    size_t prefix_first = 0, prefix_last = 0; // Inclusive row indices.
    double tail_prefix = 0; // Scenario-zero prefix at time zero.
    size_t scenarios = 0;
    size_t cursor = 0;
    double root_log_likelihood = 0;
    bool impossible = false;
    Step step(size_t boundary, size_t scenario) const;
    void set_step(size_t boundary, size_t scenario, const Step& value);
    double prefix_at(size_t boundary, size_t scenario) const;
};
class Ensemble {
public:
    Ensemble(const tsk_treeseq_t&, const std::vector<double>&,
        const std::vector<int>& populations, int admixed, int outgroup,
        size_t epochs, size_t tree_limit = 0, double memory_limit_mb = 0,
        const std::function<void(size_t,size_t,size_t)>& progress = {},
        const std::vector<double>& epoch_boundaries = {},
        size_t base_epochs = 0, const std::vector<double>& base_epoch_boundaries = {},
        std::array<double,2> time_bounds = {0, std::numeric_limits<double>::infinity()},
        bool even_tree_indices = false);
    // Append independently selected trees; call set_sizes again before evaluating.
    void append(const tsk_treeseq_t&, const std::vector<double>&,
        const std::vector<int>& populations, int admixed, int outgroup,
        size_t tree_limit = 0, double memory_limit_mb = 0,
        const std::function<void(size_t,size_t,size_t)>& progress = {},
        bool even_tree_indices = false);
    void set_sizes(const std::vector<double>&);
    void set_sizes(const std::vector<double>& base, const std::vector<double>& ghost);
    void set_linear_sizes(const std::vector<double>& base_left,
        const std::vector<double>& base_right,const std::vector<double>& ghost_left,
        const std::vector<double>& ghost_right);
    void set_split_sizes(const std::vector<double>& base_left,const std::vector<double>& base_right,
        const std::vector<double>& base_cuts,const std::vector<double>& ghost_left,
        const std::vector<double>& ghost_right,const std::vector<double>& ghost_cuts);
    void set_split_base_sizes(const std::vector<double>& base_left,
        const std::vector<double>& base_right,const std::vector<double>& base_cuts,
        const std::vector<double>& ghost);
    bool linear_rates() const { return linear_rates_; }
    // Zero selects automatic thread count; one forces serial evaluation.
    void set_threads(size_t threads);
    double evaluate(double time, double probability);
    // Observed-likelihood derivatives at fixed N. Requires 0 < probability < 1
    // and time strictly between cached event/epoch boundaries.
    LikelihoodDerivatives derivatives(double time, double probability);
    std::vector<LikelihoodDerivatives> derivatives_by_tree(double time, double probability);
    Sufficient expectation(double time, double probability);
    std::array<std::vector<KernelEstimate>,2> kernel_profile(
        double time, double probability, double bandwidth, size_t points);
    std::array<std::vector<std::pair<double,double>>,2>
        coalescence_masses(double time, double probability);
    size_t bytes() const;
    size_t tree_count() const { return trees_.size(); }
    size_t scenario_count() const;
    const std::vector<double>& sizes() const { return sizes_; }
    const std::vector<double>& base_sizes() const { return base_sizes_; }
    const std::vector<double>& boundaries() const { return boundaries_; }
    const std::vector<double>& base_boundaries() const { return base_boundaries_; }
    double split2() const { return split2_; }
    double migration_limit() const { return has_outgroup_ ? split2_ : split1_; }
    bool has_outgroup() const { return has_outgroup_; }
private:
    std::vector<TreeCache> trees_;
    std::vector<double> boundaries_, base_boundaries_, sizes_, base_sizes_;
    std::vector<double> right_sizes_, base_right_sizes_;
    std::vector<double> cuts_, base_cuts_;
    size_t resident_bytes_ = 0;
    size_t threads_ = 0;
    double base_size_, split2_, split1_;
    std::array<double,2> time_bounds_{};
    int admixed_index_ = 0, outgroup_index_ = 1;
    bool has_outgroup_ = false;
    bool linear_rates_ = false;
    bool split_base_ = false, split_ghost_ = false;
    void set_profile(const std::vector<double>& base_left,const std::vector<double>& base_right,
        const std::vector<double>& ghost_left,const std::vector<double>& ghost_right,
        const std::vector<double>& base_cuts,const std::vector<double>& ghost_cuts,
        bool linear,bool split_base,bool split_ghost);
    double rate(size_t pop,size_t epoch,double time) const;
    double rate_integral(size_t pop,size_t epoch,double low,double high) const;
    size_t locate(TreeCache&, double);
    double partial(const TreeCache&, size_t, size_t, double) const;
    double score(TreeCache&, size_t, size_t, double, double) const;
    void validate_query(double,double) const;
};
}
