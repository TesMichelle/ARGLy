#pragma once
#include <treellh.h>
#include <array>
#include <functional>
#include <utility>
#include <vector>

namespace epochfit {
struct Parameters {
    double time = 0, probability = 0;
    std::vector<double> ghost_sizes;
    std::vector<double> base_sizes; // Empty means fixed Nbase from --params.
};
struct Sufficient {
    std::vector<double> exposure, coalescences;
    std::vector<double> base_exposure, base_coalescences;
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
    std::vector<Step> steps; // Boundary-major; events at the upper endpoint, lineages below it.
    std::vector<double> prefix; // Log likelihood just ABOVE each event, rebuilt when N changes.
    size_t scenarios = 0;
    size_t cursor = 0;
    double root_log_likelihood = 0;
    bool impossible = false;
};
class Ensemble {
public:
    Ensemble(const tsk_treeseq_t&, const std::vector<double>&,
        const std::vector<int>& populations, int admixed, int outgroup,
        size_t epochs, size_t tree_limit = 0, double memory_limit_mb = 0,
        const std::function<void(size_t,size_t,size_t)>& progress = {},
        const std::vector<double>& epoch_boundaries = {},
        size_t base_epochs = 0, const std::vector<double>& base_epoch_boundaries = {});
    void set_sizes(const std::vector<double>&);
    void set_sizes(const std::vector<double>& base, const std::vector<double>& ghost);
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
    size_t resident_bytes_ = 0;
    double base_size_, split2_, split1_;
    bool has_outgroup_ = false;
    size_t locate(TreeCache&, double);
    double partial(const TreeCache&, size_t, size_t, double) const;
    double score(TreeCache&, size_t, size_t, double, double) const;
    void validate_query(double,double) const;
};
}
