#pragma once

#include "treellh.h"
#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>

namespace treellh {

// Statistics immediately above an event: its coalescence belongs to the tail.
struct LikelihoodStatistics {
    std::array<double, 2> exposure{};
    std::array<uint32_t, 2> coalescences{};
    std::array<int, 2> lineages{};
};

struct CacheMemoryEstimate {
    size_t trees = 0;
    size_t scenarios = 0;
    size_t records = 0;
    size_t statistics_bytes = 0;
    size_t resident_bytes = 0;
};

using CacheProgress = std::function<void(std::string_view, size_t, size_t)>;

class TreeLikelihoodCache {
public:
    double left = 0;
    double right = 0;
    size_t scenario_count() const { return scenarios_; }
    size_t record_count() const;
    size_t metadata_bytes() const;
    double evaluate(std::span<const double> parameters);

private:
    friend class TreeEnsemble;
    void prepare(const tsk_tree_t&, Scenario_Computer&);
    void build(const tsk_tree_t&, Scenario_Computer&);
    size_t locate(double time);
    LikelihoodStatistics at(size_t scenario, size_t boundary, double time) const;

    // Local tree times, descending; include split boundaries and zero.
    // All equal-time vertices share a boundary record.
    std::vector<double> times_;
    std::vector<LikelihoodStatistics> statistics_;
    LikelihoodStatistics root_;
    double split1_ = 0;
    double split2_ = 0;
    size_t scenarios_ = 0;
    size_t cursor_ = 0;
    bool impossible_ = false;
    bool ready_ = false;
};

// The tree sequence must outlive this object. Construction only estimates and
// prepares metadata; build() explicitly allocates and fills the scenario cache.
// Queries may change migration time/probability and N, but not split times.
class TreeEnsemble {
public:
    TreeEnsemble(const tsk_treeseq_t&, std::span<const double> parameters,
        const std::vector<int>& sample_population, int admixed_index, int outgroup_index,
        bool defer_prepare = false, size_t tree_limit = 0);
    void prepare(const CacheProgress& progress = {});
    CacheMemoryEstimate memory_estimate() const;
    void build(const CacheProgress& progress = {});
    void save(const std::string& path, const std::string& signature,
        const CacheProgress& progress = {}) const;
    void load(const std::string& path, const std::string& signature, size_t max_bytes = 0,
        const CacheProgress& progress = {});
    double evaluate(std::span<const double> parameters);
    std::vector<double> evaluate_trees(std::span<const double> parameters);
    const std::vector<TreeLikelihoodCache>& trees() const { return trees_; }

private:
    void validate_query(std::span<const double>) const;
    const tsk_treeseq_t& ts_;
    std::vector<double> parameters_;
    std::vector<int> sample_population_;
    int admixed_index_;
    int outgroup_index_;
    size_t selected_count_ = 0;
    std::vector<size_t> selected_indices_;
    size_t selected_index(size_t ordinal) const;
    bool ready_ = false;
    bool prepared_ = false;
    std::vector<TreeLikelihoodCache> trees_;
};

}
