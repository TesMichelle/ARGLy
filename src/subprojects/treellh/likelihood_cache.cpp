#include "likelihood_cache.h"

#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace treellh {
namespace {
size_t checked_add(size_t a, size_t b) {
    if (b > std::numeric_limits<size_t>::max() - a)
        throw std::overflow_error("likelihood cache size overflows size_t");
    return a + b;
}
size_t checked_mul(size_t a, size_t b) {
    if (a && b > std::numeric_limits<size_t>::max() / a)
        throw std::overflow_error("likelihood cache size overflows size_t");
    return a * b;
}
double pairs(int k) { return static_cast<double>(k) * (k - 1) / 2; }
void advance(LikelihoodStatistics& s, double dt) {
    for (size_t p = 0; p < 2; ++p) s.exposure[p] += pairs(s.lineages[p]) * dt;
}
double log_add(double a, double b) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    if (a == neg_inf) return b;
    if (b == neg_inf) return a;
    return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b)));
}
void validate_parameters(std::span<const double> p) {
    if (p.size() != 9) throw std::invalid_argument("cache requires nine parameters");
    for (double x : p)
        if (!std::isfinite(x)) throw std::invalid_argument("parameters must be finite");
    if (!(0 <= p[0] && p[0] < p[1] && p[1] < p[2]))
        throw std::invalid_argument("cache requires 0 <= migration time < split2 < split1");
    if (p[3] < 0 || p[3] > 1 || p[6] <= 0 || p[7] <= 0 || p[8] <= 0)
        throw std::invalid_argument("migration probability must be in [0,1] and N positive");
}
struct TreeGuard {
    tsk_tree_t tree{};
    explicit TreeGuard(const tsk_treeseq_t& ts) {
        int ret = tsk_tree_init(&tree, &ts, 0);
        if (ret < 0) throw std::runtime_error(tsk_strerror(ret));
    }
    ~TreeGuard() { tsk_tree_free(&tree); }
};
void check(int ret) {
    if (ret < 0) throw std::runtime_error(tsk_strerror(ret));
}

constexpr char cache_magic[8] = {'A', 'R', 'G', 'L', 'Y', 'T', 'C', '3'};
static_assert(sizeof(LikelihoodStatistics) == 32);
static_assert(std::is_trivially_copyable_v<LikelihoodStatistics>);

struct CacheWriter {
    std::ofstream& out;
    uint64_t hash = 14695981039346656037ULL;
    void bytes(const void* data, size_t count) {
        out.write(static_cast<const char*>(data), static_cast<std::streamsize>(count));
        if (!out) throw std::runtime_error("failed writing tree statistics cache");
        const auto* ptr = static_cast<const unsigned char*>(data);
        for (size_t i = 0; i < count; ++i) hash = (hash ^ ptr[i]) * 1099511628211ULL;
    }
    template<class T> void pod(const T& value) { bytes(&value, sizeof(value)); }
};
struct CacheReader {
    std::ifstream& in;
    uint64_t hash = 14695981039346656037ULL;
    void bytes(void* data, size_t count) {
        in.read(static_cast<char*>(data), static_cast<std::streamsize>(count));
        if (!in) throw std::runtime_error("truncated tree statistics cache");
        const auto* ptr = static_cast<const unsigned char*>(data);
        for (size_t i = 0; i < count; ++i) hash = (hash ^ ptr[i]) * 1099511628211ULL;
    }
    template<class T> T pod() { T value{}; bytes(&value, sizeof(value)); return value; }
};
}

size_t TreeLikelihoodCache::record_count() const {
    return checked_mul(scenarios_, times_.size());
}
size_t TreeLikelihoodCache::metadata_bytes() const {
    return checked_mul(times_.capacity(), sizeof(double));
}

void TreeLikelihoodCache::prepare(const tsk_tree_t& tree, Scenario_Computer& c) {
    if (tsk_tree_get_num_roots(&tree) != 1)
        throw std::invalid_argument("cache requires a single-root tree");
    split1_ = c.split_1_time_;
    split2_ = c.split_2_time_;
    left = tree.interval.left;
    right = tree.interval.right;
    for (auto& node : c.nodes_) node.population = -1;
    c.sort_by_time(tree);
    for (auto id : c.time_ordered_tree_nodes_ids_) {
        int children = 0;
        for (auto child = tree.left_child[id]; child != TSK_NULL; child = tree.right_sib[child])
            ++children;
        const bool sample = c.ts_.tables->nodes.flags[id] & TSK_NODE_IS_SAMPLE;
        if ((children == 0 && (!sample || c.nodes_[id].time != 0))
            || (children != 0 && (children != 2 || sample)))
            throw std::invalid_argument("cache requires binary trees with present-day sample leaves");
    }
    c.compute_lower_nodes_split_1(tree);
    c.compute_fixed_nodes(tree);
    impossible_ = c.check_imposible_split_2(tree) != 0;
    if (impossible_) { scenarios_ = 0; return; }
    const int bits = c.compute_log_scenarios_num();
    if (bits >= static_cast<int>(std::numeric_limits<size_t>::digits) || bits >= 63)
        throw std::overflow_error("too many migration scenarios for cache");
    scenarios_ = size_t{1} << bits;
    times_ = {split1_, split2_, 0};
    for (auto id : c.time_ordered_tree_nodes_ids_)
        if (c.nodes_[id].time > 0 && c.nodes_[id].time <= split1_)
            times_.push_back(c.nodes_[id].time);
    std::sort(times_.begin(), times_.end(), std::greater<double>());
    times_.erase(std::unique(times_.begin(), times_.end()), times_.end());
    checked_mul(record_count(), sizeof(LikelihoodStatistics));
}

void TreeLikelihoodCache::build(const tsk_tree_t& tree, Scenario_Computer& c) {
    if (ready_) return;
    if (impossible_) { ready_ = true; return; }
    // Colour all the way to zero once per scenario. Query-time migration is
    // represented by stitching a structured prefix to the all-base suffix.
    c.migration_time_ = 0;
    for (auto& node : c.nodes_) node.population = -1;
    c.sort_by_time(tree);
    c.compute_lower_nodes_split_1(tree);
    c.compute_fixed_nodes(tree);
    c.compute_log_scenarios_num();
    root_ = {};
    root_.lineages[0] = 1;
    double previous = std::max(split1_, c.nodes_[c.time_ordered_tree_nodes_ids_.front()].time);
    for (auto id : c.time_ordered_tree_nodes_ids_) {
        const double time = c.nodes_[id].time;
        if (time <= split1_) break;
        advance(root_, previous - time);
        ++root_.coalescences[0];
        ++root_.lineages[0];
        previous = time;
    }
    advance(root_, previous - split1_);
    statistics_.resize(record_count());
    for (size_t scenario = 0; scenario < scenarios_; ++scenario) {
        c.set_scenario(tree, scenario);
        LikelihoodStatistics current;
        current.lineages = {c.sc_lineages_num_[0], c.sc_lineages_num_[1]};
        previous = split1_;
        size_t node_i = 0;
        while (node_i < c.time_ordered_tree_nodes_ids_.size()
            && c.nodes_[c.time_ordered_tree_nodes_ids_[node_i]].time > split1_) ++node_i;
        for (size_t i = 0; i < times_.size(); ++i) {
            const double time = times_[i];
            advance(current, previous - time);
            statistics_[i * scenarios_ + scenario] = current;
            if (time == split2_) current.lineages[0] -= c.lower_split_2_out_linages_num_;
            while (node_i < c.time_ordered_tree_nodes_ids_.size()
                && c.nodes_[c.time_ordered_tree_nodes_ids_[node_i]].time == time) {
                const auto id = c.time_ordered_tree_nodes_ids_[node_i++];
                const int pop = c.nodes_[id].population;
                if (time > 0 && pop != c.outgroup_pop_index_) {
                    if (pop < 0 || pop > 1) throw std::logic_error("uncoloured cache event");
                    ++current.coalescences[pop];
                    ++current.lineages[pop];
                }
            }
            previous = time;
        }
    }
    ready_ = true;
}

size_t TreeLikelihoodCache::locate(double time) {
    auto contains = [&](size_t i) {
        return i < times_.size() && times_[i] <= time && (i == 0 || time < times_[i - 1]);
    };
    if (contains(cursor_)) return cursor_;
    if (cursor_ > 0 && contains(cursor_ - 1)) return --cursor_;
    if (contains(cursor_ + 1)) return ++cursor_;
    cursor_ = std::lower_bound(times_.begin(), times_.end(), time, std::greater<double>()) - times_.begin();
    return cursor_;
}

LikelihoodStatistics TreeLikelihoodCache::at(size_t scenario, size_t boundary, double time) const {
    auto result = statistics_[boundary * scenarios_ + scenario];
    advance(result, times_[boundary] - time);
    return result;
}

double TreeLikelihoodCache::evaluate(std::span<const double> parameters) {
    validate_parameters(parameters);
    if (!ready_) throw std::logic_error("build likelihood cache before querying it");
    if (parameters[1] != split2_ || parameters[2] != split1_)
        throw std::invalid_argument("split times changed: rebuild likelihood cache");
    if (impossible_) return -std::numeric_limits<double>::infinity();
    const size_t boundary = locate(parameters[0]);
    const auto baseline = at(0, boundary, parameters[0]);
    const auto& total = statistics_[(times_.size() - 1) * scenarios_];
    LikelihoodStatistics tail;
    for (size_t p = 0; p < 2; ++p) {
        tail.exposure[p] = total.exposure[p] - baseline.exposure[p];
        tail.coalescences[p] = total.coalescences[p] - baseline.coalescences[p];
    }
    double sum = -std::numeric_limits<double>::infinity();
    const double p = parameters[3];
    const double log_base = p < 1 ? std::log1p(-p) : 0;
    const double log_ghost = p > 0 ? std::log(p) : 0;
    const std::array<double, 2> inv = {0.5 / parameters[6], 0.5 / parameters[7]};
    const std::array<double, 2> log_inv = {std::log(inv[0]), std::log(inv[1])};
    auto value = [&](const LikelihoodStatistics& s) {
        double result = 0;
        for (size_t pop = 0; pop < 2; ++pop)
            result += s.coalescences[pop] * log_inv[pop] - s.exposure[pop] * inv[pop];
        return result;
    };
    for (size_t scenario = 0; scenario < scenarios_; ++scenario) {
        const auto s = at(scenario, boundary, parameters[0]);
        if ((p == 0 && s.lineages[1] != 0) || (p == 1 && s.lineages[0] != 0)) continue;
        double weight = value(s);
        if (p > 0 && p < 1)
            weight += s.lineages[0] * log_base + s.lineages[1] * log_ghost;
        sum = log_add(sum, weight);
    }
    return value(root_) + value(tail) + sum;
}

TreeEnsemble::TreeEnsemble(const tsk_treeseq_t& ts, std::span<const double> parameters,
    const std::vector<int>& sample_population, int admixed_index, int outgroup_index,
    bool defer_prepare, size_t tree_limit)
    : ts_(ts), parameters_(parameters.begin(), parameters.end()), sample_population_(sample_population),
      admixed_index_(admixed_index), outgroup_index_(outgroup_index),
      selected_count_(ts.num_trees) {
    validate_parameters(parameters);
    const auto* samples = tsk_treeseq_get_samples(&ts);
    for (tsk_size_t i = 0; i < tsk_treeseq_get_num_samples(&ts); ++i) {
        if (static_cast<size_t>(samples[i]) >= sample_population.size())
            throw std::invalid_argument("sample population missing for a sample node ID");
        if (sample_population[samples[i]] != admixed_index && sample_population[samples[i]] != outgroup_index)
            throw std::invalid_argument("cache supports admixed and outgroup sample populations only");
    }
    if (admixed_index == outgroup_index) throw std::invalid_argument("population indices must differ");
    if (tree_limit && tree_limit < size_t(ts.num_trees)) {
        selected_indices_.reserve(tree_limit);
        TreeGuard guard(ts);
        size_t probe = 0, index = 0;
        const long double length = tsk_treeseq_get_sequence_length(&ts);
        int ret;
        for (ret = tsk_tree_first(&guard.tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&guard.tree)) {
            while (probe < tree_limit) {
                const double position = static_cast<double>(
                    (static_cast<long double>(probe) + 0.5L) * length / tree_limit);
                if (position >= guard.tree.interval.right) break;
                if (position >= guard.tree.interval.left
                    && (selected_indices_.empty() || selected_indices_.back() != index))
                    selected_indices_.push_back(index);
                ++probe;
            }
            ++index;
        }
        check(ret);
        if (probe != tree_limit) throw std::runtime_error("failed to select trees at genomic positions");
        selected_count_ = selected_indices_.size();
    }
    if (!defer_prepare) prepare();
}

size_t TreeEnsemble::selected_index(size_t ordinal) const {
    return selected_indices_.empty() ? ordinal : selected_indices_[ordinal];
}

void TreeEnsemble::prepare(const CacheProgress& progress) {
    if (prepared_) return;
    if (progress) progress("prepare", 0, selected_count_);
    trees_.reserve(selected_count_);
    Scenario_Computer computer(ts_, parameters_, sample_population_, admixed_index_, outgroup_index_);
    TreeGuard guard(ts_);
    int ret;
    size_t done = 0, index = 0;
    for (ret = tsk_tree_first(&guard.tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&guard.tree)) {
        if (done < selected_count_ && index == selected_index(done)) {
            trees_.emplace_back();
            trees_.back().prepare(guard.tree, computer);
            if (progress) progress("prepare", ++done, selected_count_);
            else ++done;
        }
        ++index;
    }
    check(ret);
    prepared_ = true;
    memory_estimate(); // Fail on overflow before allocating scenario records.
}

CacheMemoryEstimate TreeEnsemble::memory_estimate() const {
    if (!prepared_) throw std::logic_error("prepare or load likelihood cache before estimating memory");
    CacheMemoryEstimate estimate;
    estimate.trees = trees_.size();
    estimate.resident_bytes = checked_add(sizeof(*this), checked_mul(trees_.capacity(), sizeof(TreeLikelihoodCache)));
    estimate.resident_bytes = checked_add(estimate.resident_bytes, checked_mul(parameters_.capacity(), sizeof(double)));
    estimate.resident_bytes = checked_add(estimate.resident_bytes, checked_mul(sample_population_.capacity(), sizeof(int)));
    for (const auto& tree : trees_) {
        estimate.scenarios = checked_add(estimate.scenarios, tree.scenario_count());
        estimate.records = checked_add(estimate.records, tree.record_count());
        estimate.resident_bytes = checked_add(estimate.resident_bytes, tree.metadata_bytes());
    }
    estimate.statistics_bytes = checked_mul(estimate.records, sizeof(LikelihoodStatistics));
    estimate.resident_bytes = checked_add(estimate.resident_bytes, estimate.statistics_bytes);
    return estimate;
}

void TreeEnsemble::build(const CacheProgress& progress) {
    if (ready_) return;
    prepare(progress);
    if (progress) progress("build", 0, selected_count_);
    Scenario_Computer computer(ts_, parameters_, sample_population_, admixed_index_, outgroup_index_);
    TreeGuard guard(ts_);
    int ret;
    size_t i = 0, index = 0;
    for (ret = tsk_tree_first(&guard.tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&guard.tree)) {
        if (i < selected_count_ && index == selected_index(i)) {
            trees_.at(i).build(guard.tree, computer);
            if (progress) progress("build", ++i, selected_count_);
            else ++i;
        }
        ++index;
    }
    check(ret);
    ready_ = true;
}

void TreeEnsemble::validate_query(std::span<const double> p) const {
    validate_parameters(p);
    if (!ready_) throw std::logic_error("build likelihood cache before querying it");
    if (p[1] != parameters_[1] || p[2] != parameters_[2])
        throw std::invalid_argument("split times changed: rebuild likelihood cache");
}
std::vector<double> TreeEnsemble::evaluate_trees(std::span<const double> parameters) {
    validate_query(parameters);
    std::vector<double> result;
    result.reserve(trees_.size());
    for (auto& tree : trees_) result.push_back(tree.evaluate(parameters));
    return result;
}
double TreeEnsemble::evaluate(std::span<const double> parameters) {
    validate_query(parameters);
    double result = 0;
    for (auto& tree : trees_) result += tree.evaluate(parameters);
    return result;
}

void TreeEnsemble::save(const std::string& path, const std::string& signature,
    const CacheProgress& progress) const {
    if (!ready_) throw std::logic_error("build likelihood cache before saving it");
    if (signature.size() != 16) throw std::invalid_argument("tree cache signature must have 16 hex characters");
    const std::string temporary = path + ".tmp." + std::to_string(
        std::chrono::steady_clock::now().time_since_epoch().count());
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot create tree statistics cache: " + temporary);
    CacheWriter writer{out};
    writer.bytes(cache_magic, sizeof(cache_magic));
    writer.bytes(signature.data(), signature.size());
    writer.pod(uint64_t(sizeof(LikelihoodStatistics)));
    writer.pod(uint64_t(sizeof(double)));
    writer.pod(uint64_t(sizeof(int)));
    writer.pod(uint64_t(alignof(LikelihoodStatistics)));
    writer.pod(uint64_t(offsetof(LikelihoodStatistics, exposure)));
    writer.pod(uint64_t(offsetof(LikelihoodStatistics, coalescences)));
    writer.pod(uint64_t(offsetof(LikelihoodStatistics, lineages)));
    writer.pod(uint64_t(0x0102030405060708ULL));
    writer.pod(uint64_t(trees_.size()));
    writer.pod(uint64_t(memory_estimate().resident_bytes));
    if (progress) progress("save", 0, trees_.size());
    size_t done = 0;
    for (const auto& tree : trees_) {
        writer.pod(uint64_t(tree.times_.size()));
        writer.pod(uint64_t(tree.scenarios_));
        writer.pod(uint64_t(tree.statistics_.size()));
        writer.pod(uint64_t(tree.impossible_));
        writer.pod(tree.left);
        writer.pod(tree.right);
        writer.pod(tree.split1_);
        writer.pod(tree.split2_);
        writer.pod(tree.root_);
        writer.bytes(tree.times_.data(), tree.times_.size() * sizeof(double));
        writer.bytes(tree.statistics_.data(), tree.statistics_.size() * sizeof(LikelihoodStatistics));
        if (progress && ++done < trees_.size()) progress("save", done, trees_.size());
    }
    out.write(reinterpret_cast<const char*>(&writer.hash), sizeof(writer.hash));
    out.close();
    if (!out) throw std::runtime_error("failed finishing tree statistics cache: " + temporary);
    std::filesystem::rename(temporary, path);
    if (progress) progress("save", trees_.size(), trees_.size());
}

void TreeEnsemble::load(const std::string& path, const std::string& signature, size_t max_bytes,
    const CacheProgress& progress) {
    if (prepared_ || ready_) throw std::logic_error("load tree cache into a fresh, unprepared ensemble");
    if (signature.size() != 16) throw std::invalid_argument("tree cache signature must have 16 hex characters");
    const uintmax_t file_bytes = std::filesystem::file_size(path);
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read tree statistics cache: " + path);
    CacheReader reader{in};
    char magic[8];
    reader.bytes(magic, sizeof(magic));
    char stored_signature[16];
    reader.bytes(stored_signature, sizeof(stored_signature));
    if (std::memcmp(magic, cache_magic, sizeof(magic)) != 0
        || std::memcmp(stored_signature, signature.data(), signature.size()) != 0)
        throw std::runtime_error("tree statistics cache signature/version mismatch: " + path);
    if (reader.pod<uint64_t>() != sizeof(LikelihoodStatistics)
        || reader.pod<uint64_t>() != sizeof(double)
        || reader.pod<uint64_t>() != sizeof(int)
        || reader.pod<uint64_t>() != alignof(LikelihoodStatistics)
        || reader.pod<uint64_t>() != offsetof(LikelihoodStatistics, exposure)
        || reader.pod<uint64_t>() != offsetof(LikelihoodStatistics, coalescences)
        || reader.pod<uint64_t>() != offsetof(LikelihoodStatistics, lineages)
        || reader.pod<uint64_t>() != 0x0102030405060708ULL)
        throw std::runtime_error("tree statistics cache ABI mismatch: " + path);
    if (reader.pod<uint64_t>() != selected_count_)
        throw std::runtime_error("tree statistics cache tree count mismatch: " + path);
    const uint64_t saved_estimate = reader.pod<uint64_t>();
    if (max_bytes && saved_estimate > max_bytes)
        throw std::runtime_error("estimated cache exceeds --cache-limit-mb");
    if (progress) progress("load", 0, selected_count_);
    trees_.reserve(selected_count_);
    for (size_t i = 0; i < selected_count_; ++i) {
        auto& tree = trees_.emplace_back();
        const uint64_t times = reader.pod<uint64_t>();
        const uint64_t scenarios = reader.pod<uint64_t>();
        const uint64_t records = reader.pod<uint64_t>();
        const uint64_t impossible = reader.pod<uint64_t>();
        if (times > file_bytes / sizeof(double) || scenarios > file_bytes / sizeof(LikelihoodStatistics)
            || records > file_bytes / sizeof(LikelihoodStatistics)
            || scenarios > std::numeric_limits<size_t>::max()
            || times > std::numeric_limits<size_t>::max()
            || records > std::numeric_limits<size_t>::max()
            || (times && records != checked_mul(size_t(times), size_t(scenarios)))
            || (!times && records != 0) || impossible > 1
            || (!impossible && (!times || !scenarios))
            || (impossible && (times || scenarios || records)))
            throw std::runtime_error("invalid tree statistics cache dimensions: " + path);
        tree.scenarios_ = size_t(scenarios);
        tree.impossible_ = impossible;
        tree.left = reader.pod<double>();
        tree.right = reader.pod<double>();
        tree.split1_ = reader.pod<double>();
        tree.split2_ = reader.pod<double>();
        tree.root_ = reader.pod<LikelihoodStatistics>();
        if (tree.split1_ != parameters_[2] || tree.split2_ != parameters_[1])
            throw std::runtime_error("tree statistics cache split time mismatch: " + path);
        tree.times_.resize(size_t(times));
        tree.statistics_.resize(size_t(records));
        reader.bytes(tree.times_.data(), tree.times_.size() * sizeof(double));
        reader.bytes(tree.statistics_.data(), tree.statistics_.size() * sizeof(LikelihoodStatistics));
        tree.ready_ = true;
        if (progress && i + 1 < selected_count_) progress("load", i + 1, selected_count_);
    }
    uint64_t stored_hash = 0;
    in.read(reinterpret_cast<char*>(&stored_hash), sizeof(stored_hash));
    if (!in || stored_hash != reader.hash || in.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("tree statistics cache checksum/trailing data mismatch: " + path);
    prepared_ = true;
    ready_ = true;
    memory_estimate();
    if (progress) progress("load", selected_count_, selected_count_);
}
}
