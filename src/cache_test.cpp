#include <likelihood_cache.h>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void check(int code) { if (code < 0) throw std::runtime_error(tsk_strerror(code)); }
void close(double a, double b) {
    if (a == b) return;
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > 1e-9) {
        std::cerr << "cached=" << a << " reference=" << b << '\n';
        throw std::runtime_error("likelihood mismatch");
    }
}
template<class F> void rejects(F f) {
    try { f(); } catch (const std::exception&) { return; }
    throw std::runtime_error("expected an exception");
}
double logsum(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) return b;
    if (b == -std::numeric_limits<double>::infinity()) return a;
    return std::max(a, b) + std::log1p(std::exp(-std::abs(a - b)));
}
}

int main() {
    try {
        tsk_table_collection_t tables{};
        check(tsk_table_collection_init(&tables, 0));
        tables.sequence_length = 3;
        for (int i = 0; i < 6; ++i)
            check(tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, TSK_NULL, nullptr, 0));
        for (double time : {300., 700., 400., 2200., 3000., 1300., 1500.})
            check(tsk_node_table_add_row(&tables.nodes, 0, time, TSK_NULL, TSK_NULL, nullptr, 0));
        // Shared node IDs, different parents/children across genomic intervals.
        std::vector<std::vector<std::pair<int,int>>> edges = {
            {{6,0},{6,1},{7,2},{7,3},{8,4},{8,5},{9,6},{9,7},{10,9},{10,8}},
            {{6,0},{6,2},{7,1},{7,3},{8,4},{8,5},{11,6},{11,8},{10,11},{10,7}},
            {{6,0},{6,1},{7,2},{7,3},{8,4},{8,5},{12,6},{12,7},{10,12},{10,8}}
        };
        for (size_t i = 0; i < edges.size(); ++i)
            for (auto [parent, child] : edges[i])
                check(tsk_edge_table_add_row(&tables.edges, i, i+1, parent, child, nullptr, 0));
        check(tsk_table_collection_sort(&tables, nullptr, 0));
        tsk_treeseq_t ts{};
        check(tsk_treeseq_init(&ts, &tables, TSK_TS_INIT_BUILD_INDEXES));
        std::vector<double> p = {500,1000,2000,0.2,0.5,0.5,10000,12000,9000};
        std::vector<int> pops = {0,0,0,0,1,1};
        treellh::TreeEnsemble ensemble(ts, p, pops, 0, 1);
        const auto estimate = ensemble.memory_estimate();
        require(estimate.trees == 3 && estimate.statistics_bytes > 0, "bad memory estimate");
        require(estimate.statistics_bytes == estimate.records * sizeof(treellh::LikelihoodStatistics), "bad record size");
        rejects([&] { ensemble.evaluate(p); });
        ensemble.build();
        ensemble.build(); // Idempotent.
        size_t comparisons = 0;
        double max_grid_error = 0;
        // Reuse the SAME cache: local moves, both neighbouring intervals,
        // large jumps, exact event times, endpoint probabilities, changing N.
        for (double n : {10000., 120., 30000.}) {
            for (double probability : {0.2, 0., 1., 0.8}) {
                for (double time : {500.,501.,701.,699.,100.,900.,700.,300.,0.,999.,500.}) {
                    p[0] = time; p[3] = probability; p[6] = n; p[7] = 1.7*n;
                    const auto cached = ensemble.evaluate_trees(p);
                    // Independent formula for tree 0 at p=0: four included
                    // coalescences, pair-time exposure 7700; outgroup omitted.
                    if (probability == 0)
                        close(cached[0], 4 * std::log(0.5 / n) - 7700 * (0.5 / n));
                    double sum = 0;
                    tsk_tree_t tree{};
                    check(tsk_tree_init(&tree, &ts, 0));
                    treellh::Scenario_Computer direct(ts, p, pops, 0, 1);
                    size_t i = 0;
                    int ret;
                    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
                        const auto parts = direct.compute_llh_ratio(tree, p, 0);
                        if (std::abs(cached[i] - logsum(parts.first, parts.second)) > 1e-9)
                            std::cerr << "direct mismatch tree=" << i << " t=" << time << " p=" << probability << " N=" << n << '\n';
                        close(cached[i], logsum(parts.first, parts.second));
                        // Grid is independently compared at non-event times.
                        if (time != 0 && time != 300 && time != 700) {
                            direct.set_parameters(p);
                            const auto grid = direct.compute_grid_fast(tree, {time}, {probability});
                            const auto ratio_grid = direct.compute_grid_ratio_fast(
                                tree, {time}, {probability}, {}, 0);
                            close(ratio_grid.first(0,0,0), parts.first);
                            close(ratio_grid.second(0,0,0), parts.second);
                            close(grid(0,0,0), logsum(ratio_grid.first(0,0,0), ratio_grid.second(0,0,0)));
                            if (std::abs(cached[i] - grid(0,0,0)) > 1e-9)
                                std::cerr << "grid mismatch tree=" << i << " t=" << time << " p=" << probability << " N=" << n << '\n';
                            close(cached[i], grid(0,0,0));
                            max_grid_error = std::max(max_grid_error, std::abs(cached[i] - grid(0,0,0)));
                        }
                        sum += cached[i++];
                        ++comparisons;
                    }
                    check(ret);
                    tsk_tree_free(&tree);
                    close(sum, ensemble.evaluate(p));
                }
            }
        }
        {
            tsk_tree_t tree{};
            check(tsk_tree_init(&tree, &ts, 0));
            check(tsk_tree_first(&tree));
            std::vector<double> reference = {500,1000,2000,0.2,0.5,0.5,10000,12000,9000};
            treellh::Scenario_Computer direct(ts, reference, pops, 0, 1);
            const std::vector<double> times = {900,501,250};
            const std::vector<double> probabilities = {0,0.2,1};
            const std::vector<std::vector<double>> sizes = {{10000,20000},{12000,24000},{9000}};
            const auto all_grid = direct.compute_grid_fast(tree, times, probabilities, sizes);
            for (tsk_id_t sample : {tsk_id_t{0}, tsk_id_t{2}}) {
                auto [base, ghost] = direct.compute_grid_ratio_fast(
                    tree, times, probabilities, sizes, sample);
                for (size_t t = 0; t < times.size(); ++t)
                    for (size_t n = 0; n < 2; ++n)
                        for (size_t q = 0; q < probabilities.size(); ++q) {
                            reference[0] = times[t]; reference[3] = probabilities[q];
                            reference[6] = sizes[0][n]; reference[7] = sizes[1][n];
                            const auto scalar = direct.compute_llh_ratio(tree, reference, sample);
                            close(base(t,n,q), scalar.first);
                            close(ghost(t,n,q), scalar.second);
                            close(logsum(base(t,n,q), ghost(t,n,q)), all_grid(t,n,q));
                        }
            }
            tsk_tree_free(&tree);
        }
        p[1] += 1;
        rejects([&] { ensemble.evaluate(p); });
        p[1] -= 1;
        p[6] = 0;
        rejects([&] { ensemble.evaluate(p); });
        p[6] = 10000;
        const auto cache_path = std::filesystem::temp_directory_path()
            / ("argly_tree_cache_" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count()) + ".bin");
        const std::string signature = "0123456789abcdef";
        ensemble.save(cache_path.string(), signature);
        treellh::TreeEnsemble restored(ts, p, pops, 0, 1, true);
        restored.load(cache_path.string(), signature);
        for (double time : {100., 500., 900.}) {
            p[0] = time;
            for (double probability : {0., 0.2, 1.}) {
                p[3] = probability;
                p[6] = 15000;
                p[7] = 8000;
                const auto before = ensemble.evaluate_trees(p);
                const auto after = restored.evaluate_trees(p);
                require(before.size() == after.size(), "loaded cache tree count mismatch");
                for (size_t i = 0; i < before.size(); ++i) close(before[i], after[i]);
            }
        }
        rejects([&] {
            treellh::TreeEnsemble mismatched(ts, p, pops, 0, 1, true);
            mismatched.load(cache_path.string(), "fedcba9876543210");
        });
        rejects([&] {
            treellh::TreeEnsemble limited(ts, p, pops, 0, 1, true);
            limited.load(cache_path.string(), signature, 64);
        });
        {
            std::fstream corrupted(cache_path, std::ios::binary | std::ios::in | std::ios::out);
            corrupted.seekg(-static_cast<std::streamoff>(sizeof(uint64_t) + 1), std::ios::end);
            char byte = 0;
            corrupted.read(&byte, 1);
            corrupted.seekp(-static_cast<std::streamoff>(sizeof(uint64_t) + 1), std::ios::end);
            byte ^= 1;
            corrupted.write(&byte, 1);
        }
        rejects([&] {
            treellh::TreeEnsemble corrupted(ts, p, pops, 0, 1, true);
            corrupted.load(cache_path.string(), signature);
        });
        std::filesystem::remove(cache_path);
        auto bad_pops = pops;
        bad_pops[1] = 1; // Admixed/outgroup coalescence at 300, before split2.
        treellh::TreeEnsemble impossible(ts, p, bad_pops, 0, 1);
        impossible.build();
        require(impossible.evaluate_trees(p)[0] == -std::numeric_limits<double>::infinity(),
            "impossible tree must have zero likelihood");
        tsk_treeseq_free(&ts);
        tsk_table_collection_free(&tables);
        std::cout << "Passed " << comparisons << " cached/direct comparisons; max grid error=" << max_grid_error << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
