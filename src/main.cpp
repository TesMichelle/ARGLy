extern "C" {
#include <tskit.h>
}

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <string>
#include <vector>
#include <omp.h>

#ifdef TIME
#include <chrono>
#endif

// ============================================================================
// Параметры популяции
// ============================================================================
static double N_AFR    = 10000.0;
static double N_SAMPLE = 10000.0;
static double N_ANC    = 10000.0;
static double N_GHOST  = 10000.0;
static double N_ARCH   = 10000.0;

static double T_MIG   = 10000.0;
static double T_SEP   = 30000.0;
static double T_INTRO  = 50000.0;

static const double ADMIX_RATE = 0.2;

// ============================================================================
// Предвычисленные константы (пересчитываются при изменении N_GHOST)
// ============================================================================
static double LOG_ADMIX;
static double LOG_1_ADMIX;
static double LOG_INV_2N_AFR;
static double LOG_INV_2N_ANC;
static double LOG_INV_2N_GHOST;
static double LOG_INV_2N_SAMPLE;
static double LOG_INV_2N_ARCH;
static double INV_4N_AFR;
static double INV_4N_ANC;
static double INV_4N_GHOST;
static double INV_4N_SAMPLE;
static double INV_4N_ARCH;

static void precompute_constants() {
    LOG_ADMIX       = std::log(ADMIX_RATE);
    LOG_1_ADMIX     = std::log(1.0 - ADMIX_RATE);
    LOG_INV_2N_AFR  = std::log(1.0 / (2.0 * N_AFR));
    LOG_INV_2N_ANC  = std::log(1.0 / (2.0 * N_ANC));
    LOG_INV_2N_GHOST= std::log(1.0 / (2.0 * N_GHOST));
    LOG_INV_2N_SAMPLE = std::log(1.0 / (2.0 * N_SAMPLE));
    LOG_INV_2N_ARCH = std::log(1.0 / (2.0 * N_ARCH));
    INV_4N_AFR      = 1.0 / (4.0 * N_AFR);
    INV_4N_ANC      = 1.0 / (4.0 * N_ANC);
    INV_4N_GHOST    = 1.0 / (4.0 * N_GHOST);
    INV_4N_SAMPLE   = 1.0 / (4.0 * N_SAMPLE);
    INV_4N_ARCH     = 1.0 / (4.0 * N_ARCH);
}

// ============================================================================
// cast: определение bf_nodes (линии, живые в T_INTRO и не ведущие к SAM)
// ============================================================================
static bool cast(tsk_id_t start, tsk_tree_t *tree,
                 std::vector<tsk_id_t> &bf_nodes) {
    const tsk_node_table_t &nodes = tree->tree_sequence->tables->nodes;

    if (nodes.population[start] == 1 && nodes.time[start] <= T_SEP) {
        return true;
    }
    if (tree->left_child[start] == TSK_NULL) {
        return false;
    }

    tsk_id_t v = tree->left_child[start];
    bool ok = false;

    while (v != TSK_NULL) {
        bool is_sam_anc = cast(v, tree, bf_nodes);

        if (nodes.time[tree->parent[v]] >= T_INTRO && nodes.time[v] < T_INTRO) {
            if (!is_sam_anc) {
                bf_nodes.emplace_back(v);
            }
        }
        if (is_sam_anc) {
            ok = true;
        }
        v = tree->right_sib[v];
    }
    return ok;
}

// ============================================================================
// calculate_p_mig_fix: вклад coalescence SAM + AFR до миграции
// ============================================================================
static void calculate_p_mig_fix(tsk_tree_t *tree,
                                std::vector<tsk_id_t> &sam,
                                std::vector<tsk_id_t> &afr,
                                std::vector<int> &k,
                                double &p_mig_fix) {
    const tsk_node_table_t &node = tree->tree_sequence->tables->nodes;
    double t = 0.0;

    for (tsk_id_t u : sam) {
        if (node.time[u] == 0.0 || k[1] <= 1) continue;
        p_mig_fix += LOG_INV_2N_SAMPLE -
                     (static_cast<double>(k[1]) * (k[1] - 1) * INV_4N_SAMPLE * (node.time[u] - t));
        t = node.time[u];
        --k[1];
    }
    p_mig_fix -= static_cast<double>(k[1]) * (k[1] - 1) * INV_4N_SAMPLE * (T_SEP - t);

    t = 0.0;
    for (tsk_id_t u : afr) {
        if (node.time[u] == 0.0 || k[0] <= 1) continue;
        p_mig_fix += LOG_INV_2N_AFR -
                     (static_cast<double>(k[0]) * (k[0] - 1) * INV_4N_AFR * (node.time[u] - t));
        t = node.time[u];
        --k[0];
    }
    p_mig_fix -= static_cast<double>(k[0]) * (k[0] - 1) * INV_4N_AFR * (T_MIG - t);
}

// ============================================================================
// calculate_p_intro_fix: вклад coalescence ARCHAIC после интрогрессии
// ============================================================================
static void calculate_p_intro_fix(tsk_tree_t *tree,
                                  std::vector<tsk_id_t> &arch,
                                  std::vector<int> &k,
                                  double &p_intro_fix) {
    const tsk_node_table_t &node = tree->tree_sequence->tables->nodes;
    double t = T_INTRO;

    for (tsk_id_t u : arch) {
        if (k[5] <= 1) continue;
        p_intro_fix += LOG_INV_2N_ARCH -
                       (static_cast<double>(k[5]) * (k[5] - 1) * INV_4N_ARCH * (node.time[u] - t));
        t = node.time[u];
        --k[5];
    }
}

// ============================================================================
// Struct: данные одного дерева для параллельной обработки (SoA)
// ============================================================================
struct TreeTask {
    // Phase 1: узлы в (T_MIG, T_SEP] — AFR "stay" coalescence + ghost coalescence
    // Phase 2: узлы в (T_SEP, T_INTRO] — ancestor coalescence + ghost coalescence
    // Для каждого узла: время и bf_id (-1 = fixed/not-in-bf, >=0 = active bf index)
    std::vector<double> phase1_times;
    std::vector<int>    phase1_bf_id;    // -1 = fixed (not in any bf subtree)

    std::vector<double> phase2_times;
    std::vector<int>    phase2_bf_id;

    // Предвычисленный вклад fixed nodes в phase1
    double fixed_p1;           // log-probability вклад fixed phase1 nodes
    int    fixed_k2_decrease;  // сколько coalescence events от fixed nodes в phase1

    // Active bf_nodes: delta_k3 для каждого
    std::vector<int> active_delta_k3;

    // DP-результат от inactive bf_nodes: (k3_offset, count)
    std::vector<int>    dp_k3_offsets;
    std::vector<double> dp_counts;

    // Начальные счётчики
    int k0_init;  // AFR lineages at time 0
    int k1_init;  // SAM lineages at time 0 (after p_mig_fix)
    int k3_init;  // ghost lineages already counted (should be 0 initially)

    // Base log-probability
    double p_base;

    // Число active bf_nodes
    int M_active;
};

int main(int argc, char *argv[]) {
    if (argc >= 5) {
        T_MIG   = std::atof(argv[1]);
        N_GHOST = std::atof(argv[2]);
        T_SEP   = std::atof(argv[3]);
        T_INTRO = std::atof(argv[4]);
    }

    precompute_constants();

    std::string filename = "my_arg(2).arg";
    double p_res = 0;
    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, filename.c_str(), 0);
    if (ret != 0) {
        std::cerr << "Error loading tree sequence: " << tsk_strerror(ret) << " (" << ret << ")\n";
        return 1;
    }

    tsk_tree_t tree;
    ret = tsk_tree_init(&tree, &ts, 0);
    if (ret != 0) {
        std::cerr << "Error initializing tree: " << tsk_strerror(ret) << " (" << ret << ")\n";
        return 1;
    }

    std::vector<TreeTask> tasks;

#ifdef TIME
    auto start_time = std::chrono::high_resolution_clock::now();
#endif

    for (ret = tsk_tree_first(&tree); ret == 1; ret = tsk_tree_next(&tree)) {
        tsk_id_t root = tsk_tree_get_left_root(&tree);
        const tsk_node_table_t &node = tree.tree_sequence->tables->nodes;

                double p_mig_fix = 0.0;
                double p_intro_fix = 0.0;

                std::vector<tsk_id_t> bf_nodes;
                cast(root, &tree, bf_nodes);

                // 0=afr, 1=sam, 2=stay, 3=ghost, 4=ancestors, 5=archaic
                std::vector<int> k = {0, 0, 0, 0, 0, 0};
                std::vector<tsk_id_t> nodes, arch, sam, afr_before_mig;

                const tsk_size_t num_nodes = tree.num_nodes;
                for (tsk_id_t i = 0; i < static_cast<tsk_id_t>(num_nodes); ++i) {
                    const double t_i = node.time[i];
                    const tsk_id_t pop_i = node.population[i];

                    if (t_i == 0.0) {
                        if (pop_i == 1) ++k[1];
                        else ++k[0];
                    }
                    if (pop_i == 1 && t_i <= T_SEP) {
                        sam.push_back(i);
                    }
                    if (pop_i != 1 && t_i <= T_MIG) {
                        afr_before_mig.push_back(i);
                    }
                    if (t_i > T_SEP && t_i <= T_INTRO) {
                        nodes.push_back(i);
                    }
                    if (t_i > T_MIG && t_i <= T_SEP && pop_i != 1) {
                        nodes.push_back(i);
                    }
                    if (t_i >= T_INTRO) {
                        arch.push_back(i);
                    }
                    if (t_i > T_INTRO) {
                        if (tree.left_child[i] != TSK_NULL &&
                            node.time[tree.left_child[i]] <= T_INTRO) ++k[5];
                        if (tree.right_child[i] != TSK_NULL &&
                            node.time[tree.right_child[i]] <= T_INTRO) ++k[5];
                    }
                }

                auto cmp = [&](tsk_id_t a, tsk_id_t b) {
                    return node.time[a] < node.time[b];
                };
                std::sort(sam.begin(), sam.end(), cmp);
                std::sort(afr_before_mig.begin(), afr_before_mig.end(), cmp);
                std::sort(arch.begin(), arch.end(), cmp);
                std::sort(nodes.begin(), nodes.end(), cmp);

                calculate_p_mig_fix(&tree, sam, afr_before_mig, k, p_mig_fix);
                calculate_p_intro_fix(&tree, arch, k, p_intro_fix);

                const tsk_size_t total_nodes = tree.tree_sequence->tables->nodes.num_rows;
                std::vector<int> belongs_to_bf(total_nodes, -1);
                std::vector<int> delta_k3_all(bf_nodes.size(), 0);

                for (size_t i = 0; i < bf_nodes.size(); ++i) {
                    std::queue<tsk_id_t> q;
                    q.push(bf_nodes[i]);
                    while (!q.empty()) {
                        tsk_id_t u = q.front(); q.pop();
                        belongs_to_bf[u] = static_cast<int>(i);
                        if (node.time[u] <= T_MIG &&
                            node.time[tree.parent[u]] > T_MIG) {
                            delta_k3_all[i]++;
                        }
                        tsk_id_t v = tree.left_child[u];
                        while (v != TSK_NULL) {
                            q.push(v);
                            v = tree.right_sib[v];
                        }
                    }
                }

                std::vector<bool> is_active(bf_nodes.size(), false);
                for (tsk_id_t u : nodes) {
                    if (belongs_to_bf[u] != -1) {
                        is_active[belongs_to_bf[u]] = true;
                    }
                }

                std::vector<int> active_indices;
                std::vector<int> inactive_indices;
                for (size_t i = 0; i < bf_nodes.size(); ++i) {
                    if (is_active[i]) active_indices.push_back(static_cast<int>(i));
                    else inactive_indices.push_back(static_cast<int>(i));
                }

                // --- DP по inactive (subset-sum для k3_offset) ---
                int max_k3_inactive = 0;
                for (int idx : inactive_indices) max_k3_inactive += delta_k3_all[idx];

                std::vector<double> DP(max_k3_inactive + 1, 0.0);
                DP[0] = 1.0;
                for (int idx : inactive_indices) {
                    int d = delta_k3_all[idx];
                    for (int mk = max_k3_inactive; mk >= d; --mk) {
                        DP[mk] += DP[mk - d];
                    }
                }

                TreeTask task;
                task.k0_init = k[0];
                task.k1_init = k[1];
                task.k3_init = k[3];
                task.p_base  = p_mig_fix + p_intro_fix;
                task.M_active = static_cast<int>(active_indices.size());

                for (int idx : active_indices)
                    task.active_delta_k3.push_back(delta_k3_all[idx]);

                for (int mk = 0; mk <= max_k3_inactive; ++mk) {
                    if (DP[mk] > 0.0) {
                        task.dp_k3_offsets.push_back(mk);
                        task.dp_counts.push_back(DP[mk]);
                    }
                }

                std::vector<int> orig_to_active(bf_nodes.size(), -1);
                for (int j = 0; j < static_cast<int>(active_indices.size()); ++j) {
                    orig_to_active[active_indices[j]] = j;
                }

                std::vector<double> fixed_p1_times;
                std::vector<double> ghost_p1_times;
                std::vector<int>    ghost_p1_bf_id;

                for (tsk_id_t u : nodes) {
                    if (node.time[u] > T_SEP) continue;
                    int orig_bf = belongs_to_bf[u];
                    if (orig_bf == -1) {
                        fixed_p1_times.push_back(node.time[u]);
                    } else {
                        int act_idx = orig_to_active[orig_bf];
                        if (act_idx != -1) {
                            ghost_p1_times.push_back(node.time[u]);
                            ghost_p1_bf_id.push_back(act_idx);
                        }
                    }
                }

                for (tsk_id_t u : nodes) {
                    if (node.time[u] <= T_SEP) continue;
                    int orig_bf = belongs_to_bf[u];
                    int act_idx = -1;
                    if (orig_bf != -1) {
                        act_idx = orig_to_active[orig_bf];
                    }
                    task.phase2_times.push_back(node.time[u]);
                    task.phase2_bf_id.push_back(act_idx);
                }

                task.phase1_times = std::move(ghost_p1_times);
                task.phase1_bf_id = std::move(ghost_p1_bf_id);

                task.phase1_times.clear();
                task.phase1_bf_id.clear();
                for (tsk_id_t u : nodes) {
                    if (node.time[u] > T_SEP) continue;
                    int orig_bf = belongs_to_bf[u];
                    int act_idx = -1;
                    if (orig_bf != -1) {
                        act_idx = orig_to_active[orig_bf];
                    }
                    task.phase1_times.push_back(node.time[u]);
                    task.phase1_bf_id.push_back(act_idx);
                }

                task.fixed_p1 = 0.0;
                task.fixed_k2_decrease = 0;

                tasks.push_back(std::move(task));
    }

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);

    double final_p = 0.0;

    #pragma omp parallel for reduction(+:final_p) schedule(dynamic)
    for (size_t t = 0; t < tasks.size(); ++t) {
        const TreeTask &task = tasks[t];
        const int M_act = task.M_active;
        const uint64_t max_mask = (1ULL << M_act);
        const size_t n_dp = task.dp_k3_offsets.size();
        const size_t n_p1 = task.phase1_times.size();
        const size_t n_p2 = task.phase2_times.size();



        std::vector<double> log_probs;
        log_probs.reserve(max_mask * n_dp);

        for (uint64_t mask = 0; mask < max_mask; ++mask) {
            int active_k3 = 0;
            for (int b = 0; b < M_act; ++b) {
                if (mask & (1ULL << b)) {
                    active_k3 += task.active_delta_k3[b];
                }
            }

            for (size_t di = 0; di < n_dp; ++di) {
                const int mk = task.dp_k3_offsets[di];
                const double prob_weight = task.dp_counts[di];

                const int total_k3_offset = active_k3 + mk;

                double p = total_k3_offset * LOG_ADMIX +
                        (task.k0_init - total_k3_offset) * LOG_1_ADMIX;

                double t_ghost = T_MIG;
                double t_stay  = T_MIG;

                double k2 = static_cast<double>(task.k0_init - (task.k3_init + total_k3_offset));
                double k3 = static_cast<double>(task.k3_init + total_k3_offset);

                for (size_t i = 0; i < n_p1; ++i) {
                    const int bf_id = task.phase1_bf_id[i];
                    const double ntime = task.phase1_times[i];
                    if (bf_id >= 0 && (mask & (1ULL << bf_id))) {
                        if (k3 > 1.0) {
                            p += LOG_INV_2N_GHOST -
                                (k3 * (k3 - 1.0) * INV_4N_GHOST * (ntime - t_ghost));
                            k3 -= 1.0;
                            t_ghost = ntime;
                        }
                    } else {
                        if (k2 > 1.0) {
                            p += LOG_INV_2N_AFR -
                                (k2 * (k2 - 1.0) * INV_4N_AFR * (ntime - t_stay));
                            k2 -= 1.0;
                            t_stay = ntime;
                        }
                    }
                }

                if (k2 > 1.0) {
                    p -= k2 * (k2 - 1.0) * INV_4N_AFR * (T_SEP - t_stay);
                }
                double k4 = task.k1_init + k2;
                double t_anc = T_SEP;

                for (size_t i = 0; i < n_p2; ++i) {
                    const int bf_id = task.phase2_bf_id[i];
                    const double ntime = task.phase2_times[i];

                    if (bf_id >= 0 && (mask & (1ULL << bf_id))) {
                        if (k3 > 1.0) {
                            p += LOG_INV_2N_GHOST -
                                (k3 * (k3 - 1.0) * INV_4N_GHOST * (ntime - t_ghost));
                            k3 -= 1.0;
                            t_ghost = ntime;
                        }
                    } else {
                        if (k4 > 1.0) {
                            p += LOG_INV_2N_ANC -
                                (k4 * (k4 - 1.0) * INV_4N_ANC * (ntime - t_anc));
                            k4 -= 1.0;
                            t_anc = ntime;
                        }
                    }
                }

                if (k4 > 1.0)
                    p -= k4 * (k4 - 1.0) * INV_4N_ANC * (T_INTRO - t_anc);
                if (k3 > 1.0)
                    p -= k3 * (k3 - 1.0) * INV_4N_GHOST * (T_INTRO - t_ghost);
                p += std::log(prob_weight);

                log_probs.push_back(p);
            }
        }

        if (log_probs.empty()) {
            continue;
        }

        double max_lp = log_probs[0];
        for (double lp : log_probs) {
            if (lp > max_lp) max_lp = lp;
        }

        double sum_exp = 0.0;
        for (double lp : log_probs) {
            sum_exp += std::exp(lp - max_lp);
        }

        final_p += task.p_base + max_lp + std::log(sum_exp);
    }

    #ifdef TIME
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time);
        // std::cerr << "Execution time: " << duration.count() << " ms" << std::endl;
    #endif

    // std::cout << final_p << '\n';
    p_res = final_p;
    std::cout << p_res << '\n';
    return 0;
}
