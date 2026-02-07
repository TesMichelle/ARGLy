#include <vector>
#include <span>
#include <map>
#include <cmath>
#include <algorithm>
#include <tskit.h>

#include <iostream>

namespace treellh
{
    class CoalTable
    {
        private:
            std::vector<double> time_;
            int pit_;
        public:
            CoalTable();
            CoalTable(const std::vector<double>&);
            CoalTable(const int&, const int&);

            int size() const;

            double get_likelyhood() const;

            CoalTable operator+(const CoalTable &B);
            CoalTable operator-(const CoalTable &B);
    };

    class Scenario_Computer
    {
        private:
            const tsk_treeseq_t &ts_;

            std::vector<bool> is_fixed_;
            std::vector<int> population_;
            std::vector<tsk_id_t> lower_nodes_split_1_;
            std::vector<tsk_id_t> time_ordered_tree_nodes_;

            int lower_split_2_out_linages_num_;

            int node_below_split_1_index_;

            int log_scenarios_num_;
            int sc_lineages_num_[2];
            int lineages_num_lower_migration_;

            double migration_time_;
            double migration_prob_;
            double split_1_time_;
            double split_2_time_;
            double split_1_prop_;
            double split_2_prop_;

            std::vector<double> N_;

            int outgroup_data_index_;
            int admixed_data_index_;
            int base_pop_index_ = 0;
            int ghost_pop_index_ = 1;
            int outgroup_pop_index_ = 2;
        public:
            Scenario_Computer(const tsk_treeseq_t&, std::span<const double>, int, int);

            int compute_lower_nodes_split_1(const tsk_tree_t&);
            int compute_fixed_nodes(const tsk_tree_t&);
            int compute_log_scenarios_num();

            int compute_scenarious();

            std::vector<uint64_t> generate_abacaba(uint64_t);

            int set_subtree_population(const tsk_tree_t&, tsk_id_t, int, double);
            int set_outgroup(const tsk_tree_t&, tsk_id_t*, tsk_size_t);
            int set_scenario(const tsk_tree_t&, int);
            int check_scenario(const tsk_tree_t&);
            int check_imposible_split_2(const tsk_tree_t&);

            int sort_by_time(const tsk_tree_t&);
            int debug_sort_by_time(const tsk_tree_t&);

            double smart_llh(const tsk_tree_t &, double, double, 
                std::span<const int>, std::span<int>);

            double compute_llh(const tsk_tree_t&);

            void set_parameters(std::span<const double>);
            std::vector<tsk_id_t> get_lower_nodes_split_1();

            void print_lower_nodes_split_1();
    };

    std::vector<int> lower_node(const tsk_tree_t&, double split_time);
}