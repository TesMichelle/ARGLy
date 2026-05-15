#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <mdspan>
#include <tskit.h>

#include <iostream>

namespace treellh
{
    using DoubleMatrixView = std::mdspan<double, std::dextents<size_t, 2>>;
    using IntMatrixView = std::mdspan<int, std::dextents<size_t, 2>>;
    using DoubleTensor3DView = std::mdspan<double, std::dextents<size_t, 3>>;

    class DenseMatrix
    {
        private:
            size_t rows_ = 0;
            size_t cols_ = 0;
            std::vector<double> data_;
        public:
            DenseMatrix() = default;
            DenseMatrix(size_t rows, size_t cols, double value = 0)
                : rows_(rows), cols_(cols), data_(rows * cols, value) {}

            size_t rows() const { return rows_; }
            size_t cols() const { return cols_; }
            size_t size() const { return data_.size(); }
            const std::vector<double>& data() const { return data_; }
            std::vector<double>& data() { return data_; }

            DoubleMatrixView view() { return DoubleMatrixView(data_.data(), rows_, cols_); }
            std::mdspan<const double, std::dextents<size_t, 2>> view() const
            {
                return std::mdspan<const double, std::dextents<size_t, 2>>(data_.data(), rows_, cols_);
            }

            double& operator()(size_t row, size_t col) { return data_[row * cols_ + col]; }
            double operator()(size_t row, size_t col) const { return data_[row * cols_ + col]; }
    };

    class DenseTensor3D
    {
        private:
            size_t dim0_ = 0;
            size_t dim1_ = 0;
            size_t dim2_ = 0;
            std::vector<double> data_;
        public:
            DenseTensor3D() = default;
            DenseTensor3D(size_t dim0, size_t dim1, size_t dim2, double value = 0)
                : dim0_(dim0), dim1_(dim1), dim2_(dim2), data_(dim0 * dim1 * dim2, value) {}

            size_t dim0() const { return dim0_; }
            size_t dim1() const { return dim1_; }
            size_t dim2() const { return dim2_; }
            size_t size() const { return data_.size(); }
            const std::vector<double>& data() const { return data_; }
            std::vector<double>& data() { return data_; }

            DoubleTensor3DView view() { return DoubleTensor3DView(data_.data(), dim0_, dim1_, dim2_); }
            std::mdspan<const double, std::dextents<size_t, 3>> view() const
            {
                return std::mdspan<const double, std::dextents<size_t, 3>>(data_.data(), dim0_, dim1_, dim2_);
            }

            double& operator()(size_t i, size_t j, size_t k)
            {
                return data_[(i * dim1_ + j) * dim2_ + k];
            }
            double operator()(size_t i, size_t j, size_t k) const
            {
                return data_[(i * dim1_ + j) * dim2_ + k];
            }
    };

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

    struct Node
    {
        int population = -1;
        double time = -1;
        bool is_fixed = false;
        bool is_skipped = false;
    };

    class Scenario_Computer
    {
        private:
            const tsk_treeseq_t &ts_;
            std::vector<int> sample_population_;

            std::vector<Node> nodes_;

            std::vector<bool> is_fixed_;
            std::vector<int> population_;
            std::vector<tsk_id_t> lower_nodes_split_1_;
            std::vector<tsk_id_t> time_ordered_tree_nodes_ids_;

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
            Scenario_Computer(const tsk_treeseq_t&, std::span<const double>, const std::vector<int>&, int, int);

            int compute_lower_nodes_split_1(const tsk_tree_t&);
            int compute_fixed_nodes(const tsk_tree_t&);
            int compute_skipped_nodes(const tsk_tree_t&);
            int compute_log_scenarios_num();

            int compute_scenarious();

            std::vector<uint64_t> generate_abacaba(uint64_t);

            int set_subtree_population(const tsk_tree_t&, tsk_id_t, int, double);
            int set_outgroup(const tsk_tree_t&, tsk_id_t*, tsk_size_t);
            int set_scenario(const tsk_tree_t&, uint64_t);
            int check_scenario(const tsk_tree_t&);
            int check_imposible_split_2(const tsk_tree_t&);

            double find_lowest_coal(const tsk_tree_t&, std::span<const tsk_id_t>, std::span<const tsk_id_t>);

            int sort_by_time(const tsk_tree_t&);
            int debug_sort_by_time(const tsk_tree_t&);

            double smart_llh(double, double, 
                std::span<const int>, std::span<int>);
            double smart_llh_to_nodes(double, double, 
                std::span<const int>, std::span<int>, DoubleMatrixView, IntMatrixView,
                std::mdspan<const double, std::dextents<size_t, 2>>,
                std::mdspan<const double, std::dextents<size_t, 2>>,
                std::span<double>, size_t);

            DenseTensor3D compute_grid_fast(
                const tsk_tree_t &, std::vector<double>, std::vector<double>,
                std::vector<std::vector<double>> = {});

            double compute_llh(const tsk_tree_t&);

            void set_parameters(std::span<const double>);
            std::vector<tsk_id_t> get_lower_nodes_split_1();

            void print_lower_nodes_split_1();
    };

    std::vector<int> lower_node(const tsk_tree_t&, double split_time);
}
