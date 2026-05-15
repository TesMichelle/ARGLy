#include "treellh.h"

#include <limits>
#include <stdexcept>

namespace treellh
{
    namespace
    {
        double grid_N_value(
            const std::vector<std::vector<double>>& N_grid,
            int population,
            size_t n_i
        )
        {
            const auto &values = N_grid.at(static_cast<size_t>(population));
            return values.size() == 1 ? values[0] : values[n_i];
        }
    }

    // CoalTable implemantation

    CoalTable::CoalTable()
    {
        time_ = std::vector<double>();
        pit_ = 0;
    }

    CoalTable::CoalTable(const std::vector<double> &time)
    {
        time_ = time;
        pit_ = time_.size() + 1;
    }

    CoalTable::CoalTable(const int& size, const int& pit)
    {
        time_ = std::vector<double>(size);
        pit_ = pit;
    }

    int CoalTable::size() const
    {
        return this->time_.size();
    }

    double CoalTable::get_likelyhood() const
    {
        double log_llh = 0;
        double logC = 0;
        double power = 0;
        int branch_number = pit_;
        double time = 0;
        for (int i = 0; i < (*this).size(); ++i)
        {
            power += 2 / (branch_number * (branch_number-1)) * (time_[i] - time);
            logC += 2 * std::log(i + 1); 
            time = time_[i];
        }
        logC -= std::log((*this).size());
        logC -= (*this).size()*std::log(2);
        log_llh = -logC - power;
        return log_llh;
    }

    CoalTable CoalTable::operator+(const CoalTable &B)
    {
        CoalTable result(B.size() + (*this).size(), B.pit_ + (*this).pit_);
        int i = 0;
        int j = 0;
        double t1 = 0;
        double t2 = 0;
        while (i + j < result.size())
        {
            t1 = B.time_[i];
            t2 = (*this).time_[j];
            int k = i + j;
            if (t1 < t2)
            {
                result.time_[k] = (*this).time_[i];
                i++;
            }
            else
            {
                result.time_[k] = B.time_[j];
                j++;
            }
        }
        return result;
    }

    CoalTable CoalTable::operator-(const CoalTable &B)
    {
        CoalTable result((*this).size() - B.size(), (*this).pit_ - B.pit_);
        int i = 0;
        int j = 0;
        int k = 0;
        while (i < (*this).size())
        {
            if ((j < B.size()) && ((*this).time_[i] == B.time_[j]))
            {
                j++;
            }
            else
            {
                result.time_[k] = (*this).time_[i];
                k++;
            }
            i++;
        }
        return result;
    }











    // Scenario_Computer implementation

    Scenario_Computer::Scenario_Computer(
        const tsk_treeseq_t& ts, std::span<const double> parameters, 
        const std::vector<int>& sample_population,
        int admixed_data_index, int outgroup_data_index
    ) : ts_(ts)
    {
        sample_population_ = sample_population;
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts);

        nodes_ = std::vector<Node>(number_of_nodes + 1);

        set_parameters(parameters);

        outgroup_data_index_ = outgroup_data_index;
        admixed_data_index_ = admixed_data_index;
    }

    int Scenario_Computer::compute_lower_nodes_split_1(const tsk_tree_t& tree)
    {
        tsk_id_t root_id = tree.virtual_root;
        std::vector<tsk_id_t> stack{}; 
        stack.push_back(root_id);
        lower_nodes_split_1_.clear();

        tsk_id_t node_id;
        tsk_id_t child_id;
        while (!stack.empty())
        {
            node_id = stack.back(); // removes last element
            nodes_[node_id].population = base_pop_index_; // крашу от корня до разделения
            stack.pop_back();
            for (child_id = tree.left_child[node_id]; child_id != TSK_NULL; child_id = tree.right_sib[child_id])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child_id];
                if (child_time > split_1_time_)
                {
                    stack.push_back(child_id);
                }
                else
                {
                    lower_nodes_split_1_.push_back(child_id);
                }
            }
        }
        return 0;
    }

    int Scenario_Computer::compute_fixed_nodes(const tsk_tree_t& tree)
    {
        for (auto &node : nodes_)
        {
            node.is_fixed = 0;
        }
        // get samples from ts object
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        lower_split_2_out_linages_num_ = 0;
        // loop over all samples
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // for each sample from outgroup we mark parents until first (ancient) split. 
            if (sample_population_[samples[i]] == outgroup_data_index_)// Outgroup 
            {
                node_id = samples[i];   
                while ( 
                    (node_id != TSK_NULL) &&
                    (nodes_[node_id].is_fixed == 0)
                    // (tree.tree_sequence->tables->nodes.time[node_id] < split_1_time_) &&
                )
                {
                    // mark population for all outgroup nodes
                    if (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_)
                        nodes_[node_id].population = outgroup_pop_index_;
                    nodes_[node_id].is_fixed = true;
                    node_id = tree.parent[node_id];
                }
                // compute the number of linages at the moment of split 2 in outgrou
                if (
                    (node_id == TSK_NULL) ||
                    (tree.tree_sequence->tables->nodes.time[node_id] >= split_2_time_)
                    
                ) 
                    lower_split_2_out_linages_num_ += 1;
            }
        }
        return 0;
    }

    int Scenario_Computer::compute_skipped_nodes(const tsk_tree_t& tree)
    {
        int c = 0;
        for (auto &node : nodes_)
        {
            node.is_skipped = 0;
        }

        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        // loop over all samples
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // for each sample from outgroup we mark parents until first (ancient) split. 
            if (sample_population_[samples[i]] == outgroup_data_index_)// Outgroup 
            {
                node_id = samples[i];   
                while ( 
                    (node_id != TSK_NULL) &&
                    (nodes_[node_id].is_skipped == false) &&
                    (
                        (node_id == samples[i]) || 
                        (nodes_[tree.left_child[node_id]].is_fixed == nodes_[tree.right_child[node_id]].is_fixed)
                    )
                )
                {
                    // mark population for all outgroup nodes
                    nodes_[node_id].is_skipped = true;
                    node_id = tree.parent[node_id];
                }
            }
        }

        return c;
    }

    /** 
    * @brief Сгенерировать патерн 1213121 длины 2^n, оканчивающуюся нулем.
    */
    std::vector<uint64_t> Scenario_Computer::generate_abacaba(uint64_t n)
    {
        std::vector<uint64_t> seq(1ULL << n, 0);

        uint64_t k = 0;
        int i = 0;
        while (k < n)
        {
            k++;
            seq[i] = (1ULL << (k - 1));
            for (int j = 0; j < i; j++)
            {
                seq[j + i + 1] = seq[j];
            }
            i = 2*i + 1;
        }
        return seq;
    }

    int Scenario_Computer::compute_log_scenarios_num()
    {
        log_scenarios_num_ = 0;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            if (nodes_[id].is_fixed == false)
            {
                log_scenarios_num_++;
            }
        }
        return log_scenarios_num_;
    }

    int Scenario_Computer::set_subtree_population(
        const tsk_tree_t &tree, tsk_id_t subtree_root_id, int population_id, double end_time)
    {
        std::vector<tsk_id_t> stack{}; 
        stack.push_back(subtree_root_id);

        int c = 0;
        tsk_id_t node_id;
        tsk_id_t child_id;
        double time;
        while (!stack.empty())
        {
            node_id = stack.back(); // removes last element
            time = tree.tree_sequence->tables->nodes.time[node_id];
            if ((nodes_[node_id].population != outgroup_pop_index_) && (time >= end_time))
                nodes_[node_id].population = population_id;
            c++;
            stack.pop_back();
            
            for (child_id = tree.left_child[node_id]; child_id != TSK_NULL; child_id = tree.right_sib[child_id])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child_id];
                if ((child_time >= end_time) && (nodes_[child_id].population != outgroup_pop_index_))
                {
                    stack.push_back(child_id);
                }
                else if (
                    (population_id == ghost_pop_index_) && 
                    (child_time < end_time) && 
                    (nodes_[child_id].population != outgroup_pop_index_))
                {
                    c += set_subtree_population(tree, child_id, base_pop_index_, 0);
                }
            }
        } 
        return c;
    }

    int Scenario_Computer::set_scenario(const tsk_tree_t& tree, uint64_t scenario_id)
    {
        if (scenario_id >= (1ULL << log_scenarios_num_))
            return -1;
        
        sc_lineages_num_[0] = 0;
        sc_lineages_num_[1] = 0;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            if (nodes_[id].is_fixed == false)
            {
                if ((scenario_id%2 == 0) && (nodes_[id].population != base_pop_index_))
                {
                    set_subtree_population(tree, id, base_pop_index_, 0);
                }
                else if ((scenario_id%2 == 1) && (nodes_[id].population != ghost_pop_index_))
                {
                    set_subtree_population(tree, id, ghost_pop_index_, migration_time_); // это смущает, как я крашу остаток?
                }
                sc_lineages_num_[scenario_id%2] += 1; 
                scenario_id /= 2;
            }
            else
            {
                if (nodes_[id].population != 0)
                {
                    set_subtree_population(tree, id, base_pop_index_, 0); // внутри функции мы не идем красить аутгруппу -- все ок!
                }
                sc_lineages_num_[0] += 1;
            }
        }
        return sc_lineages_num_[0];
    }

    int Scenario_Computer::check_scenario(const tsk_tree_t& tree)
    {
        for (auto node_id : time_ordered_tree_nodes_ids_)
        {
            if ((nodes_[node_id].population == -1) && (tree.tree_sequence->tables->nodes.time[node_id] != 0))  
                return node_id;
        }
        return 0;
    }

    int Scenario_Computer::check_imposible_split_2(const tsk_tree_t &tree)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<bool> flag(number_of_nodes + 1, 0);
        // get samples from ts object
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        // loop over all samples
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // for each sample from outgroup we mark parents until second (early) split. 
            if (sample_population_[samples[i]] == admixed_data_index_)// Admixed 
            {
                node_id = samples[i];   
                while ( 
                    (node_id     != TSK_NULL) &&
                    (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_) &&
                    (flag[node_id] == 0)
                )
                {
                    if (nodes_[node_id].population == outgroup_pop_index_)
                    {
                        return node_id;
                    }
                    flag[node_id] = 1;
                    node_id = tree.parent[node_id];
                }
            }
        }
        return 0;   
    }


    double Scenario_Computer::find_lowest_coal(const tsk_tree_t& tree, std::span<const tsk_id_t> A, std::span<const tsk_id_t> B)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<bool> flag_A(number_of_nodes + 1, 0);
        std::vector<bool> flag_B(number_of_nodes + 1, 0);
        // get samples from ts object

        for (tsk_id_t node_id : A)
        {
            // for each sample from outgroup we mark parents until second (early) split. 
            while ( 
                (node_id != TSK_NULL) &&
                (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_) &&
                (flag_A[node_id] == 0)
            )
            {
                flag_A[node_id] = 1;
                node_id = tree.parent[node_id];
            }
        }

        double time = split_2_time_;
        for (tsk_id_t node_id : B)
        {
            while ( 
                (node_id != TSK_NULL) &&
                (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_) &&
                (flag_A[node_id] == 0) &&
                (flag_B[node_id] == 0)
            )
            {
                flag_B[node_id] = 1;
                node_id = tree.parent[node_id];
                if ((flag_A[node_id] == 1) && (tree.tree_sequence->tables->nodes.time[node_id] < time))
                {
                    time = tree.tree_sequence->tables->nodes.time[node_id];
                }
            }
        }
        return time;
    }


    int Scenario_Computer::sort_by_time(const tsk_tree_t& tree)
    {
        tsk_id_t virtual_root_id = tree.virtual_root;
        std::vector<tsk_id_t> part_ordered{}; 
        time_ordered_tree_nodes_ids_.clear();
        part_ordered.push_back(virtual_root_id);

        tsk_id_t node_id;
        tsk_id_t child;
        std::vector<tsk_id_t>::iterator it;
        while (!part_ordered.empty())
        {
            node_id = part_ordered.back(); // removes last element
            part_ordered.pop_back();
            if (node_id != virtual_root_id)
            {
                time_ordered_tree_nodes_ids_.push_back(node_id);
                nodes_[node_id].time = tree.tree_sequence->tables->nodes.time[node_id];
            }

            for (child = tree.left_child[node_id]; child != TSK_NULL; child = tree.right_sib[child])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child];
                it = part_ordered.begin();
                while (
                    (it != part_ordered.end()) 
                    && (tree.tree_sequence->tables->nodes.time[*it] < child_time))
                {
                    it++;
                }
                part_ordered.insert(it, child);
            }
        }
        return 0;     
    }

    /** 
    * @brief Вычисляет правдоподобие дерева на заданном участке.
    * 
    * @param tree объект дерева в формате `tsk_tree_t`.
    * @param time_start верхняя граница временного интервала.
    * @param time_end нижняя граница временного интервала.
    * @param populations вектор из популяций, для которых нужно вычислять правдоподобие.
    * @param init_lineages_num количество ветвей в дереве в каждой популяции в момент времени `time_start`
    * 
    * @return вычисленное правдоподобие. 
    */
    double Scenario_Computer::smart_llh(
        double time_start, double time_end, 
        std::span<const int> populations, std::span<int> lineages_num)
    {
        double llh = 0;

        int pop_num = populations.size();
        
        auto population_index = [populations](int population) {
            for (size_t i = 0; i < populations.size(); ++i)
            {
                if (populations[i] == population)
                    return static_cast<int>(i);
            }
            return -1;
        };

        std::vector<Node> current_node(pop_num, {-1, -1, 0, 0});
        Node next_node = {-1, -1, 0, 0};
        
        int i = 0;
        next_node = nodes_[time_ordered_tree_nodes_ids_[i]];
        while (((next_node.time > time_start) || (population_index(next_node.population) == -1)) && next_node.time != 0)
        {
            i++;
            next_node = nodes_[time_ordered_tree_nodes_ids_[i]];
        }

        int pop_i = population_index(next_node.population);
        while (next_node.time > time_end)
        {
            pop_i = population_index(next_node.population);
            if (pop_i != -1)
            {
                if (current_node[pop_i].time == -1)
                {
                    // llh without coal
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                        * (time_start - next_node.time)
                    );
                }
                else
                {
                    // llh with coal 
                    if (lineages_num[pop_i] > 1)
                    {
                        llh += std::log(1 / N_[pop_i] / 2.);
                        llh += (
                            -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                            * (current_node[pop_i].time - next_node.time)
                        );
                    } 
                }
                current_node[pop_i] = next_node;
                lineages_num[pop_i] += 1;
            }
            i++;
            next_node = nodes_[time_ordered_tree_nodes_ids_[i]];
        }
        // llh end
        for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
        {
            if (lineages_num[pop_i] > 1)
            {
                if (current_node[pop_i].time == -1)
                {
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2. 
                        * (time_start - time_end)
                    );  
                }
                else
                {
                    llh += std::log(1 / N_[pop_i] / 2.); // exact pair, earlier: (lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2.)
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                        * (current_node[pop_i].time - time_end)
                    );   
                }
            }
        }
        return llh;
    }

    /** 
    * @brief Вычисляет правдоподобие дерева на заданном участке.
    * 
    * @param tree объект дерева в формате `tsk_tree_t`.
    * @param time_start верхняя граница временного интервала.
    * @param time_end нижняя граница временного интервала.
    * @param populations вектор из популяций, для которых нужно вычислять правдоподобие.
    * @param init_lineages_num количество ветвей в дереве в каждой популяции в момент времени `time_start`
    * 
    * @return вычисленное правдоподобие. 
    */
    double Scenario_Computer::smart_llh_to_nodes(
        double time_start, double time_end, 
        std::span<const int> populations, 
        std::span<int> lineages_num, 
        DoubleMatrixView nodes_llh,
        IntMatrixView nodes_lineage_num,
        std::mdspan<const double, std::dextents<size_t, 2>> inv_2N_grid,
        std::mdspan<const double, std::dextents<size_t, 2>> log_inv_2N_grid,
        std::span<double> llh_step,
        size_t N_grid_size
    )
    {
        auto population_index = [populations](int population) {
            for (size_t i = 0; i < populations.size(); ++i)
            {
                if (populations[i] == population)
                    return static_cast<int>(i);
            }
            return -1;
        };

        tsk_id_t next_node_id = -1;
        tsk_id_t cur_node_id = -1;
        tsk_id_t cur_pop_i = -1;

        double current_time = time_start;
    
        
        size_t i = 0;
        next_node_id = time_ordered_tree_nodes_ids_[i];

        while (((nodes_[next_node_id].time > time_start) || (population_index(nodes_[next_node_id].population) == -1)) && nodes_[next_node_id].time != 0)
        {
            i++;
            next_node_id = time_ordered_tree_nodes_ids_[i];
        }

        while (nodes_[next_node_id].time > time_end)
        {
            const int next_pop_i = population_index(nodes_[next_node_id].population);
            if (next_pop_i != -1)
            {
                std::fill(llh_step.begin(), llh_step.end(), 0);
                for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
                {
                    if (cur_node_id != -1)
                    {
                        llh_step[n_i] = log_inv_2N_grid[populations[cur_pop_i], n_i]
                            + nodes_llh[cur_node_id, n_i];
                    }
                    for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                    {
                        llh_step[n_i] += (
                            -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. * inv_2N_grid[populations[pop_i], n_i]
                            * (current_time - nodes_[next_node_id].time)
                        );
                    }

                    nodes_llh[next_node_id, n_i] += llh_step[n_i];
                }

                for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                {
                    nodes_lineage_num[next_node_id, pop_i] = lineages_num[pop_i];
                }

                cur_node_id = next_node_id;
                cur_pop_i = next_pop_i;
                current_time = nodes_[next_node_id].time;
                lineages_num[cur_pop_i] += 1;
            }
            i++;
            next_node_id = time_ordered_tree_nodes_ids_[i];
        }
        // llh end
        // std::cout << "ln: " << lineages_num[pop_i] << std::endl;

        while ((i < time_ordered_tree_nodes_ids_.size()) && (population_index(nodes_[next_node_id].population) == -1))
        {
            // std::cout << next_pop << std::endl;
            i++;
            next_node_id = time_ordered_tree_nodes_ids_[i];
        }

        std::fill(llh_step.begin(), llh_step.end(), 0);
        for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
        {
            if (cur_node_id != -1) // it means we have coalescence in some population
            {
                llh_step[n_i] = log_inv_2N_grid[populations[cur_pop_i], n_i]
                    + nodes_llh[cur_node_id, n_i];
            }
            for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
            {
                llh_step[n_i] += (
                    -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. * inv_2N_grid[populations[pop_i], n_i]
                    * (current_time - time_end) // current_time = start_time if no coalescence obtained.
                );  
            }
            nodes_llh[next_node_id, n_i] += llh_step[n_i];
        }

        for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
        {
            // std::cout << pop_i << " " << lineages_num[pop_i] << ", ";
            nodes_lineage_num[next_node_id, pop_i] = lineages_num[pop_i];
        }
        // std::cout << std::endl;
        return nodes_llh[next_node_id, 0];
    }

    int Scenario_Computer::debug_sort_by_time(const tsk_tree_t &tree)
    {
        double t2, t1;
        t2 = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_ids_[0]];
        for (size_t i = 1; i < time_ordered_tree_nodes_ids_.size(); i++)
        {
            t1 = t2;
            t2 = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_ids_[i]];
            if (t2 > t1)
                return -1;
        }
        return 0;
    }


DenseTensor3D Scenario_Computer::compute_grid_fast(
        const tsk_tree_t& tree, 
        std::vector<double> migration_time,
        std::vector<double> migration_prob,
        std::vector<std::vector<double>> N
    ) 
    {
        if (N.empty())
        {
            N = {{N_[0]}, {N_[1]}, {N_[2]}};
        }
        if (N.size() > N_.size())
        {
            throw std::invalid_argument("N grid must contain at most one vector per population");
        }
        while (N.size() < N_.size())
        {
            N.push_back({N_[N.size()]});
        }

        size_t N_grid_size = 1;
        for (const auto &values : N)
        {
            N_grid_size = std::max(N_grid_size, values.size());
        }
        for (size_t pop_i = 0; pop_i < N.size(); ++pop_i)
        {
            if (N[pop_i].empty())
            {
                N[pop_i] = {N_[pop_i]};
            }
            if (N[pop_i].size() != 1 && N[pop_i].size() != N_grid_size)
            {
                throw std::invalid_argument("each N grid vector must have size 1 or the common N grid size");
            }
        }

        DenseMatrix inv_2N_matrix(N.size(), N_grid_size, 0);
        DenseMatrix log_inv_2N_matrix(N.size(), N_grid_size, 0);
        auto inv_2N = inv_2N_matrix.view();
        auto log_inv_2N = log_inv_2N_matrix.view();
        for (size_t pop_i = 0; pop_i < N.size(); ++pop_i)
        {
            for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
            {
                const double N_value = grid_N_value(N, static_cast<int>(pop_i), n_i);
                inv_2N[pop_i, n_i] = 0.5 / N_value;
                log_inv_2N[pop_i, n_i] = std::log(inv_2N[pop_i, n_i]);
            }
        }

        std::vector<double> log_migration_prob(migration_prob.size(), 0);
        std::vector<double> log_no_migration_prob(migration_prob.size(), 0);
        std::vector<int> migration_prob_kind(migration_prob.size(), 0);
        for (size_t j = 0; j < migration_prob.size(); ++j)
        {
            const double p = migration_prob[j];
            if (p == 0)
            {
                migration_prob_kind[j] = 1;
                log_migration_prob[j] = -std::numeric_limits<double>::infinity();
                log_no_migration_prob[j] = 0;
            }
            else if (p == 1)
            {
                migration_prob_kind[j] = 2;
                log_migration_prob[j] = 0;
                log_no_migration_prob[j] = -std::numeric_limits<double>::infinity();
            }
            else
            {
                log_migration_prob[j] = std::log(p);
                log_no_migration_prob[j] = std::log(1.0 - p);
            }
        }
    
        DenseTensor3D result_tensor(migration_time.size(), N_grid_size, migration_prob.size(), 0);
        DoubleTensor3DView result = result_tensor.view();

        tsk_size_t roots_num = tsk_tree_get_num_roots(&tree);

        if (roots_num != 1)
        {
            std::cout << "Roots number is not equal to 1. Skipping the tree." << std::endl;
            return result_tensor;
        }

        double migration_time_lb = migration_time.back();
        migration_time_ = migration_time_lb;
        std::vector<double> llh_root(N_grid_size, 0);

        for (auto &node : nodes_)
        {
            node.population = -1;
        }
        
        sort_by_time(tree);
        compute_lower_nodes_split_1(tree);
        compute_fixed_nodes(tree);

        // print_lower_nodes_split_1();

        int impossible_coal_node = check_imposible_split_2(tree);
        if (impossible_coal_node != 0)
        {
            std::cout << "impossible scenario (coalescence between two isolated populations)." << std::endl;
            std::cout << "Coalescence node (" << impossible_coal_node << ")." << std::endl;
            return result_tensor;
        }

        int log_scenario_num = compute_log_scenarios_num();

        std::vector<int> lineages_num;
        std::vector<int> populations;

        double root_time = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_ids_[0]];
        populations = {base_pop_index_}; 
        int linages_num_root[] = {1};

        if (root_time > split_1_time_)
        {
            std::vector<double> saved_N = N_;
            for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
            {
                for (size_t pop_i = 0; pop_i < N_.size(); ++pop_i)
                {
                    N_[pop_i] = grid_N_value(N, static_cast<int>(pop_i), n_i);
                }
                int root_lineages[] = {linages_num_root[0]};
                llh_root[n_i] += smart_llh(
                    root_time, split_1_time_, populations, root_lineages
                );
            }
            N_ = saved_N;
        }

        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<double> no_structure_llh_data(number_of_nodes * N_grid_size, 0);
        std::vector<int> no_structure_lineage_num_data(number_of_nodes, 0);
        DoubleMatrixView no_structure_llh(no_structure_llh_data.data(), number_of_nodes, N_grid_size);
        IntMatrixView no_structure_lineage_num(no_structure_lineage_num_data.data(), number_of_nodes, 1);
        std::vector<double> llh_step_buffer(N_grid_size, 0);

        uint64_t scenario = 0;
        set_scenario(tree, scenario);
        populations = {base_pop_index_}; 
        smart_llh_to_nodes(
            split_1_time_, split_2_time_, populations, sc_lineages_num_, 
            no_structure_llh, no_structure_lineage_num, inv_2N, log_inv_2N,
            llh_step_buffer, N_grid_size
        );
        lineages_num = { 
                sc_lineages_num_[0] - lower_split_2_out_linages_num_
            };
        smart_llh_to_nodes(
            split_2_time_, 0, populations, lineages_num, 
            no_structure_llh, no_structure_lineage_num, inv_2N, log_inv_2N,
            llh_step_buffer, N_grid_size
        );

        for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
        {
            size_t min_llh_i = 0;
            for (size_t i = 0; i < time_ordered_tree_nodes_ids_.size(); ++i)
            {
                const tsk_id_t node_id = time_ordered_tree_nodes_ids_[i];
                const tsk_id_t min_node_id = time_ordered_tree_nodes_ids_[min_llh_i];
                if (no_structure_llh[node_id, n_i] != 0)
                {
                    if (no_structure_llh[node_id, n_i] < no_structure_llh[min_node_id, n_i])
                        min_llh_i = i;
                }
            }

            double A = no_structure_llh[time_ordered_tree_nodes_ids_[min_llh_i], n_i];

            for (size_t i = min_llh_i + 1; i-- > 0;)
            {
                const tsk_id_t node_id = time_ordered_tree_nodes_ids_[i];
                no_structure_llh[node_id, n_i] = - (
                    no_structure_llh[node_id, n_i] - A
                );
            }
        }


        if (log_scenario_num >= 20)
        {
            std::cout << "Attention! Memory need: " << (1 << (log_scenario_num - 10 - 10)) * sizeof(double) * 4 * number_of_nodes / 1024. << " gB." << std::endl;
        }

#ifndef NDEBUG
        tsk_id_t err_id;
#endif
        std::vector<double> structure_llh_data(number_of_nodes * N_grid_size, 0);
        std::vector<int> structure_lineage_num_data(number_of_nodes * 2, 0);
        DoubleMatrixView structure_llh(structure_llh_data.data(), number_of_nodes, N_grid_size);
        IntMatrixView structure_lineage_num(structure_lineage_num_data.data(), number_of_nodes, 2);
        

        std::vector<uint64_t> abacaba = generate_abacaba(log_scenario_num);
        scenario = 0;
        int lineages_num_total;
        int lineages_between[2];
        size_t node_i;
        tsk_id_t node_id;
        tsk_id_t prev_node_id = -1;
        double llh_step_r = 0;
        double llh_step = 0;
        double llh_for_sc = 0;
        double mig_prob = 0;
        double t = 0;
        for (size_t sc_i = 0; sc_i < abacaba.size(); sc_i++)
        {
            std::fill(structure_llh_data.begin(), structure_llh_data.end(), 0);
            std::fill(structure_lineage_num_data.begin(), structure_lineage_num_data.end(), 0);
            // std::cout << std::bitset<6>(scenario) << std::endl;
            set_scenario(tree, scenario);
#ifndef NDEBUG
            err_id = check_scenario(tree);
            if (err_id != 0)
            {
                std::cout << "Error scenario: " << err_id << std::endl;
                std::cout << nodes_[err_id].is_fixed << std::endl;
                std::cout << nodes_[err_id].population << std::endl;
                std::cout << ts_.tables->nodes.population[err_id] << std::endl;
                std::cout << tree.tree_sequence->tables->nodes.time[err_id] << std::endl; 
                std::cout << nodes_[tree.parent[err_id]].population << std::endl;
                std::cout << ts_.tables->nodes.time[tree.parent[err_id]] << std::endl;
                break;
            }
#endif

            populations = {base_pop_index_, ghost_pop_index_}; 
            smart_llh_to_nodes(
                    split_1_time_, split_2_time_, populations, sc_lineages_num_, 
                    structure_llh, structure_lineage_num, inv_2N, log_inv_2N,
                    llh_step_buffer, N_grid_size
            );

            populations = {base_pop_index_, ghost_pop_index_}; 
            lineages_num = { 
                sc_lineages_num_[0] - lower_split_2_out_linages_num_,
                sc_lineages_num_[1]
            };
            // ||
            smart_llh_to_nodes(
                    split_2_time_, migration_time_lb, populations, lineages_num, 
                    structure_llh, structure_lineage_num, inv_2N, log_inv_2N,
                    llh_step_buffer, N_grid_size
            );

            // std::cout << lineages_num[0] << " / " << lineages_num[1] << std::endl;
            // std::cout << sc_i << " of " << pow(2, log_scenario_num) << " : " << log_prob << std::endl;
            node_i = 0;
            node_id = time_ordered_tree_nodes_ids_[node_i];
            prev_node_id = - 1;
            for (size_t i = 0; i < migration_time.size(); ++i)
            {
                t = migration_time[i];
                while ((nodes_[node_id].time > t) || (nodes_[node_id].population == outgroup_pop_index_))
                {
                    node_i++;
                    if (nodes_[node_id].population != outgroup_pop_index_)
                    {
                        prev_node_id = node_id;
                    }
                    node_id = time_ordered_tree_nodes_ids_[node_i];
                }

                if (prev_node_id == -1)
                {
                    if (sc_i == abacaba.size() - 1)
                    {
                        for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
                        {
                            for (size_t j = 0; j < migration_prob.size(); ++j)
                            {
                                result[i, n_i, j] = llh_root[n_i] + no_structure_llh[node_id, n_i];
                            }
                        }
                    }
                    continue;
                }

                lineages_num_total = 0;
                lineages_between[0] = structure_lineage_num[node_id, 0]; 
                lineages_between[1] = structure_lineage_num[node_id, 1]; 
                if (split_2_time_ <= nodes_[prev_node_id].time) // prev_time <= split_2 <= migration <= node_time
                    lineages_between[0] += lower_split_2_out_linages_num_;
                for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                {
                    lineages_num_total += structure_lineage_num[node_id, pop_i];
                }

                for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
                {
                    llh_step = log_inv_2N[nodes_[prev_node_id].population, n_i];
                    for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                    {
                        llh_step += (
                            -lineages_between[pop_i] * (lineages_between[pop_i] - 1) / 2. * inv_2N[populations[pop_i], n_i] 
                            * (nodes_[prev_node_id].time - split_2_time_) // between prev_time and split.
                        );
                        llh_step += (
                            -structure_lineage_num[node_id, pop_i] * (structure_lineage_num[node_id, pop_i] - 1) / 2. * inv_2N[populations[pop_i], n_i] 
                            * (split_2_time_ - t) // between split and migration removing splited lineages.
                        );
                    }

                    llh_for_sc = structure_llh[prev_node_id, n_i] + llh_step;
                    const int base_lineages = structure_lineage_num[node_id, 0];
                    const int ghost_lineages = structure_lineage_num[node_id, 1];
                    for (size_t j = 0; j < migration_prob.size(); ++j)
                    {
                        if (migration_prob_kind[j] == 1)
                        {
                            mig_prob = ghost_lineages == 0 ? 0 : -std::numeric_limits<double>::infinity();
                        }
                        else if (migration_prob_kind[j] == 2)
                        {
                            mig_prob = base_lineages == 0 ? 0 : -std::numeric_limits<double>::infinity();
                        }
                        else
                        {
                            mig_prob = base_lineages * log_no_migration_prob[j]
                                + ghost_lineages * log_migration_prob[j];
                        }
                        result[i, n_i, j] += std::exp(llh_for_sc + mig_prob);
                    }
                    if (sc_i == abacaba.size() - 1)
                    {
                        llh_step_r = (
                            -lineages_num_total * (lineages_num_total - 1) / 2. * inv_2N[base_pop_index_, n_i]
                            * (t - nodes_[node_id].time) // current_time = start_time if no coalescence obtained.
                        );
                        for (size_t j = 0; j < migration_prob.size(); ++j)
                        {
                            result[i, n_i, j] = llh_root[n_i] + std::log(result[i, n_i, j])
                                + no_structure_llh[node_id, n_i]
                                + llh_step_r;
                        }
                    }
                }
            }
            scenario = scenario ^ abacaba[sc_i];
        }
        return result_tensor;
    }


    double Scenario_Computer::compute_llh(const tsk_tree_t& tree)
    {
        double llh = 0;

        for (auto &node : nodes_)
        {
            node.population = -1;
        }

        sort_by_time(tree);
        compute_lower_nodes_split_1(tree);


        compute_fixed_nodes(tree);

        // std::cout << lower_nodes_split_1_.size() << std::endl;
        // for (auto node : lower_nodes_split_1_)
        // {
        //     std:: cout << node << ": " << nodes_[node].is_fixed << std::endl;
        // }

        int log_scenario_num = compute_log_scenarios_num();

        std::vector<int> lineages_num;
        std::vector<int> populations;

        double root_time = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_ids_[0]];
        populations = {base_pop_index_}; 
        int linages_num_root[] = {1};
        int lineage_num_lower_migration = 0;

        if (root_time > split_1_time_)
        {
            llh += smart_llh(
                root_time, split_1_time_, populations, linages_num_root
            );
        }
        
        // std::cout << llh << std::endl;

        double llh_temp = 0;
        double log_prob = 0;
        double prob_sum = 0;

        tsk_id_t err_id;
        // tsk_id_t impossible_coal_node;
        std::vector<uint64_t> abacaba = generate_abacaba(log_scenario_num);

        uint64_t scenario = 0;
        for (size_t sc_i = 0; sc_i < abacaba.size(); sc_i++)
        {
            // std::cout << std::bitset<6>(scenario) << std::endl;
            log_prob = 0;
            set_scenario(tree, scenario);
            // impossible_coal_node = check_imposible_split_2(tree);
            // if (impossible_coal_node != 0)
            // {
            //     std::cout << "impossible scenario (coalescence between two isolated populations)." << std::endl;
            //     std::cout << "Coalescence node (" << impossible_coal_node << ")." << std::endl;
            //     return std::log(0);
            // }
            err_id = check_scenario(tree);
            if (err_id != 0)
            {
                std::cout << "Error scenario: " << err_id << std::endl;
                std::cout << nodes_[err_id].is_fixed << std::endl;
                std::cout << nodes_[err_id].population << std::endl;
                std::cout << ts_.tables->nodes.population[err_id] << std::endl;
                std::cout << tree.tree_sequence->tables->nodes.time[err_id] << std::endl; 
                std::cout << nodes_[tree.parent[err_id]].population << std::endl;
                std::cout << ts_.tables->nodes.time[tree.parent[err_id]] << std::endl;
                break;
            }

            // fisrt split probabiltiy
            // log_prob += sc_lineages_num_[0] * std::log(split_1_prop_)
            //           + sc_lineages_num_[1] * std::log(1 - split_1_prop_);

            // /|
            // std::cout << sc_lineages_num_[0] << " / " << sc_lineages_num_[1] << std::endl;
            populations = {base_pop_index_, ghost_pop_index_}; 
            llh_temp = smart_llh(
                    split_1_time_, split_2_time_, populations, sc_lineages_num_
            );
            log_prob += llh_temp;

            // std::cout << llh_temp << std::endl;
            // second split probability
            // log_prob += lower_split_2_out_linages_num_ * std::log(1 - split_2_prop_)
            //           + (sc_lineages_num_[0] - lower_split_2_out_linages_num_) * std::log(split_2_prop_);
            // std::cout << sc_lineages_num_[0] - lower_split_2_out_linages_num_ << " / " << sc_lineages_num_[1] << std::endl;
            // std::cout << "lower_split_2_out_linages_num_: " << lower_split_2_out_linages_num_ << std::endl;
            populations = {base_pop_index_, ghost_pop_index_}; 
            lineages_num = { 
                sc_lineages_num_[0] - lower_split_2_out_linages_num_,
                sc_lineages_num_[1]
            };
            // ||
            llh_temp = smart_llh(
                 split_2_time_, migration_time_, populations, lineages_num
            );
            log_prob += llh_temp;
            // std::cout << llh_temp << std::endl;
            // migration probability
            // std::cout << migration_prob_ <<  std::endl;
            if (migration_prob_ != 1)
            {
                log_prob += lineages_num[0] * std::log(1 - migration_prob_) \
                        + lineages_num[1] * std::log(migration_prob_);
            }

            prob_sum += std::exp(log_prob);

            if (sc_i == 0)
                lineage_num_lower_migration = lineages_num[0] + lineages_num[1];
            // std::cout << lineages_num[0] << " / " << lineages_num[1] << std::endl;
            // std::cout << sc_i << " of " << pow(2, log_scenario_num) << " : " << log_prob << std::endl;
            scenario = scenario ^ abacaba[sc_i];
        }
        // std::cout << "between " << std::log(prob_sum) << std::endl;
        llh += std::log(prob_sum);
        // populations = {outgroup_pop_index_};
        // lineages_num = {lower_split_2_out_linages_num_};
        // llh += smart_llh(tree, split_2_time_, migration_time_, populations, lineages_num);

        populations = {base_pop_index_}; 
        lineages_num = {
            lineage_num_lower_migration
            // lineages_num[0]
        };
        // std::cout << lineages_num[0] << " " << lineages_num[1] << std::endl;
        llh_temp = smart_llh(migration_time_, 0, populations, lineages_num);
        // std::cout << "tail " << llh_temp << std::endl;
        llh += llh_temp;
        return llh;
    }

    void Scenario_Computer::set_parameters(std::span<const double> parameters)
    {
        migration_time_ = parameters[0];
        split_2_time_ = parameters[1];
        split_1_time_ = parameters[2];
        
        migration_prob_ = parameters[3];
        split_2_prop_ = parameters[4];
        split_1_prop_ = parameters[5];

        N_ = {parameters[6], parameters[7], parameters[8]};
    }

    std::vector<tsk_id_t> Scenario_Computer::get_lower_nodes_split_1()
    {
        return lower_nodes_split_1_;
    }

    void Scenario_Computer::print_lower_nodes_split_1()
    {
        for (tsk_id_t node : lower_nodes_split_1_)
        {
            std::cout << node << ": " << nodes_[node].is_fixed << std::endl;
        }
        std::cout << std::endl;
    }
}
