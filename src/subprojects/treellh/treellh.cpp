#include "treellh.h"


namespace treellh
{

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
            result.time_[i+j] = (t1 < t2) ? (*this).time_[i++] : B.time_[j++];
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
        const tsk_treeseq_t& ts, std::span<const double> parameters, int admixed_data_index, int outgroup_data_index) : ts_(ts)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts);
        is_fixed_ = std::vector<bool>(number_of_nodes, 0); 
        population_ = std::vector<int>(number_of_nodes, -1);

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

        tsk_id_t node;
        tsk_id_t child;
        while (!stack.empty())
        {
            node = stack.back(); // removes last element
            population_[node] = base_pop_index_; // крашу от корня до разделения
            stack.pop_back();
            for (child = tree.left_child[node]; child != TSK_NULL; child = tree.right_sib[child])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child];
                if (child_time > split_1_time_)
                {
                    stack.push_back(child);
                }
                else
                {
                    lower_nodes_split_1_.push_back(child);
                }
            }
        }
        return 0;
    }

    int Scenario_Computer::compute_fixed_nodes(const tsk_tree_t& tree)
    {
        std::fill(is_fixed_.begin(), is_fixed_.end(), 0);
        // get samples from ts object
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        lower_split_2_out_linages_num_ = 0;
        // loop over all samples
        tsk_id_t node;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // for each sample from outgroup we mark parents until first (ancient) split. 
            if (tree.tree_sequence->tables->nodes.population[samples[i]] == outgroup_data_index_)// Outgroup 
            {
                node = samples[i];   
                while ( 
                    (node != TSK_NULL) &&
                    (tree.tree_sequence->tables->nodes.time[node] < split_1_time_) &&
                    (is_fixed_[node] == 0)
                )
                {
                    // mark population for all outgroup nodes
                    if (tree.tree_sequence->tables->nodes.time[node] < split_2_time_)
                        population_[node] = outgroup_pop_index_;
                    is_fixed_[node] = 1;
                    node = tree.parent[node];
                }
                // compute the number of linages at the moment of split 2 in outgrou
                if (tree.tree_sequence->tables->nodes.time[node] >= split_2_time_)
                    lower_split_2_out_linages_num_ += 1;
            }
        }
        return 0;
    }

    /** 
    * @brief Сгенерировать патерн 1213121 длины 2^n, оканчивающуюся нулем.
    */
    std::vector<uint64_t> Scenario_Computer::generate_abacaba(uint64_t n)
    {
        std::vector<uint64_t> seq(std::pow(2, n), 0);

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
            if (is_fixed_[id] == false)
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
        tsk_id_t node;
        tsk_id_t child;
        double time;
        while (!stack.empty())
        {
            node = stack.back(); // removes last element
            time = tree.tree_sequence->tables->nodes.time[node];
            if ((population_[node] != outgroup_pop_index_) && (time > end_time))
                population_[node] = population_id;
            c++;
            stack.pop_back();
            
            for (child = tree.left_child[node]; child != TSK_NULL; child = tree.right_sib[child])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child];
                if ((child_time > end_time) && (population_[child] != outgroup_pop_index_))
                {
                    stack.push_back(child);
                }
                else if (
                    (population_id == ghost_pop_index_) && 
                    (child_time <= end_time) && 
                    (population_[child] != outgroup_pop_index_))
                {
                    c += set_subtree_population(tree, child, base_pop_index_, 0);
                }
            }
        } 
        return c;
    }

    int Scenario_Computer::set_scenario(const tsk_tree_t& tree, int scenario_id)
    {
        if (scenario_id >= pow(2, log_scenarios_num_))
            return -1;
        
        sc_lineages_num_[0] = 0;
        sc_lineages_num_[1] = 0;
        // std::cout << "scenario: " << scenario_id << std::endl;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            if (is_fixed_[id] == false)
            {
                // std::cout << scenario_id%2;
                if ((scenario_id%2 == 0) && (population_[id] != base_pop_index_))
                {
                    set_subtree_population(tree, id, base_pop_index_, 0);
                }
                else if ((scenario_id%2 == 1) && (population_[id] != ghost_pop_index_))
                {
                    set_subtree_population(tree, id, ghost_pop_index_, migration_time_); // это смущает, как я крашу остаток?
                }
                sc_lineages_num_[scenario_id%2] += 1; 
                scenario_id /= 2;
            }
            else
            {
                if (population_[id] != 0)
                {
                    set_subtree_population(tree, id, base_pop_index_, 0); // внутри функции мы не идем красить аутгруппу -- все ок!
                }
                sc_lineages_num_[0] += 1;
            }
        }
        // std::cout << std::endl;
        return sc_lineages_num_[0];
    }

    int Scenario_Computer::check_scenario(const tsk_tree_t& tree)
    {
        for (auto node : time_ordered_tree_nodes_)
        {
            if ((population_[node] == -1) && (tree.tree_sequence->tables->nodes.time[node] != 0))  
                return node;
        }
        return 0;
    }

    int Scenario_Computer::check_imposible_split_2(const tsk_tree_t &tree)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<bool> flag(number_of_nodes, 0);
        // get samples from ts object
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        // loop over all samples
        tsk_id_t node;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // for each sample from outgroup we mark parents until second (early) split. 
            if (tree.tree_sequence->tables->nodes.population[samples[i]] == admixed_data_index_)// Admixed 
            {
                node = samples[i];   
                while ( 
                    (node != TSK_NULL) &&
                    (tree.tree_sequence->tables->nodes.time[node] < split_2_time_) &&
                    (flag[node] == 0)
                )
                {
                    if (population_[node] == outgroup_pop_index_)
                    {
                        return node;
                    }
                    flag[node] = 1;
                    node = tree.parent[node];
                }
            }
        }
        return 0;   
    }


    int Scenario_Computer::sort_by_time(const tsk_tree_t& tree)
    {
        tsk_id_t virtual_root_id = tree.virtual_root;
        std::vector<tsk_id_t> part_ordered{}; 
        time_ordered_tree_nodes_.clear();
        part_ordered.push_back(virtual_root_id);

        tsk_id_t node;
        tsk_id_t child;
        std::vector<tsk_id_t>::iterator it;
        while (!part_ordered.empty())
        {
            node = part_ordered.back(); // removes last element
            part_ordered.pop_back();
            if (node != virtual_root_id)
                time_ordered_tree_nodes_.push_back(node);

            for (child = tree.left_child[node]; child != TSK_NULL; child = tree.right_sib[child])
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
    double Scenario_Computer::smart_llh(const tsk_tree_t &tree, 
        double time_start, double time_end, 
        std::span<const int> populations, std::span<int> lineages_num)
    {
        double llh = 0;

        int pop_num = populations.size();
        

        std::map<int, int> pop_to_i; 
        for (size_t i = 0; i < populations.size(); ++i)
        {
            pop_to_i[populations[i]] = i;
        }

        std::vector<tsk_id_t> current_node(pop_num, -1);
        tsk_id_t next_node = -1;

        std::vector<double> current_time(pop_num);
        
        int i = 0;
        next_node = time_ordered_tree_nodes_[i];
        double next_time = tree.tree_sequence->tables->nodes.time[next_node];
        int next_pop = population_[next_node];

        while (((next_time > time_start) || (pop_to_i.find(next_pop) == pop_to_i.end())) && next_time != 0)
        {
            i++;
            next_node = time_ordered_tree_nodes_[i];
            next_time = tree.tree_sequence->tables->nodes.time[next_node];
            next_pop = population_[next_node];
        }

        int pop_i = pop_to_i[next_pop];
        while (next_time > time_end)
        {
            if (pop_to_i.find(next_pop) != pop_to_i.end())
            {
                pop_i = pop_to_i[next_pop];
                if (current_node[pop_i] == -1)
                {
                    // llh without coal
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                        * (time_start - next_time)
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
                            * (current_time[pop_i] - next_time)
                        );
                    } 
                }
                current_node[pop_i] = next_node;
                current_time[pop_i] = next_time;
                lineages_num[pop_i] += 1;
            }
            i++;
            next_node = time_ordered_tree_nodes_[i];
            next_time = tree.tree_sequence->tables->nodes.time[next_node];
            next_pop = population_[next_node];
            // std::cout << next_node << " pop: " << next_pop << std::endl;
        }
        // llh end
        for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
        {
            // std::cout << "ln: " << lineages_num[pop_i] << std::endl;
            if (lineages_num[pop_i] > 1)
            {
                if (current_node[pop_i] == -1)
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
                        * (current_time[pop_i] - time_end)
                    );   
                }
            }
        }
        return llh;
    }

    int Scenario_Computer::debug_sort_by_time(const tsk_tree_t &tree)
    {
        double t2, t1;
        t2 = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_[0]];
        for (size_t i = 1; i < time_ordered_tree_nodes_.size(); i++)
        {
            t1 = t2;
            t2 = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_[i]];
            if (t2 > t1)
                return -1;
        }
        return 0;
    }

    double Scenario_Computer::compute_llh(const tsk_tree_t& tree)
    {
        double llh = 0;

        std::fill(population_.begin(), population_.end(), -1);

        sort_by_time(tree);
        compute_lower_nodes_split_1(tree);


        compute_fixed_nodes(tree);

        // std::cout << lower_nodes_split_1_.size() << std::endl;
        // for (auto node : lower_nodes_split_1_)
        // {
        //     std:: cout << node << ": " << is_fixed_[node] << std::endl;
        // }

        int log_scenario_num = compute_log_scenarios_num();

        std::vector<int> lineages_num;
        std::vector<int> populations;

        double root_time = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_[0]];
        populations = {base_pop_index_}; 
        int linages_num_root[] = {1};
        int lineage_num_lower_migration = 0;

        if (root_time > split_1_time_)
        {
            llh += smart_llh(
                tree, root_time, split_1_time_, populations, linages_num_root
            );
        }
        // std::cout << llh << std::endl;

        double llh_temp = 0;
        double log_prob = 0;
        double prob_sum = 0;

        tsk_id_t err;
        // tsk_id_t impossible_coal_node;
        std::vector<uint64_t> abacaba = generate_abacaba(log_scenario_num);

        uint64_t scenario = 0;
        for (size_t sc_i = 0; sc_i < abacaba.size(); sc_i++)
        {
            std::cout << std::bitset<6>(scenario) << std::endl;
            log_prob = 0;
            set_scenario(tree, scenario);
            // impossible_coal_node = check_imposible_split_2(tree);
            // if (impossible_coal_node != 0)
            // {
            //     std::cout << "impossible scenario (coalescence between two isolated populations)." << std::endl;
            //     std::cout << "Coalescence node (" << impossible_coal_node << ")." << std::endl;
            //     return std::log(0);
            // }
            err = check_scenario(tree);
            if (err != 0)
            {
                std::cout << "Error scenario: " << err << std::endl;
                std::cout << is_fixed_[err] << std::endl;
                std::cout << population_[err] << std::endl;
                std::cout << ts_.tables->nodes.population[err] << std::endl;
                std::cout << tree.tree_sequence->tables->nodes.time[err] << std::endl; 
                std::cout << population_[tree.parent[err]] << std::endl;
                std::cout << ts_.tables->nodes.time[tree.parent[err]] << std::endl;
                break;
            }

            // fisrt split probabiltiy
            // log_prob += sc_lineages_num_[0] * std::log(split_1_prop_) \
            //           + sc_lineages_num_[1] * std::log(1 - split_1_prop_);

            // /|
            // std::cout << sc_lineages_num_[0] << " / " << sc_lineages_num_[1] << std::endl;
            populations = {base_pop_index_, ghost_pop_index_}; 
            llh_temp = smart_llh(
                    tree, split_1_time_, split_2_time_, populations, sc_lineages_num_
            );
            log_prob += llh_temp;

            // second split probability
            // log_prob += lower_split_2_out_linages_num_ * std::log(1 - split_2_prop_) \
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
                    tree, split_2_time_, migration_time_, populations, lineages_num
            );
            log_prob += llh_temp;

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
        llh += std::log(prob_sum);
        populations = {outgroup_pop_index_};
        lineages_num = {lower_split_2_out_linages_num_};
        llh += smart_llh(tree, split_2_time_, migration_time_, populations, lineages_num);

        populations = {base_pop_index_, outgroup_pop_index_}; 
        lineages_num = {
            lineage_num_lower_migration,
            lineages_num[0]
        };
        // std::cout << lineages_num[0] << " " << lineages_num[1] << std::endl;
        llh += smart_llh(tree, migration_time_, 0, populations, lineages_num);
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
            std::cout << node << std::endl;
        }
        std::cout << std::endl;
    }
}