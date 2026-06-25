#include "treellh.h"

#include <limits>
#include <stdexcept>

namespace treellh
{
    namespace
    {
        // Вспомогательная функция для получения значения Ne (эффективного размера популяции)
        // из сетки N_grid для заданной популяции и индекса сетки n_i.
        // Если размер вектора равен 1, то значение Ne фиксировано.
        // Иначе возвращается значение Ne на шаге n_i сетки.
        double grid_N_value(
            const std::vector<std::vector<double>>& N_grid,
            int population,
            size_t n_i
        )
        {
            const auto &values = N_grid.at(static_cast<size_t>(population));
            return values.size() == 1 ? values[0] : values[n_i];
        }

        // Вспомогательная функция для расчета сочетаний C_n^k
        uint64_t choose(int n, int k)
        {
            if (k < 0 || k > n) return 0;
            if (k == 0 || k == n) return 1;
            if (k > n / 2) k = n - k;
            uint64_t res = 1;
            for (int i = 1; i <= k; ++i)
            {
                res = res * (n - k + i) / i;
            }
            return res;
        }
    }

    // =========================================================================
    // Реализация класса CoalTable (Вспомогательный класс для расчета коалесценции)
    // Представляет собой упорядоченный по времени список событий коалесценции.
    // =========================================================================

    // Конструктор по умолчанию
    CoalTable::CoalTable()
    {
        time_ = std::vector<double>();
        pit_ = 0;
    }

    // Конструктор от вектора времен коалесценции
    CoalTable::CoalTable(const std::vector<double> &time)
    {
        time_ = time;
        pit_ = time_.size() + 1;
    }

    // Конструктор от размера и начального числа линий
    CoalTable::CoalTable(const int& size, const int& pit)
    {
        time_ = std::vector<double>(size);
        pit_ = pit;
    }

    // Возвращает количество событий коалесценции в таблице
    int CoalTable::size() const
    {
        return this->time_.size();
    }

    // Вычисляет логарифм правдоподобия Kingman coalescent на основе времен в таблице.
    // Для интервалов между событиями добавляет вероятность отсутствия коалесценции,
    // а для каждого события коалесценции — плотность вероятности самого события коалесценции.
    double CoalTable::get_likelyhood() const
    {
        double log_llh = 0;
        double logC = 0;
        double power = 0;
        int branch_number = pit_;
        double time = 0;
        
        // Проходим по всем интервалам времени коалесценции
        for (int i = 0; i < (*this).size(); ++i)
        {
            // Формула: сумма по интервалам dt * (k*(k-1)/2)
            power += 2 / (branch_number * (branch_number-1)) * (time_[i] - time);
            logC += 2 * std::log(i + 1); 
            time = time_[i];
        }
        logC -= std::log((*this).size());
        logC -= (*this).size()*std::log(2);
        log_llh = -logC - power;
        return log_llh;
    }

    // Оператор слияния двух таблиц времен (сортировка слиянием)
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

    // Оператор разности двух таблиц времен (удаляет из текущей таблицы времена, которые есть в B)
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


    // =========================================================================
    // Реализация класса Scenario_Computer (Основной вычислитель правдоподобия)
    // =========================================================================

    // Конструктор вычислителя. Принимает ARG, вектор начальных параметров,
    // вектор популяций для каждого образца генома и индексы исследуемых популяций.
    Scenario_Computer::Scenario_Computer(
        const tsk_treeseq_t& ts, std::span<const double> parameters, 
        const std::vector<int>& sample_population,
        int admixed_data_index, int outgroup_data_index
    ) : ts_(ts)
    {
        sample_population_ = sample_population;
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts);

        // Инициализируем массив внутренних узлов (+1 для виртуального корня)
        nodes_ = std::vector<Node>(number_of_nodes + 1);

        set_parameters(parameters);

        outgroup_data_index_ = outgroup_data_index; // Индекс аутгруппы (SAM) в исходных данных
        admixed_data_index_ = admixed_data_index;   // Индекс африканской популяции (AFR) в исходных данных
    }

    // Находит узлы дерева, которые находятся непосредственно ниже разделения split_1_time_.
    // Проводится обход в ширину (BFS) от корня дерева. Узлы выше split_1_time_ сразу окрашиваются
    // в предковую популяцию base_pop_index_ (популяция ARCH).
    // Первые встреченные узлы, чье время жизни меньше или равно split_1_time_, сохраняются в lower_nodes_split_1_.
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
            node_id = stack.back();
            nodes_[node_id].population = base_pop_index_; // Красим в предковую популяцию выше split 1
            stack.pop_back();
            
            // Проходим по всем потомкам текущего узла
            for (child_id = tree.left_child[node_id]; child_id != TSK_NULL; child_id = tree.right_sib[child_id])
            {
                double child_time = tree.tree_sequence->tables->nodes.time[child_id];
                if (child_time > split_1_time_)
                {
                    // Если потомок старше split_1_time_, продолжаем обход вглубь
                    stack.push_back(child_id);
                }
                else
                {
                    // Если потомок моложе split_1_time_, это граница сценариев. Сохраняем его.
                    lower_nodes_split_1_.push_back(child_id);
                }
            }
        }
        return 0;
    }

    // Помечает фиксированные узлы, однозначно принадлежащие аутгруппе (SAM).
    // Для каждого образца из аутгруппы поднимаемся вверх к предкам. Пока время жизни узлов
    // меньше split_2_time_, они красятся в популяцию аутгруппы (outgroup_pop_index_).
    // Все такие узлы помечаются как is_fixed = true, так как они не могут участвовать
    // в миграции от призрачной популяции.
    // Также подсчитывается количество независимых линий аутгруппы (SAM) на границе split_2_time_.
    int Scenario_Computer::compute_fixed_nodes(const tsk_tree_t& tree)
    {
        for (auto &node : nodes_)
        {
            node.is_fixed = 0;
        }
        
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        lower_split_2_out_linages_num_ = 0;
        
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            // Проверяем, принадлежит ли образец аутгруппе SAM
            if (sample_population_[samples[i]] == outgroup_data_index_) 
            {
                node_id = samples[i];   
                while ( 
                    (node_id != TSK_NULL) &&
                    (nodes_[node_id].is_fixed == 0)
                )
                {
                    // Красим узлы ниже разделения split_2 в популяцию аутгруппы
                    if (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_)
                        nodes_[node_id].population = outgroup_pop_index_;
                    nodes_[node_id].is_fixed = true;
                    node_id = tree.parent[node_id]; // Поднимаемся к предку
                }
                
                // Считаем число линий аутгруппы в момент времени split_2
                if (
                    (node_id == TSK_NULL) ||
                    (tree.tree_sequence->tables->nodes.time[node_id] >= split_2_time_)
                ) 
                    lower_split_2_out_linages_num_ += 1;
            }
        }
        return 0;
    }

    // Находит узлы, которые можно пропустить при расчете сценариев, так как они полностью
    // лежат внутри изолированной ветви аутгруппы (SAM) и не пересекаются с африканской ветвью.
    int Scenario_Computer::compute_skipped_nodes(const tsk_tree_t& tree)
    {
        int c = 0;
        for (auto &node : nodes_)
        {
            node.is_skipped = 0;
        }

        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            if (sample_population_[samples[i]] == outgroup_data_index_)
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
                    nodes_[node_id].is_skipped = true;
                    node_id = tree.parent[node_id];
                }
            }
        }

        return c;
    }

    // Генерирует последовательность кода Грея «abacaba» (например, для n=3: 1, 2, 1, 4, 1, 2, 1, 0).
    // Позволяет перебрать все 2^n сценариев распределения популяций для свободных линий.
    // На каждом шаге меняется состояние ровно одного бита (одного узла), что позволяет пересчитывать
    // правдоподобие инкрементально, а не обходить всё дерево с нуля.
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

    // Вычисляет число свободных узлов (линий, чья эволюционная судьба не фиксирована аутгруппой)
    // на уровне ниже split_1. Это число определяет размерность пространства сценариев: всего 2^n сценариев.
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

    // Рекурсивно окрашивает поддерево с корнем в subtree_root_id в заданную популяцию population_id
    // вглубь по времени до достижения end_time.
    // Обрабатывает переключение популяций на границах времени миграции (Ghost переходит в Base/AFR ниже migration_time).
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
            node_id = stack.back();
            time = tree.tree_sequence->tables->nodes.time[node_id];
            
            // Если узел не принадлежит аутгруппе и его время больше end_time, присваиваем ему популяцию
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
                // Если мы красили в Ghost, но дошли до границы времени интрогрессии (migration_time_):
                // Линии, идущие дальше в прошлое (моложе migration_time_), должны быть перекрашены в base (AFR)
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

    // Назначает аутгруппе (SAM) соответствующий индекс популяции
    int Scenario_Computer::set_outgroup(const tsk_tree_t& tree, tsk_id_t* samples, tsk_size_t samples_num)
    {
        return 0;
    }

    // Применяет сценарий с номером scenario_id к дереву.
    // Свободные линии из lower_nodes_split_1_ красятся в base (0) или в ghost (1)
    // в зависимости от битов числа scenario_id.
    // Возвращает количество линий, отнесенных к базовой популяции.
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
                // Если нулевой бит равен 0, красим поддерево в base (AFR)
                if ((scenario_id%2 == 0) && (nodes_[id].population != base_pop_index_))
                {
                    set_subtree_population(tree, id, base_pop_index_, 0);
                }
                // Если нулевой бит равен 1, красим поддерево в ghost (призрачную) популяцию
                else if ((scenario_id%2 == 1) && (nodes_[id].population != ghost_pop_index_))
                {
                    set_subtree_population(tree, id, ghost_pop_index_, migration_time_);
                }
                sc_lineages_num_[scenario_id%2] += 1; 
                scenario_id /= 2; // Сдвиг к следующему биту
            }
            else
            {
                // Фиксированные узлы по умолчанию красятся в базовую популяцию
                if (nodes_[id].population != 0)
                {
                    set_subtree_population(tree, id, base_pop_index_, 0);
                }
                sc_lineages_num_[0] += 1;
            }
        }
        return sc_lineages_num_[0];
    }

    // Проверяет корректность разметки сценария.
    // Каждый узел с ненулевым временем должен иметь назначенную популяцию (population != -1).
    // Возвращает id некорректного узла или 0, если все в порядке.
    int Scenario_Computer::check_scenario(const tsk_tree_t& tree)
    {
        for (auto node_id : time_ordered_tree_nodes_ids_)
        {
            if ((nodes_[node_id].population == -1) && (tree.tree_sequence->tables->nodes.time[node_id] != 0))  
                return node_id;
        }
        return 0;
    }

    // Проверяет невозможные с биологической точки зрения коалесценции.
    // Если две популяции физически изолированы друг от друга в определенный интервал времени
    // (например, AFR и SAM моложе split_2_time_), они не могут коалесцировать (иметь общего предка).
    // Если такая коалесценция обнаружена, данный сценарий имеет нулевую вероятность.
    int Scenario_Computer::check_imposible_split_2(const tsk_tree_t &tree)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<bool> flag(number_of_nodes + 1, 0);
        const tsk_id_t *samples = tsk_treeseq_get_samples(tree.tree_sequence);
        const tsk_size_t samples_num = tsk_treeseq_get_num_samples(tree.tree_sequence);
        
        tsk_id_t node_id;
        for (tsk_size_t i = 0; i < samples_num; i++)
        {
            if (sample_population_[samples[i]] == admixed_data_index_)
            {
                node_id = samples[i];   
                while ( 
                    (node_id     != TSK_NULL) &&
                    (tree.tree_sequence->tables->nodes.time[node_id] < split_2_time_) &&
                    (flag[node_id] == 0)
                )
                {
                    // Если линия AFR пересеклась с аутгруппой SAM ниже split_2
                    if (nodes_[node_id].population == outgroup_pop_index_)
                    {
                        return node_id; // Возвращаем проблемный узел
                    }
                    flag[node_id] = 1;
                    node_id = tree.parent[node_id];
                }
            }
        }
        return 0;   
    }

    // Вспомогательная функция поиска самого раннего времени коалесценции
    // между линиями множества A и множества B ниже split_2_time_.
    double Scenario_Computer::find_lowest_coal(const tsk_tree_t& tree, std::span<const tsk_id_t> A, std::span<const tsk_id_t> B)
    {
        tsk_size_t number_of_nodes = tsk_treeseq_get_num_nodes(&ts_);
        std::vector<bool> flag_A(number_of_nodes + 1, 0);
        std::vector<bool> flag_B(number_of_nodes + 1, 0);

        for (tsk_id_t node_id : A)
        {
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

    // Сортирует узлы дерева по времени их возникновения в порядке убывания (от корня к листьям).
    // Результат сохраняется в time_ordered_tree_nodes_ids_. Это необходимо, чтобы при расчете
    // правдоподобия мы могли двигаться по дереву строго назад во времени (снизу вверх)
    // и последовательно обрабатывать события коалесценции.
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
            node_id = part_ordered.back();
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

    // Вычисляет логарифм правдоподобия коалесценции на временном интервале [time_end, time_start]
    // для заданного набора популяций и начального числа линий.
    // Реализует модель Кингмана:
    // - Коалесценция k линий дает вклад log(1 / 2Ne)
    // - Не-коалесценция k линий на интервале dt дает вклад -k*(k-1)/2 * dt / 2Ne
    double Scenario_Computer::smart_llh(
        double time_start, double time_end, 
        std::span<const int> populations, std::span<int> lineages_num)
    {
        double llh = 0;
        int pop_num = populations.size();
        
        // Лямбда-функция для поиска индекса популяции в локальном массиве
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
        
        // Пропускаем узлы старше верхней границы интервала time_start
        while (((next_node.time > time_start) || (population_index(next_node.population) == -1)) && next_node.time != 0)
        {
            i++;
            next_node = nodes_[time_ordered_tree_nodes_ids_[i]];
        }

        int pop_i = population_index(next_node.population);
        
        // Основной цикл: идем вниз по времени (в прошлое) до нижней границы time_end
        while (next_node.time > time_end)
        {
            pop_i = population_index(next_node.population);
            if (pop_i != -1)
            {
                if (current_node[pop_i].time == -1)
                {
                    // Если это первое событие в этой популяции на интервале,
                    // считаем вероятность не-коалесценции от начала интервала time_start до этого события
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                        * (time_start - next_node.time)
                    );
                }
                else
                {
                    // Если в популяции уже были события коалесценции:
                    if (lineages_num[pop_i] > 1)
                    {
                        // Добавляем вероятность самого события слияния: log(1 / 2Ne)
                        llh += std::log(1 / N_[pop_i] / 2.);
                        
                        // Добавляем вероятность не-коалесценции на интервале между соседними событиями
                        llh += (
                            -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                            * (current_node[pop_i].time - next_node.time)
                        );
                    } 
                }
                current_node[pop_i] = next_node;
                lineages_num[pop_i] += 1; // Увеличиваем число линий (так как идем назад во времени, число линий увеличивается при ветвлении)
            }
            i++;
            next_node = nodes_[time_ordered_tree_nodes_ids_[i]];
        }

        // Хвостовой расчет: добавляем вероятность не-коалесценции от последнего события до конца интервала time_end
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
                    llh += std::log(1 / N_[pop_i] / 2.);
                    llh += (
                        -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. / N_[pop_i] / 2.
                        * (current_node[pop_i].time - time_end)
                    );   
                }
            }
        }
        return llh;
    }

    // Векторизованная версия smart_llh. Вычисляет правдоподобие коалесценции одновременно
    // для всей сетки эффективных размеров популяций Ne (размера N_grid_size).
    // Позволяет избежать многократного прохода по структуре дерева для разных Ne.
    double Scenario_Computer::smart_llh_to_nodes(
        double time_start, double time_end, 
        std::span<const int> populations, 
        std::span<int> lineages_num, 
        DoubleMatrixView nodes_llh,
        IntMatrixView nodes_lineage_num,
        DoubleMatrixView inv_2N_grid,
        DoubleMatrixView log_inv_2N_grid,
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
                        llh_step[n_i] = log_inv_2N_grid(populations[cur_pop_i], n_i)
                            + nodes_llh(cur_node_id, n_i);
                    }
                    for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                    {
                        llh_step[n_i] += (
                            -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. * inv_2N_grid(populations[pop_i], n_i)
                            * (current_time - nodes_[next_node_id].time)
                        );
                    }

                    nodes_llh(next_node_id, n_i) += llh_step[n_i];
                }

                for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                {
                    nodes_lineage_num(next_node_id, pop_i) = lineages_num[pop_i];
                }

                cur_node_id = next_node_id;
                cur_pop_i = next_pop_i;
                current_time = nodes_[next_node_id].time;
                lineages_num[cur_pop_i] += 1;
            }
            i++;
            next_node_id = time_ordered_tree_nodes_ids_[i];
        }

        while ((i < time_ordered_tree_nodes_ids_.size()) && (population_index(nodes_[next_node_id].population) == -1))
        {
            i++;
            next_node_id = time_ordered_tree_nodes_ids_[i];
        }

        std::fill(llh_step.begin(), llh_step.end(), 0);
        for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
        {
            if (cur_node_id != -1)
            {
                llh_step[n_i] = log_inv_2N_grid(populations[cur_pop_i], n_i)
                    + nodes_llh(cur_node_id, n_i);
            }
            for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
            {
                llh_step[n_i] += (
                    -lineages_num[pop_i] * (lineages_num[pop_i] - 1) / 2. * inv_2N_grid(populations[pop_i], n_i)
                    * (current_time - time_end)
                );  
            }
            nodes_llh(next_node_id, n_i) += llh_step[n_i];
        }

        for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
        {
            nodes_lineage_num(next_node_id, pop_i) = lineages_num[pop_i];
        }
        return nodes_llh(next_node_id, 0);
    }

    // Проверка корректности сортировки узлов по времени
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

    // =========================================================================
    // Основная оптимизированная функция расчета правдоподобия на сетке параметров
    // =========================================================================
    DenseTensor3D Scenario_Computer::compute_grid_fast(
        const tsk_tree_t& tree, 
        std::vector<double> migration_time,
        std::vector<double> migration_prob,
        std::vector<std::vector<double>> N
    ) 
    {
        // Инициализируем сетку Ne если она пуста
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

        // Предварительно вычисляем обратные величины 1/(2Ne) для быстрого умножения вместо деления
        DenseMatrix inv_2N_matrix(N.size(), N_grid_size, 0);
        DenseMatrix log_inv_2N_matrix(N.size(), N_grid_size, 0);
        auto inv_2N = inv_2N_matrix.view();
        auto log_inv_2N = log_inv_2N_matrix.view();
        for (size_t pop_i = 0; pop_i < N.size(); ++pop_i)
        {
            for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
            {
                const double N_value = grid_N_value(N, static_cast<int>(pop_i), n_i);
                inv_2N(pop_i, n_i) = 0.5 / N_value;
                log_inv_2N(pop_i, n_i) = std::log(inv_2N(pop_i, n_i));
            }
        }

        // Предварительно логарифмируем вероятности миграции для ускорения расчетов
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

        // Если у дерева не один корень (мульти-дерево), пропускаем его
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
        
        // Подготовка дерева: сортировка узлов, поиск нижних узлов, определение фиксированных
        sort_by_time(tree);
        compute_fixed_nodes(tree);
        compute_lower_nodes_split_1(tree);

        // Отбраковываем дерево, если в нем есть недопустимая коалесценция ниже split 2
        int impossible_coal_node = check_imposible_split_2(tree);
        if (impossible_coal_node != 0)
        {
            std::cout << "impossible scenario (coalescence between two isolated populations)." << std::endl;
            std::cout << "Coalescence node (" << impossible_coal_node << ")." << std::endl;
            return result_tensor;
        }

        // Идентификация неактивных и активных узлов
        inactive_nodes_.clear();
        std::vector<tsk_id_t> active_nodes;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            if (nodes_[id].is_fixed == false)
            {
                if (tree.tree_sequence->tables->nodes.time[id] <= migration_time_lb)
                {
                    inactive_nodes_.push_back(id);
                }
                else
                {
                    active_nodes.push_back(id);
                }
            }
        }

        int log_scenario_num = active_nodes.size();

        std::vector<int> lineages_num;
        std::vector<int> populations;

        // Вычисляем правдоподобие предковой популяции ARCH (выше split 1)
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
        
        // Матрицы для хранения промежуточных результатов безструктурного (базового) правдоподобия
        std::vector<double> no_structure_llh_data(number_of_nodes * N_grid_size, 0);
        std::vector<int> no_structure_lineage_num_data(number_of_nodes, 0);
        DoubleMatrixView no_structure_llh(no_structure_llh_data.data(), number_of_nodes, N_grid_size);
        IntMatrixView no_structure_lineage_num(no_structure_lineage_num_data.data(), number_of_nodes, 1);
        std::vector<double> llh_step_buffer(N_grid_size, 0);

        // Рассчитываем правдоподобие для базового "бесструктурного" сценария (scenario = 0)
        uint64_t scenario = 0;
        set_scenario(tree, scenario);
        populations = {base_pop_index_}; 
        
        // От split 1 до split 2
        int lineages_num_copy[] = {sc_lineages_num_[0], sc_lineages_num_[1]};
        smart_llh_to_nodes(
            split_1_time_, split_2_time_, populations, lineages_num_copy, 
            no_structure_llh, no_structure_lineage_num, inv_2N, log_inv_2N,
            llh_step_buffer, N_grid_size
        );
        // От split 2 до 0
        lineages_num = { 
                lineages_num_copy[0] - lower_split_2_out_linages_num_
            };
        smart_llh_to_nodes(
            split_2_time_, 0, populations, lineages_num, 
            no_structure_llh, no_structure_lineage_num, inv_2N, log_inv_2N,
            llh_step_buffer, N_grid_size
        );

        // Нормировка правдоподобий во избежание численного переполнения
        for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
        {
            size_t min_llh_i = 0;
            for (size_t i = 0; i < time_ordered_tree_nodes_ids_.size(); ++i)
            {
                const tsk_id_t node_id = time_ordered_tree_nodes_ids_[i];
                const tsk_id_t min_node_id = time_ordered_tree_nodes_ids_[min_llh_i];
                if (no_structure_llh(node_id, n_i) != 0)
                {
                    if (no_structure_llh(node_id, n_i) < no_structure_llh(min_node_id, n_i))
                        min_llh_i = i;
                }
            }

            double A = no_structure_llh(time_ordered_tree_nodes_ids_[min_llh_i], n_i);

            for (size_t i = min_llh_i + 1; i-- > 0;)
            {
                const tsk_id_t node_id = time_ordered_tree_nodes_ids_[i];
                no_structure_llh(node_id, n_i) = - (
                    no_structure_llh(node_id, n_i) - A
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
        // Матрицы для структурированного сценария (с разделением на Ghost и Base)
        std::vector<double> structure_llh_data(number_of_nodes * N_grid_size, 0);
        std::vector<int> structure_lineage_num_data(number_of_nodes * 2, 0);
        DoubleMatrixView structure_llh(structure_llh_data.data(), number_of_nodes, N_grid_size);
        IntMatrixView structure_lineage_num(structure_lineage_num_data.data(), number_of_nodes, 2);
        
        // Генерируем Gray-код abacaba для эффективного перебора всех 2^n сценариев
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

        auto flip_active_node = [&](tsk_id_t n_id) {
            if (nodes_[n_id].population == base_pop_index_)
            {
                set_subtree_population(tree, n_id, ghost_pop_index_, migration_time_lb);
                sc_lineages_num_[0]--;
                sc_lineages_num_[1]++;
            }
            else
            {
                set_subtree_population(tree, n_id, base_pop_index_, 0);
                sc_lineages_num_[0]++;
                sc_lineages_num_[1]--;
            }
        };

        size_t n = inactive_nodes_.size();
        std::vector<bool> inactive_is_ghost(n, false);
        auto flip_inactive_node = [&](size_t idx) {
            if (!inactive_is_ghost[idx])
            {
                inactive_is_ghost[idx] = true;
                sc_lineages_num_[0]--;
                sc_lineages_num_[1]++;
            }
            else
            {
                inactive_is_ghost[idx] = false;
                sc_lineages_num_[0]++;
                sc_lineages_num_[1]--;
            }
        };

        // Инициализация состояний для змейки
        sc_lineages_num_[0] = lower_nodes_split_1_.size();
        sc_lineages_num_[1] = 0;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            set_subtree_population(tree, id, base_pop_index_, 0);
        }

        // Основной цикл по сценариям
        for (size_t sc_i = 0; sc_i < abacaba.size(); sc_i++)
        {
            bool moving_up = (sc_i % 2 == 0);
            for (size_t step = 0; step <= n; ++step)
            {
                size_t k = moving_up ? step : (n - step);
                double comb_prob = choose(n, k);

                std::fill(structure_llh_data.begin(), structure_llh_data.end(), 0);
                std::fill(structure_lineage_num_data.begin(), structure_lineage_num_data.end(), 0);
#ifndef NDEBUG
                err_id = check_scenario(tree);
                if (err_id != 0)
                {
                    std::cout << "Error scenario: " << err_id << std::endl;
                    break;
                }
#endif

                // Считаем коалесценцию на первом интервале (от split 1 до split 2)
                populations = {base_pop_index_, ghost_pop_index_}; 
                int lineages_num_copy[] = {sc_lineages_num_[0], sc_lineages_num_[1]};
                smart_llh_to_nodes(
                        split_1_time_, split_2_time_, populations, lineages_num_copy, 
                        structure_llh, structure_lineage_num, inv_2N, log_inv_2N,
                        llh_step_buffer, N_grid_size
                );

                // Считаем коалесценцию на втором интервале (от split 2 до нижней границы сетки времени миграции)
                populations = {base_pop_index_, ghost_pop_index_}; 
                lineages_num = { 
                    lineages_num_copy[0] - lower_split_2_out_linages_num_,
                    lineages_num_copy[1]
                };
                smart_llh_to_nodes(
                        split_2_time_, migration_time_lb, populations, lineages_num, 
                        structure_llh, structure_lineage_num, inv_2N, log_inv_2N,
                        llh_step_buffer, N_grid_size
                );

                node_i = 0;
                node_id = time_ordered_tree_nodes_ids_[node_i];
                prev_node_id = - 1;
                
                // Итерируемся по сетке времен миграции
                for (size_t i = 0; i < migration_time.size(); ++i)
                {
                    t = migration_time[i];
                    // Ищем узел дерева, непосредственно предшествующий времени миграции t
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
                        if (sc_i == abacaba.size() - 1 && step == n)
                        {
                            for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
                            {
                                for (size_t j = 0; j < migration_prob.size(); ++j)
                                {
                                    result(i, n_i, j) = llh_root[n_i] + no_structure_llh(node_id, n_i);
                                }
                            }
                        }
                        continue;
                    }

                    lineages_num_total = 0;
                    lineages_between[0] = structure_lineage_num(node_id, 0); 
                    lineages_between[1] = structure_lineage_num(node_id, 1); 
                    if (split_2_time_ <= nodes_[prev_node_id].time)
                        lineages_between[0] += lower_split_2_out_linages_num_;
                    for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                    {
                        lineages_num_total += structure_lineage_num(node_id, pop_i);
                    }

                    for (size_t n_i = 0; n_i < N_grid_size; ++n_i)
                    {
                        llh_step = log_inv_2N(nodes_[prev_node_id].population, n_i);
                        for (size_t pop_i = 0; pop_i < populations.size(); pop_i++)
                        {
                            llh_step += (
                                -lineages_between[pop_i] * (lineages_between[pop_i] - 1) / 2. * inv_2N(populations[pop_i], n_i) 
                                * (nodes_[prev_node_id].time - split_2_time_)
                            );
                            llh_step += (
                                -structure_lineage_num(node_id, pop_i) * (structure_lineage_num(node_id, pop_i) - 1) / 2. * inv_2N(populations[pop_i], n_i) 
                                * (split_2_time_ - t)
                            );
                        }

                        llh_for_sc = structure_llh(prev_node_id, n_i) + llh_step;
                        const int base_lineages = structure_lineage_num(node_id, 0);
                        const int ghost_lineages = structure_lineage_num(node_id, 1);
                        
                        // Умножаем на вероятности миграции для текущего сценария (доля Ghost линий)
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
                            result(i, n_i, j) += comb_prob * std::exp(llh_for_sc + mig_prob);
                        }
                        
                        // По окончании перебора всех сценариев переводим в логарифмический масштаб и добавляем хвост
                        if (sc_i == abacaba.size() - 1 && step == n)
                        {
                            llh_step_r = (
                                -lineages_num_total * (lineages_num_total - 1) / 2. * inv_2N(base_pop_index_, n_i)
                                * (t - nodes_[node_id].time)
                            );
                            for (size_t j = 0; j < migration_prob.size(); ++j)
                            {
                                result(i, n_i, j) = llh_root[n_i] + std::log(result(i, n_i, j))
                                    + no_structure_llh(node_id, n_i)
                                    + llh_step_r;
                            }
                        }
                    }
                }
                
                // Переключение неактивной линии
                if (step < n)
                {
                    size_t flip_idx = moving_up ? step : (n - 1 - step);
                    flip_inactive_node(flip_idx);
                }
            }

            // Переход к следующему активному сценарию по коду Грея
            uint64_t diff = abacaba[sc_i];
            if (diff != 0)
            {
                int bit_idx = 0;
                uint64_t temp = diff;
                while ((temp & 1) == 0) {
                    temp >>= 1;
                    bit_idx++;
                }
                flip_active_node(active_nodes[bit_idx]);
            }
        }
        return result_tensor;
    }

    // Вычисляет суммарное логарифмическое правдоподобие для одного дерева
    // при фиксированных параметрах (без использования сеток).
    // Полезно для быстрой проверки правдоподобия конкретной точки в пространстве параметров.
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

        // Идентификация неактивных и активных узлов
        inactive_nodes_.clear();
        std::vector<tsk_id_t> active_nodes;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            if (nodes_[id].is_fixed == false)
            {
                if (tree.tree_sequence->tables->nodes.time[id] <= migration_time_)
                {
                    inactive_nodes_.push_back(id);
                }
                else
                {
                    active_nodes.push_back(id);
                }
            }
        }

        int log_scenario_num = active_nodes.size();

        std::vector<int> lineages_num;
        std::vector<int> populations;

        double root_time = tree.tree_sequence->tables->nodes.time[time_ordered_tree_nodes_ids_[0]];
        populations = {base_pop_index_}; 
        int linages_num_root[] = {1};
        int lineage_num_lower_migration = 0;

        // Расчет над split 1 (ARCH предки)
        if (root_time > split_1_time_)
        {
            llh += smart_llh(
                root_time, split_1_time_, populations, linages_num_root
            );
        }
        
        double llh_temp = 0;
        double log_prob = 0;
        double prob_sum = 0;

        tsk_id_t err_id;
        std::vector<uint64_t> abacaba = generate_abacaba(log_scenario_num);

        auto flip_active_node = [&](tsk_id_t n_id) {
            if (nodes_[n_id].population == base_pop_index_)
            {
                set_subtree_population(tree, n_id, ghost_pop_index_, migration_time_);
                sc_lineages_num_[0]--;
                sc_lineages_num_[1]++;
            }
            else
            {
                set_subtree_population(tree, n_id, base_pop_index_, 0);
                sc_lineages_num_[0]++;
                sc_lineages_num_[1]--;
            }
        };

        size_t n = inactive_nodes_.size();
        std::vector<bool> inactive_is_ghost(n, false);
        auto flip_inactive_node = [&](size_t idx) {
            if (!inactive_is_ghost[idx])
            {
                inactive_is_ghost[idx] = true;
                sc_lineages_num_[0]--;
                sc_lineages_num_[1]++;
            }
            else
            {
                inactive_is_ghost[idx] = false;
                sc_lineages_num_[0]++;
                sc_lineages_num_[1]--;
            }
        };

        // Инициализация состояний для змейки
        sc_lineages_num_[0] = lower_nodes_split_1_.size();
        sc_lineages_num_[1] = 0;
        for (tsk_id_t id : lower_nodes_split_1_)
        {
            set_subtree_population(tree, id, base_pop_index_, 0);
        }

        // Суммируем по всем сценариям (с змейкой для неактивных)
        for (size_t sc_i = 0; sc_i < abacaba.size(); sc_i++)
        {
            bool moving_up = (sc_i % 2 == 0);
            for (size_t step = 0; step <= n; ++step)
            {
                size_t k = moving_up ? step : (n - step);
                double comb_prob = choose(n, k);

                log_prob = 0;
#ifndef NDEBUG
                err_id = check_scenario(tree);
                if (err_id != 0)
                {
                    std::cout << "Error scenario: " << err_id << std::endl;
                    break;
                }
#endif

                // Коалесценция на интервале split 1 -> split 2 (популяции AFR/Ghost)
                populations = {base_pop_index_, ghost_pop_index_}; 
                int lineages_num_copy[] = {sc_lineages_num_[0], sc_lineages_num_[1]};
                llh_temp = smart_llh(
                        split_1_time_, split_2_time_, populations, lineages_num_copy
                );
                log_prob += llh_temp;

                // Коалесценция на интервале split 2 -> migration_time
                populations = {base_pop_index_, ghost_pop_index_}; 
                lineages_num = { 
                    lineages_num_copy[0] - lower_split_2_out_linages_num_,
                    lineages_num_copy[1]
                };
                llh_temp = smart_llh(
                     split_2_time_, migration_time_, populations, lineages_num
                );
                log_prob += llh_temp;
                
                // Взвешиваем сценарий по вероятности миграции (долям линий, ушедших в Ghost)
                if (migration_prob_ != 1)
                {
                    log_prob += lineages_num[0] * std::log(1 - migration_prob_) \
                            + lineages_num[1] * std::log(migration_prob_);
                }

                prob_sum += std::exp(log_prob) * comb_prob;

                if (sc_i == 0 && step == 0)
                    lineage_num_lower_migration = lineages_num[0] + lineages_num[1];

                // Переключение неактивной линии
                if (step < n)
                {
                    size_t flip_idx = moving_up ? step : (n - 1 - step);
                    flip_inactive_node(flip_idx);
                }
            }
            
            // Переход к следующему активному сценарию по коду Грея
            uint64_t diff = abacaba[sc_i];
            if (diff != 0)
            {
                int bit_idx = 0;
                uint64_t temp = diff;
                while ((temp & 1) == 0) {
                    temp >>= 1;
                    bit_idx++;
                }
                flip_active_node(active_nodes[bit_idx]);
            }
        }
        
        llh += std::log(prob_sum);

        // Хвост ниже migration_time (все линии объединились в базовую популяцию AFR)
        populations = {base_pop_index_}; 
        lineages_num = { lineage_num_lower_migration };
        llh_temp = smart_llh(migration_time_, 0, populations, lineages_num);
        llh += llh_temp;
        
        return llh;
    }

    // Загрузка параметров модели из вектора из 9 элементов.
    void Scenario_Computer::set_parameters(std::span<const double> parameters)
    {
        migration_time_ = parameters[0];
        split_2_time_ = parameters[1];
        split_1_time_ = parameters[2];
        
        migration_prob_ = parameters[3];
        split_2_prop_ = parameters[4]; // Доля разделения 2, в текущей модели не используется
        split_1_prop_ = parameters[5]; // Доля разделения 1, в текущей модели не используется

        N_ = {parameters[6], parameters[7], parameters[8]};
    }

    // Возвращает список узлов, лежащих под split 1
    std::vector<tsk_id_t> Scenario_Computer::get_lower_nodes_split_1()
    {
        return lower_nodes_split_1_;
    }

    // Вывод в консоль информации по нижним узлам split 1
    void Scenario_Computer::print_lower_nodes_split_1()
    {
        for (tsk_id_t node : lower_nodes_split_1_)
        {
            std::cout << node << ": " << nodes_[node].is_fixed << std::endl;
        }
        std::cout << std::endl;
    }
}
