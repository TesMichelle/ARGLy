// =============================================================================
// treellh.h -- Ядро вычисления правдоподобия для модели демографии популяций
// =============================================================================
//
// Реализует вычисление log-правдоподобия L(theta | ARG) для предковых графов
// рекомбинации (ARG). Каждое локальное дерево оценивается независимо, а полное
// правдоподобие -- произведение по деревьям:
//
//   L(theta | ARG) = П_i  P(Ti | theta)
//   P(Ti | theta) = sum_A  P(Ti | theta, A) * P(A)
//
// где A -- «сценарий» (см. Scenario_Computer), theta -- параметры модели.
//
// Популяционная модель (идём назад во времени):
//
//   t=0            Современность
//   t=migration_time  Призрачная популяция Ghost «вливается» в AFR (base)
//   t=split_2_time    ANC разделяется на Ghost и AFR (второй split)
//   t=split_1_time    ARCH разделяется на ANC и SAM (первый split, outgroup)
//
// Индексы популяций:
//   base_pop_index_     = 0  (AFR, основная современная популяция)
//   ghost_pop_index_    = 1  (Ghost, призрачная)
//   outgroup_pop_index_ = 2  (SAM, аутгруппа)
//
// Вероятности коалесценции (модель Кингмана):
//   P(k линий не коалесцируют за dt) = exp( -k*(k-1)/2 * 1/(2Ne) * dt )
//   P(коалесценция k линий) = 1/(2Ne)   (для одной конкретной пары)
// =============================================================================

#include <vector>
#include <span>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <tskit.h>

#include <iostream>

namespace treellh
{
    // =========================================================================
    // MatrixView<T> -- невладеющая обёртка над двумерным массивом данных.
    // Позволяет обращаться к элементам через operator()(row, col) без копирования.
    // =========================================================================
    template <typename T>
    class MatrixView
    {
        private:
            T *data_;
            size_t rows_;
            size_t cols_;
        public:
            MatrixView(T *data, size_t rows, size_t cols)
                : data_(data), rows_(rows), cols_(cols) {}

            size_t rows() const { return rows_; }
            size_t cols() const { return cols_; }
            // Доступ к элементу (row, col) -- row-major порядок.
            T& operator()(size_t row, size_t col) const { return data_[row * cols_ + col]; }
    };

    // =========================================================================
    // Tensor3DView<T> -- невладеющая обёртка над трёхмерным массивом данных.
    // Используется для результата compute_grid_fast:
    //   dim0 = размер сетки по migration_time
    //   dim1 = размер сетки по Ne (эффективный размер)
    //   dim2 = размер сетки по migration_prob
    // =========================================================================
    template <typename T>
    class Tensor3DView
    {
        private:
            T *data_;
            size_t dim0_;
            size_t dim1_;
            size_t dim2_;
        public:
            Tensor3DView(T *data, size_t dim0, size_t dim1, size_t dim2)
                : data_(data), dim0_(dim0), dim1_(dim1), dim2_(dim2) {}

            size_t dim0() const { return dim0_; }
            size_t dim1() const { return dim1_; }
            size_t dim2() const { return dim2_; }
            // Доступ operator()(i,j,k) -- row-major, k меняется быстрее всего.
            T& operator()(size_t i, size_t j, size_t k) const
            {
                return data_[(i * dim1_ + j) * dim2_ + k];
            }
    };

    // Удобные псевдонимы типов для часто используемых специализаций.
    using DoubleMatrixView      = MatrixView<double>;
    using ConstDoubleMatrixView = MatrixView<const double>;
    using IntMatrixView         = MatrixView<int>;
    using DoubleTensor3DView    = Tensor3DView<double>;
    using ConstDoubleTensor3DView = Tensor3DView<const double>;

    // =========================================================================
    // DenseMatrix -- владеющая двумерная матрица double (row-major).
    // Предназначена для хранения промежуточных данных (inv_2N, log_inv_2N и т.п.)
    // =========================================================================
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

            // Возвращает невладеющий MatrixView для передачи во вспомогательные функции.
            DoubleMatrixView view() { return DoubleMatrixView(data_.data(), rows_, cols_); }
            ConstDoubleMatrixView view() const
            {
                return ConstDoubleMatrixView(data_.data(), rows_, cols_);
            }

            double& operator()(size_t row, size_t col) { return data_[row * cols_ + col]; }
            double operator()(size_t row, size_t col) const { return data_[row * cols_ + col]; }
    };

    // =========================================================================
    // DenseTensor3D -- владеющий трёхмерный тензор double (row-major).
    // Основной контейнер результатов compute_grid_fast:
    //   result[i][n_i][j] = суммированное log-правдоподобие для
    //     i-й точки сетки migration_time, n_i-й Ne, j-й migration_prob.
    // =========================================================================
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
            ConstDoubleTensor3DView view() const
            {
                return ConstDoubleTensor3DView(data_.data(), dim0_, dim1_, dim2_);
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

    // =========================================================================
    // CoalTable -- вспомогательный класс для хранения упорядоченных по времени
    // событий коалесценции (используется в устаревшей логике get_likelyhood).
    // pit_ -- начальное число ветвей (branches) до первой коалесценции.
    // time_ -- отсортированный вектор времён коалесцентных событий.
    // =========================================================================
    class CoalTable
    {
        private:
            std::vector<double> time_;
            int pit_;
        public:
            CoalTable();
            // Создаёт таблицу из вектора уже известных времён коалесценции.
            CoalTable(const std::vector<double>&);
            // Создаёт пустую таблицу заданного размера с указанным pit.
            CoalTable(const int&, const int&);

            int size() const;

            // Вычисляет log-правдоподобие по хранящимся временам коалесценции
            // (модель Кингмана с постоянным Ne).
            double get_likelyhood() const;

            // Объединяет две таблицы коалесценции (слияние двух популяций).
            CoalTable operator+(const CoalTable &B);
            // Вычитает таблицу B из текущей (разделение популяций).
            CoalTable operator-(const CoalTable &B);
    };

    // =========================================================================
    // Node -- состояние одного узла дерева в контексте вычисления.
    // Хранит, какой популяции принадлежит узел в данном сценарии,
    // его время, и служебные флаги.
    // =========================================================================
    struct Node
    {
        int population  = -1;   // Индекс популяции (-1 = не назначено)
        double time     = -1;   // Время узла (в поколениях)
        bool is_fixed   = false; // true = популяция узла зафиксирована (outgroup)
        bool is_skipped = false; // true = узел пропускается при суммировании сценариев
    };

    // =========================================================================
    // Scenario_Computer -- главный класс вычисления правдоподобия.
    //
    // Реализует суммирование по всем возможным «сценариям» -- назначениям
    // предковых линий (lineages) по популяциям в момент split_1_time.
    // Для каждого дерева ARG вычисляется:
    //
    //   P(T | theta) = sum_A P(T | theta, A) * P(A)
    //
    // где A -- конкретное бинарное распределение неопределённых линий
    // между base (AFR) и ghost популяциями.
    //
    // Параметры модели (вектор из 9 чисел):
    //   [0] migration_time  -- время интрогрессии Ghost -> AFR (в поколениях)
    //   [1] split_2_time    -- время разделения ANC на Ghost и AFR
    //   [2] split_1_time    -- время разделения ARCH на ANC и SAM
    //   [3] migration_prob  -- доля линий, пришедших из Ghost
    //   [4] split_2_prop    -- доля при втором split (не используется в текущей версии)
    //   [5] split_1_prop    -- доля при первом split (не используется в текущей версии)
    //   [6] N_base          -- эффективный размер base (AFR)
    //   [7] N_ghost         -- эффективный размер ghost
    //   [8] N_outgroup      -- эффективный размер outgroup (SAM)
    // =========================================================================
    class Scenario_Computer
    {
        private:
            const tsk_treeseq_t &ts_;          // Ссылка на весь tree sequence (ARG)
            std::vector<int> sample_population_; // Популяция каждого образца (по индексу образца)

            std::vector<Node> nodes_;            // Состояния всех узлов текущего дерева

            std::vector<bool> is_fixed_;         // (не используется напрямую вне nodes_)
            std::vector<int> population_;        // (не используется напрямую вне nodes_)

            // Узлы дерева, находящиеся непосредственно ниже split_1_time --
            // именно они являются «свободными» переменными сценария.
            std::vector<tsk_id_t> lower_nodes_split_1_;
            std::vector<tsk_id_t> inactive_nodes_;

            // Все узлы текущего дерева, отсортированные по убыванию времени
            // (от корня к листьям). Используется в smart_llh для прохода по событиям.
            std::vector<tsk_id_t> time_ordered_tree_nodes_ids_;

            // Число аутгрупповых (outgroup) линий ниже split_2_time.
            int lower_split_2_out_linages_num_;

            int node_below_split_1_index_;       // (вспомогательный индекс)

            // log2 числа сценариев = число «свободных» узлов в lower_nodes_split_1_.
            int log_scenarios_num_;

            // Количество линий в каждой популяции для текущего сценария:
            //   sc_lineages_num_[0] = линии base (AFR)
            //   sc_lineages_num_[1] = линии ghost
            int sc_lineages_num_[2];

            // Число линий в момент миграции (ниже split_2_time, до migration_time).
            int lineages_num_lower_migration_;

            // Параметры модели (см. set_parameters):
            double migration_time_;
            double migration_prob_;
            double split_1_time_;
            double split_2_time_;
            double split_1_prop_;
            double split_2_prop_;

            // Эффективные размеры популяций: N_[0]=base, N_[1]=ghost, N_[2]=outgroup.
            std::vector<double> N_;

            // Индексы популяций в данных (sample_population_):
            int outgroup_data_index_; // индекс SAM в входных данных
            int admixed_data_index_;  // индекс AFR (или admixed) во входных данных

            // Индексы популяций в модели (всегда фиксированные):
            int base_pop_index_     = 0; // AFR
            int ghost_pop_index_    = 1; // Ghost
            int outgroup_pop_index_ = 2; // SAM

        public:
            // Конструктор: инициализирует вычислятель для данного tree sequence.
            // ts                -- tree sequence (ARG)
            // parameters        -- вектор из 9 параметров модели
            // sample_population -- популяция каждого образца
            // admixed_data_index   -- индекс AFR (admixed) в sample_population
            // outgroup_data_index  -- индекс SAM (outgroup) в sample_population
            Scenario_Computer(const tsk_treeseq_t&, std::span<const double>, const std::vector<int>&, int, int);

            // Находит все узлы дерева, лежащие непосредственно ниже split_1_time.
            // Результат записывается в lower_nodes_split_1_.
            // Это узлы, принадлежность которых популяции неизвестна (сценарии).
            int compute_lower_nodes_split_1(const tsk_tree_t&);

            // Помечает узлы флагом is_fixed, если их популяция известна точно
            // (outgroup-линии прослеживаются до split, не имея неоднозначности).
            // Заодно считает lower_split_2_out_linages_num_.
            int compute_fixed_nodes(const tsk_tree_t&);

            // Помечает узлы флагом is_skipped -- те, что принадлежат аутгруппе
            // и не участвуют в суммировании сценариев.
            int compute_skipped_nodes(const tsk_tree_t&);

            // Возвращает число «свободных» узлов (log2 числа сценариев).
            int compute_log_scenarios_num();

            int compute_scenarious(); // (устаревшее/заглушка)

            // Генерирует последовательность «abacaba» длины 2^n --
            // последовательность XOR-масок, которая позволяет перебрать все 2^n
            // сценариев, меняя за раз только один бит (код Грея).
            std::vector<uint64_t> generate_abacaba(uint64_t);

            // Рекурсивно красит (присваивает популяцию) все узлы поддерева,
            // начиная с subtree_root_id, вплоть до времени end_time.
            int set_subtree_population(const tsk_tree_t&, tsk_id_t, int, double);

            // Присваивает SAM-популяцию аутгрупповым образцам.
            int set_outgroup(const tsk_tree_t&, tsk_id_t*, tsk_size_t);

            // Применяет конкретный сценарий с номером scenario_id:
            // распределяет «свободные» линии lower_nodes_split_1_ между
            // base (бит=0) и ghost (бит=1) популяциями и обновляет sc_lineages_num_.
            int set_scenario(const tsk_tree_t&, uint64_t);

            // Проверяет корректность сценария: у каждого не-нулевого узла должна
            // быть назначена популяция. Возвращает 0 при успехе или id ошибочного узла.
            int check_scenario(const tsk_tree_t&);

            // Проверяет, нет ли коалесценции между изолированными популяциями
            // (физически невозможный сценарий для split_2).
            int check_imposible_split_2(const tsk_tree_t&);

            // Находит минимальное время коалесценции между линиями двух групп A и B
            // ниже split_2_time. Используется для эвристики нижней границы split_2.
            double find_lowest_coal(const tsk_tree_t&, std::span<const tsk_id_t>, std::span<const tsk_id_t>);

            // Сортирует узлы текущего дерева по убыванию времени (от корня к листьям)
            // и сохраняет результат в time_ordered_tree_nodes_ids_.
            // Используется как preprocessing перед smart_llh.
            int sort_by_time(const tsk_tree_t&);

            // Отладочная проверка корректности сортировки.
            int debug_sort_by_time(const tsk_tree_t&);

            // Вычисляет log-правдоподобие коалесценции на интервале [time_end, time_start]
            // для заданного набора популяций и начального числа линий.
            // Реализует модель Кингмана: проходит по time_ordered_tree_nodes_ids_
            // и для каждого события коалесценции добавляет:
            //   log(1/(2Ne)) - k*(k-1)/2 * 1/(2Ne) * dt
            // Версия для одиночного набора Ne.
            double smart_llh(double, double, 
                std::span<const int>, std::span<int>);

            // Аналог smart_llh, но работает одновременно для сетки значений Ne (N_grid_size).
            // Заполняет матрицы nodes_llh и nodes_lineage_num для передачи результатов
            // в compute_grid_fast. Это позволяет не повторять обход дерева для каждого Ne.
            double smart_llh_to_nodes(double, double, 
                std::span<const int>, std::span<int>, DoubleMatrixView, IntMatrixView,
                DoubleMatrixView,
                DoubleMatrixView,
                std::span<double>, size_t);

            // Основная функция: вычисляет log-правдоподобие для сетки параметров
            //   migration_time x N_grid x migration_prob
            // для одного локального дерева из ARG.
            // Внутри суммирует по всем 2^n сценариям с помощью кода Грея (abacaba).
            // Возвращает DenseTensor3D: result[i][n_i][j] = лог-правдоподобие для
            // i-й migration_time, n_i-й Ne, j-й migration_prob.
            DenseTensor3D compute_grid_fast(
                const tsk_tree_t &, std::vector<double>, std::vector<double>,
                std::vector<std::vector<double>> = {});

            // Вычисляет суммарное log-правдоподобие для одного дерева
            // при фиксированных параметрах (не по сетке).
            // Суммирует по сценариям, используя smart_llh на каждом временном интервале.
            double compute_llh(const tsk_tree_t&);

            // Устанавливает параметры модели из вектора из 9 чисел.
            void set_parameters(std::span<const double>);

            // Возвращает lower_nodes_split_1_ (для отладки).
            std::vector<tsk_id_t> get_lower_nodes_split_1();

            // Выводит lower_nodes_split_1_ в stdout (для отладки).
            void print_lower_nodes_split_1();
    };

    // =========================================================================
    // lower_node -- вспомогательная свободная функция:
    // возвращает список ID узлов дерева, лежащих непосредственно ниже split_time.
    // =========================================================================
    std::vector<int> lower_node(const tsk_tree_t&, double split_time);
}
