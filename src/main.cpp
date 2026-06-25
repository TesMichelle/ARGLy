#include <treellh.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// .npy file format
// #include <cnpy.h>

#include <err.h>
#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }

// Печать текущего прогресса выполнения (для визуализации расчета по деревьям)
void printProgress(int it, int total) {
    double progress = (double) it / total;
    int barWidth = 50;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << it + 1 << "/" << total << "\r";
    std::cout.flush();
}

// Генерация равномерной сетки параметров (аналог numpy.linspace)
std::vector<double> linspace(double a, double b, size_t n)
{
    std::vector<double> r(n, a);
    for (size_t i = 1; i < n; ++i)
    {
        r[i] += ((b - a) / (n - 1)) * i;
    }
    return r;
}

// Вывод матрицы логарифма правдоподобия в консоль
void print_llh(const std::vector<std::vector<double>> &llh)
{
    for (size_t i = 0; i < llh.size(); ++i)
    {
        for (size_t j = 0; j < llh.size(); ++j)
        {
            std::cout << llh[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

// Сохраняет тензор результатов в бинарный формат
void save_matrix_bin(const std::string &filename, const treellh::DenseTensor3D &tensor)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    file.write(
        reinterpret_cast<const char *>(tensor.data().data()),
        static_cast<std::streamsize>(tensor.size() * sizeof(double))
    );
}

// Вычисление правдоподобия на выборке ARG (последовательности деревьев) из директории path.
// Итерируется по sim_n симуляциям (файлам arg0.arg, arg1.arg и т.д.).
// Для каждой симуляции строит сетку параметров и вычисляет трехмерный тензор правдоподобия,
// после чего сохраняет результаты в бинарный файл _llh.bin.
int compute_llh_for_sample(
    std::string path,
    std::vector<int> sample_population, 
    std::vector<double> fixed_parameters,
    std::vector<double> migration_time,
    std::vector<double> migratioh_prob,
    size_t sim_n = 100,
    long tree_num = 10000,
    const std::string prefix = "arg")
{
    for (size_t sim_i = 0; sim_i < sim_n; sim_i++)
    {
        std::cout << "sim: " << sim_i << std::endl;
        // Загружаю ARG
        int ret;
        tsk_treeseq_t ts;
        ret = tsk_treeseq_load(&ts, (path + prefix + std::to_string(sim_i) + ".arg").c_str(), 0);
        if (ret != 0) {
            fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
            exit(EXIT_FAILURE);
        }

        // Создаю Scenario_Computer для текущего ARG
        treellh::Scenario_Computer computer(ts, fixed_parameters, sample_population, 0, 1);
        tsk_tree_t tree;
        computer.set_parameters(fixed_parameters);

        // Инициализация структуры дерева для прохода по локальным деревьям
        ret = tsk_tree_init(&tree, &ts, 0);
        check_tsk_error(ret);

        // Тензор для результатов расчета на сетке
        treellh::DenseTensor3D result(migration_time.size(), 1, migratioh_prob.size(), 0);
        treellh::DenseTensor3D temp;
        
        // Рабочий цикл по локальным деревьям вдоль хромосомы
        double seq_len = tsk_treeseq_get_sequence_length(&ts);
        double step = seq_len / double(tree_num);
        double left_x = 0;
        int tree_c = 0;
        
        // Перебираю деревья с шагом step для прореживания
        left_x = - step - 1;
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
        {
            if (tree.interval.left - left_x > step)
            {
                printProgress(tree_c++, tree_num);
                // Вычисление правдоподобия на сетке для текущего дерева
                temp = computer.compute_grid_fast(tree, migration_time, migratioh_prob);
                std::cout << "hello" << std::endl;
                for (size_t i = 0; i < migration_time.size(); i++)
                {
                    for (size_t j = 0; j < migratioh_prob.size(); ++j)
                    {            
                        result(i, 0, j) += temp(i, 0, j);
                    }
                }
                left_x = tree.interval.left;
            }
        }
        std::cout << std::endl;
        // Сохранение итоговой матрицы
        save_matrix_bin(path + prefix + std::to_string(sim_i) + "_llh.bin", result);
    }
    return 0;
}

int main()
{ 
    // Начальные фиксированные параметры модели
    // 0: migration_time, 1: split_2_time, 2: split_1_time
    // 3: migration_prop, 4: split_2_prop, 5: split_1_prop
    // 6: N_base, 7: N_ghost, 8: N_outgroup
    std::vector<double> parameters = {
        500, 1000, 2000,
        0.2, 0.5, 0.5,
        10000, 10000, 10000
    };

    size_t migration_time_size = 1;
    size_t migration_prob_size = 1;
    std::vector<double> migration_time = linspace(500, 500, migration_time_size);
    std::vector<double> migration_prob = linspace(0.2, 0.2, migration_prob_size);

    long tree_num = 1;

    // Загружаю ARG
    int ret;
    tsk_treeseq_t ts;
    ret = tsk_treeseq_load(&ts, "my_test.trees", 0);
    if (ret != 0) {
        fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
        exit(EXIT_FAILURE);
    }

    // Задаем популяцию для 100 образцов (все 0 - AFR по умолчанию)
    std::vector<int> sample_population(100, 0);
    for (size_t i = 0; i < 100; ++i)
    {
        sample_population[i] = 0;
    }

    tsk_treeseq_free(&ts);

    compute_llh_for_sample("arg_files/", sample_population, parameters, migration_time, migration_prob, 100, 1);

    return 0;
}
