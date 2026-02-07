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

std::vector<double> linspace(double a, double b, size_t n)
{
    std::vector<double> r(n, a);
    for (size_t i = 0; i < n; ++i)
    {
        r[i] += ((b - a) / (n - 1)) * i;
    }
    return r;
}

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

void print_llh(const std::vector<double> &llh)
{
    for (size_t i = 0; i < llh.size(); ++i)
    {
        std::cout << llh[i] << std::endl;
    }
}

void save_matrix_bin(std::string filename, const std::vector<std::vector<double>> mtx)
{
    std::vector<double> flat;
    for (const auto& row : mtx) {
        flat.insert(flat.end(), row.begin(), row.end());
    }

    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(flat.data()), flat.size() * sizeof(double));
    file.close();
}


std::vector<double> compute_for_sample(
    std::string path, 
    std::vector<double> parameters, 
    const std::string prefix = "test_arg_10_sim_")
{
    std::vector<double> llh(100, 0);
    for (size_t sim_i = 0; sim_i < 100; sim_i++)
    {
        std::cout << "sim: " << sim_i << std::endl;
        // Загружаю АРГ 
        int ret;
        tsk_treeseq_t ts;
        ret = tsk_treeseq_load(&ts, (path + prefix + std::to_string(sim_i) + ".arg").c_str(), 0);
        if (ret != 0) {
            fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
            exit(EXIT_FAILURE);
        }

        // создраю вычислятор
        treellh::Scenario_Computer computer(ts, parameters, 3, 1);
        tsk_tree_t tree;

        // иницилизация дерева для нужного трисека
        ret = tsk_tree_init(&tree, &ts, 0);
        check_tsk_error(ret);


        computer.set_parameters(parameters);
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
        {
            llh[sim_i] += computer.compute_llh(tree);
        }
        std::cout << llh[sim_i] << std::endl;
    }
    return llh;
}

std::vector<double> estimate_for_sample(
    std::string path, 
    std::vector<double> parameters, 
    size_t n = 30,
    size_t sim_n = 100,
    long tree_num = 10000,
    const std::string prefix = "real_arg_admixed_only_50_sim_")
{


    // создаю сетку
    std::vector<double> time_migration = linspace(1, 7000, n);

    std::vector<double> pe(sim_n, 0);
    for (size_t sim_i = 0; sim_i < sim_n; sim_i++)
    {
        std::cout << "sim: " << sim_i << std::endl;
        // Загружаю АРГ 
        int ret;
        tsk_treeseq_t ts;
        ret = tsk_treeseq_load(&ts, (path + prefix + std::to_string(sim_i) + ".arg").c_str(), 0);
        if (ret != 0) {
            fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
            exit(EXIT_FAILURE);
        }

        // создраю вычислятор
        treellh::Scenario_Computer computer(ts, parameters, 3, 1);
        tsk_tree_t tree;

        // иницилизация дерева для нужного трисека
        ret = tsk_tree_init(&tree, &ts, 0);
        check_tsk_error(ret);

        // std::vector<double> mig_prob = linspace(0.01, 0.01, n);

        // правдоподобие на сетке
        std::vector<std::vector<double>> llh(n, std::vector<double>(n, 0));

        // рабочий цикл
        int tree_c = 0;
        double seq_len = tsk_treeseq_get_sequence_length(&ts);
        double step = seq_len / double(tree_num);
        double left_x = 0;
        size_t max_i = 0;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                std::cout << i << ", " << j << std::endl;
                parameters[0] = time_migration[i];
                computer.set_parameters(parameters);

                // перебираю деревья
                left_x = - step - 1;
                tree_c = 0;
                for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
                {
                    if (tree.interval.left - left_x > step)
                    {
                        // printProgress(tree_c, tree_num);
                        llh[i][j] += computer.compute_llh(tree);
                        left_x = tree.interval.left;
                        tree_c += 1;
                    }
                }
                std::cout << std::endl;
                if (llh[i][j] > llh[max_i][j])
                {
                    max_i = i;
                }
                std::cout << llh[i][j] << "max: " << llh[max_i][j] << std::endl;
                break;
            }
        }
        pe[sim_i] = time_migration[max_i];
        std::cout << "Point estimate: " << pe[sim_i] << std::endl;
    }
    return pe;
}

int main()
{ 
    // Параметры 
    // miration_time, split_2_time, split_1_time 
    // migration_prop, split_2_prop, split_1_prop
    // base_size, outgroup_size, ghost_size
    std::vector<double> parameters = {
        5000, 7000, 7001,
        0.03, 0.5, 0.5,
        10000, 10000, 10000
    };

    std::vector<double> pe;
    pe = estimate_for_sample("sim/", parameters, 2, 1, 1, "real_arg_30_sim_");

    // std::ofstream file("real_pe_30.bin", std::ios::binary);
    // file.write(reinterpret_cast<const char*>(pe.data()), pe.size() * sizeof(double));
    // file.close();

    // save_matrix_bin("llh.bin", llh);
    
    return 0;
}