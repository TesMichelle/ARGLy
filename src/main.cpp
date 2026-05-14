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
    for (size_t i = 1; i < n; ++i)
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

void save_matrix_bin(std::string filename, const treellh::DenseTensor3D &mtx)
{
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(mtx.data().data()), mtx.size() * sizeof(double));
    file.close();
}


std::vector<double> compute_for_sample(
    std::string path, 
    std::vector<int> sample_population,
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
        treellh::Scenario_Computer computer(ts, parameters, sample_population, 0, 1);
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

// std::vector<double> estimate_for_sample(
//     std::string path, 
//     std::vector<double> parameters, 
//     size_t n = 30,
//     size_t sim_n = 100,
//     long tree_num = 10000,
//     const std::string prefix = "real_arg_admixed_only_50_sim_")
// {


//     // создаю сетку
//     std::vector<double> time_migration = linspace(1, 7000, n);

//     std::vector<double> pe(sim_n, 0);
//     for (size_t sim_i = 0; sim_i < sim_n; sim_i++)
//     {
//         std::cout << "sim: " << sim_i << std::endl;
//         // Загружаю АРГ 
//         int ret;
//         tsk_treeseq_t ts;
//         ret = tsk_treeseq_load(&ts, (path + prefix + std::to_string(sim_i) + ".arg").c_str(), 0);
//         if (ret != 0) {
//             fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
//             exit(EXIT_FAILURE);
//         }

//         // создраю вычислятор
//         treellh::Scenario_Computer computer(ts, parameters, 3, 1);
//         tsk_tree_t tree;

//         // иницилизация дерева для нужного трисека
//         ret = tsk_tree_init(&tree, &ts, 0);
//         check_tsk_error(ret);

//         // std::vector<double> mig_prob = linspace(0.01, 0.01, n);

//         // правдоподобие на сетке
//         std::vector<std::vector<double>> llh(n, std::vector<double>(n, 0));

//         // рабочий цикл
//         double seq_len = tsk_treeseq_get_sequence_length(&ts);
//         double step = seq_len / double(tree_num);
//         double left_x = 0;
//         size_t max_i = 0;
//         for (size_t i = 0; i < n; ++i)
//         {
//             for (size_t j = 0; j < n; ++j)
//             {
//                 std::cout << i << ", " << j << std::endl;
//                 parameters[0] = time_migration[i];
//                 computer.set_parameters(parameters);

//                 // перебираю деревья
//                 left_x = - step - 1;
//                 for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
//                 {
//                     if (tree.interval.left - left_x > step)
//                     {
//                         llh[i][j] += computer.compute_llh(tree);
//                         left_x = tree.interval.left;
//                     }
//                 }
//                 std::cout << std::endl;
//                 if (llh[i][j] > llh[max_i][j])
//                 {
//                     max_i = i;
//                 }
//                 std::cout << llh[i][j] << "max: " << llh[max_i][j] << std::endl;
//                 break;
//             }
//         }
//         pe[sim_i] = time_migration[max_i];
//         std::cout << "Point estimate: " << pe[sim_i] << std::endl;
//     }
//     return pe;
// }

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
        // Загружаю АРГ 
        int ret;
        tsk_treeseq_t ts;
        ret = tsk_treeseq_load(&ts, (path + prefix + std::to_string(sim_i) + ".arg").c_str(), 0);
        if (ret != 0) {
            fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
            exit(EXIT_FAILURE);
        }

        // создраю вычислятор
        treellh::Scenario_Computer computer(ts, fixed_parameters, sample_population, 0, 1);
        tsk_tree_t tree;
        computer.set_parameters(fixed_parameters);

        // иницилизация дерева для нужного трисека
        ret = tsk_tree_init(&tree, &ts, 0);
        check_tsk_error(ret);

        // std::vector<double> mig_prob = linspace(0.01, 0.01, n);

        // правдоподобие сюда
        treellh::DenseTensor3D result(migration_time.size(), 1, migratioh_prob.size(), 0);
        treellh::DenseTensor3D temp;
        // рабочий цикл
        double seq_len = tsk_treeseq_get_sequence_length(&ts);
        double step = seq_len / double(tree_num);
        double left_x = 0;
        int tree_c = 0;
        // перебираю деревья
        left_x = - step - 1;
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
        {
            if (tree.interval.left - left_x > step)
            {
                printProgress(tree_c++, tree_num);
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
        save_matrix_bin(path + prefix + std::to_string(sim_i) + "_llh.bin", result);
    }
    return 0;
}

int main()
{ 
    // Параметры 
    // miration_time, split_2_time, split_1_time 
    // migration_prop, split_2_prop, split_1_prop
    // base_size, outgroup_size, ghost_size
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

    // Загружаю АРГ 
    int ret;
    tsk_treeseq_t ts;
    ret = tsk_treeseq_load(&ts, "my_test.trees", 0);
    // ret = tsk_treeseq_load(&ts, "RELATE_arg/real_arg_30_sim_0.relate.ts.trees", 0);
    if (ret != 0) {
        fprintf(stderr, "Load error:%s\n", tsk_strerror(ret));
        exit(EXIT_FAILURE);
    }

    std::vector<int> sample_population(100, 0);
    for (size_t i = 0; i < 100; ++i)
    {
        sample_population[i] = 0;
    }
    // создраю вычислятор
    // treellh::Scenario_Computer computer(ts, parameters, sample_population, 1, 0);
    // tsk_tree_t tree;
    // computer.set_parameters(parameters);

    // // иницилизация дерева для нужного трисека
    // ret = tsk_tree_init(&tree, &ts, 0);
    // check_tsk_error(ret);

    // // samples
    // const tsk_id_t *samples = tsk_treeseq_get_samples(&ts);
    // const tsk_size_t samples_num = tsk_treeseq_get_num_samples(&ts);

    // std::vector<tsk_id_t> samples_out(samples, samples + 8);
    // std::vector<tsk_id_t> samples_adm(samples + 8, samples + 16); 

    // // правдоподобие сюда
    // std::vector<std::vector<double>> result(migration_time.size(), std::vector<double>(migration_prob.size(), 0));
    // std::vector<std::vector<double>> temp;
    // // рабочий цикл
    // double seq_len = tsk_treeseq_get_sequence_length(&ts);
    // double step = seq_len / double(tree_num);
    // double left_x = 0;
    // int tree_c = 0;
    // // перебираю деревья
    // left_x = - step - 1;
    // for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) 
    // {
    //     if (tree.interval.left - left_x > step)
    //     {
    //         printProgress(tree_c++, tree_num);
    //         // std::cout << tree.index << std::endl;
    //         // parameters[1] = computer.find_lowest_coal(tree, samples_out, samples_adm);
    //         computer.set_parameters(parameters);
    //         temp = computer.compute_grid_fast(tree, migration_time, migration_prob);
    //         for (size_t i = 0; i < migration_time.size(); i++)
    //         {
    //             for (size_t j = 0; j < migration_prob.size(); ++j)
    //             {            
    //                 result[i][j] += temp[i][j];
    //             }
    //         }
    //         left_x = tree.interval.left;
    //     }
    // }
    // std::cout << std::endl;
    // std::cout << result[0][0] << '\n';
    // save_matrix_bin("my_test.llh.bin", result);

    compute_llh_for_sample("arg_files/", sample_population, parameters, migration_time, migration_prob, 100, 1);

    return 0;
}
