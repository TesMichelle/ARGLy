#include <treellh.h>

#include <err.h>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }

namespace {

struct Config {
    std::string mode;
    std::string input;
    std::string path;
    std::string prefix = "arg";
    std::string suffix = ".arg";
    std::string output;
    std::vector<double> parameters = {
        500, 1000, 2000,
        0.2, 0.5, 0.5,
        10000, 10000, 10000
    };
    std::vector<int> sample_population;
    std::vector<double> grid_time = {500};
    std::vector<double> grid_prob = {0.2};
    int admixed_index = 0;
    int outgroup_index = 1;
    size_t sim_count = 0;
    long tree_num = 10000;
    bool binary_output = true;
};

void print_usage(const char *program)
{
    std::cerr
        << "Usage:\n"
        << "  " << program << " --mode llh  --input FILE --sample-population LIST --params LIST [--output FILE]\n"
        << "  " << program << " --mode grid --input FILE --sample-population LIST --params LIST --grid-time SPEC --grid-prob SPEC --output FILE\n"
        << "  " << program << " --mode grid --path DIR --prefix PREFIX --sim-count N [--suffix .arg] ...\n\n"
        << "Parameters:\n"
        << "  --params mt,split2,split1,mp,split2p,split1p,Nbase,Nghost,Noutgroup\n"
        << "  --sample-population comma-separated population id per sample node id\n"
        << "  --grid-time and --grid-prob accept either a,b,n or comma-separated values\n"
        << "  --tree-num K samples approximately K trees along the sequence\n"
        << "  --text-output writes grid as whitespace text instead of raw doubles\n";
}

std::vector<std::string> split(const std::string &value, char delim)
{
    std::vector<std::string> result;
    std::stringstream ss(value);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    return result;
}

std::vector<double> parse_doubles(const std::string &value)
{
    std::vector<double> result;
    for (const auto &item : split(value, ',')) {
        result.push_back(std::stod(item));
    }
    return result;
}

std::vector<int> parse_ints(const std::string &value)
{
    std::vector<int> result;
    for (const auto &item : split(value, ',')) {
        result.push_back(std::stoi(item));
    }
    return result;
}

std::vector<double> linspace(double a, double b, size_t n)
{
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return {a};
    }
    std::vector<double> result(n, a);
    for (size_t i = 1; i < n; ++i) {
        result[i] += ((b - a) / static_cast<double>(n - 1)) * static_cast<double>(i);
    }
    return result;
}

std::vector<double> parse_grid_spec(const std::string &value)
{
    std::vector<std::string> fields = split(value, ',');
    if (fields.size() == 3) {
        char *end = nullptr;
        long n = std::strtol(fields[2].c_str(), &end, 10);
        if (end != fields[2].c_str() && *end == '\0' && n > 0) {
            return linspace(std::stod(fields[0]), std::stod(fields[1]), static_cast<size_t>(n));
        }
    }
    return parse_doubles(value);
}

std::string require_value(int &i, int argc, char **argv)
{
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + argv[i]);
    }
    return argv[++i];
}

Config parse_args(int argc, char **argv)
{
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else if (arg == "--mode") {
            cfg.mode = require_value(i, argc, argv);
        } else if (arg == "--input") {
            cfg.input = require_value(i, argc, argv);
        } else if (arg == "--path") {
            cfg.path = require_value(i, argc, argv);
        } else if (arg == "--prefix") {
            cfg.prefix = require_value(i, argc, argv);
        } else if (arg == "--suffix") {
            cfg.suffix = require_value(i, argc, argv);
        } else if (arg == "--output") {
            cfg.output = require_value(i, argc, argv);
        } else if (arg == "--params") {
            cfg.parameters = parse_doubles(require_value(i, argc, argv));
        } else if (arg == "--sample-population") {
            cfg.sample_population = parse_ints(require_value(i, argc, argv));
        } else if (arg == "--grid-time") {
            cfg.grid_time = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-prob") {
            cfg.grid_prob = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--admixed-index") {
            cfg.admixed_index = std::stoi(require_value(i, argc, argv));
        } else if (arg == "--outgroup-index") {
            cfg.outgroup_index = std::stoi(require_value(i, argc, argv));
        } else if (arg == "--sim-count") {
            cfg.sim_count = static_cast<size_t>(std::stoul(require_value(i, argc, argv)));
        } else if (arg == "--tree-num") {
            cfg.tree_num = std::stol(require_value(i, argc, argv));
        } else if (arg == "--text-output") {
            cfg.binary_output = false;
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    if (cfg.mode != "llh" && cfg.mode != "grid") {
        throw std::runtime_error("--mode must be either llh or grid");
    }
    if (cfg.input.empty() && (cfg.path.empty() || cfg.sim_count == 0)) {
        throw std::runtime_error("provide either --input FILE or --path DIR --sim-count N");
    }
    if (!cfg.input.empty() && (!cfg.path.empty() || cfg.sim_count != 0)) {
        throw std::runtime_error("--input cannot be combined with batch --path/--sim-count");
    }
    if (cfg.parameters.size() != 9) {
        throw std::runtime_error("--params must contain exactly 9 numbers");
    }
    if (cfg.sample_population.empty()) {
        throw std::runtime_error("--sample-population is required");
    }
    if (cfg.tree_num <= 0) {
        throw std::runtime_error("--tree-num must be positive");
    }
    if (cfg.mode == "grid" && (cfg.grid_time.empty() || cfg.grid_prob.empty())) {
        throw std::runtime_error("--grid-time and --grid-prob must be non-empty");
    }
    return cfg;
}

void print_progress(long it, long total)
{
    const int bar_width = 50;
    double progress = static_cast<double>(it) / static_cast<double>(total);
    int pos = static_cast<int>(bar_width * progress);

    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) {
            std::cout << "=";
        } else if (i == pos) {
            std::cout << ">";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "] " << it + 1 << "/" << total << "\r";
    std::cout.flush();
}

std::string batch_input_name(const Config &cfg, size_t sim_i)
{
    std::filesystem::path path(cfg.path);
    path /= cfg.prefix + std::to_string(sim_i) + cfg.suffix;
    return path.string();
}

std::string default_output_name(const std::string &input, const std::string &mode, bool binary)
{
    std::filesystem::path path(input);
    std::string ext = mode == "grid" ? (binary ? ".llh.bin" : ".llh.txt") : ".llh.txt";
    return path.replace_extension(ext).string();
}

void save_matrix_bin(const std::string &filename, const std::vector<std::vector<double>> &matrix)
{
    std::vector<double> flat;
    for (const auto &row : matrix) {
        flat.insert(flat.end(), row.begin(), row.end());
    }

    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    file.write(reinterpret_cast<const char *>(flat.data()), static_cast<std::streamsize>(flat.size() * sizeof(double)));
}

void save_matrix_text(const std::string &filename, const std::vector<std::vector<double>> &matrix)
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    for (const auto &row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            if (j != 0) {
                file << ' ';
            }
            file << row[j];
        }
        file << '\n';
    }
}

void save_scalar_text(const std::string &filename, double value)
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    file << value << '\n';
}

tsk_treeseq_t load_ts(const std::string &filename)
{
    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, filename.c_str(), 0);
    if (ret != 0) {
        throw std::runtime_error("load error for " + filename + ": " + tsk_strerror(ret));
    }
    return ts;
}

double compute_llh_file(const Config &cfg, const std::string &input)
{
    tsk_treeseq_t ts = load_ts(input);
    treellh::Scenario_Computer computer(ts, cfg.parameters, cfg.sample_population, cfg.admixed_index, cfg.outgroup_index);

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);

    double result = 0;
    computer.set_parameters(cfg.parameters);
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        result += computer.compute_llh(tree);
    }
    check_tsk_error(ret);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    return result;
}

std::vector<std::vector<double>> compute_grid_file(const Config &cfg, const std::string &input)
{
    tsk_treeseq_t ts = load_ts(input);
    treellh::Scenario_Computer computer(ts, cfg.parameters, cfg.sample_population, cfg.admixed_index, cfg.outgroup_index);

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);

    std::vector<std::vector<double>> result(cfg.grid_time.size(), std::vector<double>(cfg.grid_prob.size(), 0));
    double seq_len = tsk_treeseq_get_sequence_length(&ts);
    double step = seq_len / static_cast<double>(cfg.tree_num);
    double left_x = -step - 1;
    long tree_c = 0;

    computer.set_parameters(cfg.parameters);
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        if (tree.interval.left - left_x > step) {
            print_progress(tree_c++, cfg.tree_num);
            auto temp = computer.compute_grid_fast(tree, cfg.grid_time, cfg.grid_prob);
            for (size_t i = 0; i < cfg.grid_time.size(); ++i) {
                for (size_t j = 0; j < cfg.grid_prob.size(); ++j) {
                    result[i][j] += temp[i][j];
                }
            }
            left_x = tree.interval.left;
        }
    }
    check_tsk_error(ret);
    std::cout << '\n';

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    return result;
}

void run_single(const Config &cfg, const std::string &input, const std::string &output)
{
    if (cfg.mode == "llh") {
        double llh = compute_llh_file(cfg, input);
        if (output.empty()) {
            std::cout << llh << '\n';
        } else {
            save_scalar_text(output, llh);
        }
    } else {
        auto matrix = compute_grid_file(cfg, input);
        if (cfg.binary_output) {
            save_matrix_bin(output, matrix);
        } else {
            save_matrix_text(output, matrix);
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    try {
        Config cfg = parse_args(argc, argv);
        if (!cfg.input.empty()) {
            std::string output = cfg.output;
            if (output.empty() && cfg.mode == "grid") {
                output = default_output_name(cfg.input, cfg.mode, cfg.binary_output);
            }
            run_single(cfg, cfg.input, output);
            return EXIT_SUCCESS;
        }

        for (size_t sim_i = 0; sim_i < cfg.sim_count; ++sim_i) {
            std::string input = batch_input_name(cfg, sim_i);
            std::string output;
            if (!cfg.output.empty()) {
                std::filesystem::path out_dir(cfg.output);
                output = (out_dir / (cfg.prefix + std::to_string(sim_i)
                    + (cfg.mode == "grid" ? (cfg.binary_output ? "_llh.bin" : "_llh.txt") : "_llh.txt"))).string();
            } else {
                output = default_output_name(input, cfg.mode, cfg.binary_output);
            }

            std::cout << "sim: " << sim_i << " input: " << input << '\n';
            run_single(cfg, input, output);
        }
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
