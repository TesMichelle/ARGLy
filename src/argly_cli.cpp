#include <treellh.h>

#include <algorithm>
#include <err.h>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
    tsk_id_t sample_id = TSK_NULL;
    std::vector<double> grid_time = {500};
    std::vector<double> grid_prob = {0.2};
    std::vector<std::vector<double>> grid_N;
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
        << "  " << program << " --mode ratio --input FILE --sample-population LIST --params LIST --sample-id ID [--output FILE]\n"
        << "  " << program << " --mode grid --input FILE --sample-population LIST --params LIST --grid-time SPEC --grid-prob SPEC --output FILE\n"
        << "  " << program << " --mode grid --path DIR --prefix PREFIX --sim-count N [--suffix .arg] ...\n\n"
        << "  ratio also accepts --path DIR --prefix PREFIX --sim-count N [--suffix .arg]\n\n"
        << "Parameters:\n"
        << "  --params mt,split2,split1,mp,split2p,split1p,Nbase,Nghost,Noutgroup\n"
        << "  --sample-population comma-separated population id per sample node id\n"
        << "  --sample-id tskit node ID of the sample to classify (ratio mode)\n"
        << "  ratio output: one TSV row per tree: tree_index, left, right, log_base, log_ghost, log_ratio\n"
        << "  log_ratio = log_base - log_ghost; log_base/log_ghost are log likelihoods\n"
        << "  --grid-time and --grid-prob accept either a,b,n or comma-separated values\n"
        << "  --grid-N-base, --grid-N-ghost and --grid-N-outgroup add an N grid axis\n"
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
        } else if (arg == "--sample-id") {
            std::string value = require_value(i, argc, argv);
            size_t parsed = 0;
            long id = std::stol(value, &parsed);
            if (parsed != value.size() || id < 0 || id > std::numeric_limits<tsk_id_t>::max()) {
                throw std::runtime_error("--sample-id must be a nonnegative node ID");
            }
            cfg.sample_id = static_cast<tsk_id_t>(id);
        } else if (arg == "--grid-time") {
            cfg.grid_time = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-prob") {
            cfg.grid_prob = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-N-base") {
            if (cfg.grid_N.empty()) {
                cfg.grid_N = {{}, {}, {}};
            }
            cfg.grid_N[0] = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-N-ghost") {
            if (cfg.grid_N.empty()) {
                cfg.grid_N = {{}, {}, {}};
            }
            cfg.grid_N[1] = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-N-outgroup") {
            if (cfg.grid_N.empty()) {
                cfg.grid_N = {{}, {}, {}};
            }
            cfg.grid_N[2] = parse_grid_spec(require_value(i, argc, argv));
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

    if (cfg.mode != "llh" && cfg.mode != "grid" && cfg.mode != "ratio") {
        throw std::runtime_error("--mode must be llh, grid, or ratio");
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
    if (cfg.mode == "ratio" && cfg.sample_id == TSK_NULL) {
        throw std::runtime_error("--sample-id is required in ratio mode");
    }
    if (cfg.tree_num <= 0) {
        throw std::runtime_error("--tree-num must be positive");
    }
    if (cfg.mode == "grid" && (cfg.grid_time.empty() || cfg.grid_prob.empty())) {
        throw std::runtime_error("--grid-time and --grid-prob must be non-empty");
    }
    if (cfg.grid_N.empty()) {
        cfg.grid_N = {{cfg.parameters[6]}, {cfg.parameters[7]}, {cfg.parameters[8]}};
    } else {
        if (cfg.grid_N[0].empty()) {
            cfg.grid_N[0] = {cfg.parameters[6]};
        }
        if (cfg.grid_N[1].empty()) {
            cfg.grid_N[1] = {cfg.parameters[7]};
        }
        if (cfg.grid_N[2].empty()) {
            cfg.grid_N[2] = {cfg.parameters[8]};
        }
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
    std::string ext = mode == "grid" ? (binary ? ".llh.bin" : ".llh.txt")
        : (mode == "ratio" ? ".ratio.tsv" : ".llh.txt");
    return path.replace_extension(ext).string();
}

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

void save_matrix_text(const std::string &filename, const treellh::DenseTensor3D &tensor)
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    for (size_t time_i = 0; time_i < tensor.dim0(); ++time_i) {
        for (size_t n_i = 0; n_i < tensor.dim1(); ++n_i) {
            for (size_t prob_i = 0; prob_i < tensor.dim2(); ++prob_i) {
                if (prob_i != 0) {
                    file << ' ';
                }
                file << tensor(time_i, n_i, prob_i);
            }
            file << '\n';
        }
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

void compute_ratio_file(const Config &cfg, const std::string &input, std::ostream &out)
{
    tsk_treeseq_t ts = load_ts(input);
    const tsk_size_t node_count = tsk_treeseq_get_num_nodes(&ts);
    if (static_cast<tsk_size_t>(cfg.sample_id) >= node_count
        || !(ts.tables->nodes.flags[cfg.sample_id] & TSK_NODE_IS_SAMPLE)) {
        tsk_treeseq_free(&ts);
        throw std::runtime_error("--sample-id must identify a sample node in " + input);
    }
    const tsk_id_t *samples = tsk_treeseq_get_samples(&ts);
    const tsk_size_t sample_count = tsk_treeseq_get_num_samples(&ts);
    tsk_id_t max_sample_id = cfg.sample_id;
    for (tsk_size_t i = 0; i < sample_count; ++i) {
        max_sample_id = std::max(max_sample_id, samples[i]);
    }
    if (cfg.sample_population.size() <= static_cast<size_t>(max_sample_id)) {
        tsk_treeseq_free(&ts);
        throw std::runtime_error("--sample-population needs an entry for every sample node ID in " + input);
    }

    treellh::Scenario_Computer computer(
        ts, cfg.parameters, cfg.sample_population, cfg.admixed_index, cfg.outgroup_index);
    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);

    out << "tree_index\tleft\tright\tlog_base\tlog_ghost\tlog_ratio\n";
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    size_t tree_index = 0;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        const auto [log_base, log_ghost] = computer.compute_llh_ratio(
            tree, cfg.parameters, cfg.sample_id);
        double log_ratio = log_base - log_ghost;
        if (std::isinf(log_base) && std::isinf(log_ghost)
            && log_base < 0 && log_ghost < 0) {
            log_ratio = std::numeric_limits<double>::quiet_NaN();
        }
        out << tree_index++ << '\t' << tree.interval.left << '\t' << tree.interval.right
            << '\t' << log_base << '\t' << log_ghost << '\t' << log_ratio << '\n';
    }
    check_tsk_error(ret);
    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
}

treellh::DenseTensor3D compute_grid_file(const Config &cfg, const std::string &input)
{
    tsk_treeseq_t ts = load_ts(input);
    treellh::Scenario_Computer computer(ts, cfg.parameters, cfg.sample_population, cfg.admixed_index, cfg.outgroup_index);

    tsk_tree_t tree;
    int ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);

    size_t N_grid_size = 1;
    for (const auto &values : cfg.grid_N) {
        N_grid_size = std::max(N_grid_size, values.size());
    }
    treellh::DenseTensor3D result(cfg.grid_time.size(), N_grid_size, cfg.grid_prob.size(), 0);
    double seq_len = tsk_treeseq_get_sequence_length(&ts);
    double step = seq_len / static_cast<double>(cfg.tree_num);
    double left_x = -step - 1;
    long tree_c = 0;

    computer.set_parameters(cfg.parameters);
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        if (tree.interval.left - left_x > step) {
            print_progress(tree_c++, cfg.tree_num);
            auto temp = computer.compute_grid_fast(tree, cfg.grid_time, cfg.grid_prob, cfg.grid_N);
            for (size_t time_i = 0; time_i < result.dim0(); ++time_i) {
                for (size_t n_i = 0; n_i < result.dim1(); ++n_i) {
                    for (size_t prob_i = 0; prob_i < result.dim2(); ++prob_i) {
                        result(time_i, n_i, prob_i) += temp(time_i, n_i, prob_i);
                    }
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
    } else if (cfg.mode == "ratio") {
        if (output.empty()) {
            compute_ratio_file(cfg, input, std::cout);
        } else {
            std::ofstream file(output);
            if (!file) {
                throw std::runtime_error("cannot open output file: " + output);
            }
            compute_ratio_file(cfg, input, file);
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
                    + (cfg.mode == "grid" ? (cfg.binary_output ? "_llh.bin" : "_llh.txt")
                        : (cfg.mode == "ratio" ? "_ratio.tsv" : "_llh.txt")))).string();
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
