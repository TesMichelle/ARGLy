#include <treellh.h>
#include <likelihood_cache.h>
#include "fit_search.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <err.h>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/resource.h>
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
    bool grid_time_specified = false;
    bool grid_prob_specified = false;
    bool grid_requested = false;
    std::vector<std::vector<double>> grid_N;
    int admixed_index = 0;
    int outgroup_index = 1;
    size_t sim_count = 0;
    long tree_num = 10000;
    bool tree_num_specified = false;
    bool binary_output = true;
    std::string queries;
    bool estimate_only = false;
    bool per_tree = false;
    double cache_limit_mb = 0;
    std::array<std::vector<double>, 4> fit_axes;
    size_t fit_maxeval = 250;
    size_t fit_starts = 3;
    std::string tree_cache;
};

double peak_rss_mib()
{
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
    return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
}

class ProgressPrinter {
public:
    void print(const std::string& label, size_t done, size_t total, bool complete = false,
        double best = -std::numeric_limits<double>::infinity())
    {
        if (label != current_) {
            current_ = label;
            last_bucket_ = -1;
            start_ = std::chrono::steady_clock::now();
        }
        const int percent = total ? static_cast<int>(100.0 * done / total) : 100;
        const int bucket = complete ? 100 : (percent / 10) * 10;
        if (!complete && bucket <= last_bucket_) return;
        if (complete && last_bucket_ == 100) return;
        last_bucket_ = bucket;
        std::ostringstream line;
        line << label << ": " << bucket << "% (" << done << "/" << total << ")"
             << ", elapsed " << std::fixed << std::setprecision(2)
             << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count() << " s";
        if (std::isfinite(best)) line << ", best log_llh " << std::setprecision(6) << best;
        if (complete && done < total) line << ", stopped before evaluation limit";
        std::cerr << line.str() << '\n';
        std::cerr.flush();
    }
private:
    std::string current_;
    int last_bucket_ = -1;
    std::chrono::steady_clock::time_point start_{};
};

void print_usage(const char *program)
{
    std::cerr
        << "Usage:\n"
        << "  " << program << " --mode llh  --input FILE --sample-population LIST --params LIST [--output FILE]\n"
        << "  " << program << " --mode ratio --input FILE --sample-population LIST --params LIST --sample-id ID [--output FILE]\n"
        << "  " << program << " --mode ratio --input FILE --sample-population LIST --params LIST --sample-id ID --grid-time SPEC --grid-prob SPEC [--grid-N-* SPEC] [--output FILE]\n"
        << "  " << program << " --mode cached --input FILE --sample-population LIST --params LIST [--queries FILE|-] [--estimate-only] [--output FILE]\n"
        << "  " << program << " --mode grid --input FILE --sample-population LIST --params LIST --grid-time SPEC --grid-prob SPEC --output FILE\n"
        << "  " << program << " --mode grid --path DIR --prefix PREFIX --sim-count N [--suffix .arg] ...\n\n"
        << "  " << program << " --mode fit-grid|fit-opt|fit-hybrid --input FILE --sample-population LIST --params LIST --fit-time LO,HI,COUNT --fit-prob LO,HI,COUNT --fit-N-base LO,HI,COUNT --fit-N-ghost LO,HI,COUNT [--output FILE]\n"
        << "  ratio also accepts --path DIR --prefix PREFIX --sim-count N [--suffix .arg]\n\n"
        << "Parameters:\n"
        << "  --params mt,split2,split1,mp,split2p,split1p,Nbase,Nghost,Noutgroup\n"
        << "  --sample-population comma-separated population id per sample node id\n"
        << "  --sample-id tskit node ID of the sample to classify (ratio mode)\n"
        << "  ratio output: one TSV row per tree; with grid options, one row per tree and grid point\n"
        << "  log_ratio = log_base - log_ghost; log_base/log_ghost are log likelihoods\n"
        << "  --grid-time and --grid-prob accept either a,b,n or comma-separated values\n"
        << "  --grid-N-base, --grid-N-ghost and --grid-N-outgroup add an N grid axis\n"
        << "  --tree-num K: fit/cached sample K evenly spaced genomic positions (unique trees); grid samples approximately K along the sequence\n"
        << "  --text-output writes grid as whitespace text instead of raw doubles\n";
    std::cerr
        << "  cached: builds statistics for selected trees once, then evaluates arbitrary t,p,N queries\n"
        << "  --queries FILE|-: nine comma/space-separated parameters per line; '-' reads stdin\n"
        << "  split times must stay fixed; requires 0 <= mt < split2 < split1\n"
        << "  --estimate-only: report cache bytes without building the scenario cache\n"
        << "  --cache-limit-mb MB: refuse a cache larger than this estimate (MiB; 0 means unlimited)\n"
        << "  --tree-cache FILE: save/load selected tree scenario statistics between processes; cached mode opts in, fit modes use an automatic path\n"
        << "  --per-tree: cached output per tree instead of the sum of tree log likelihoods\n"
        << "  cached memory/build/query timings go to stderr; results go to stdout or --output\n";
    std::cerr
        << "  fit modes maximize the sum of tree log likelihoods using the scenario cache\n"
        << "  fit-grid scans the Cartesian grid; fit-opt starts NLopt BOBYQA at --params; fit-hybrid adds local searches from separated grid points and --params\n"
        << "  --fit-maxeval K: objective calls per local optimization (default 250)\n"
        << "  --fit-starts K: number of hybrid starts (default 3)\n"
        << "  fit axes contain inclusive lower/upper bounds and grid count; N axes use logarithmic spacing\n"
        << "  split times and Noutgroup stay fixed; the score ignores split proportions and Noutgroup\n";
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
        } else if (arg == "--queries") {
            cfg.queries = require_value(i, argc, argv);
        } else if (arg == "--estimate-only") {
            cfg.estimate_only = true;
        } else if (arg == "--per-tree") {
            cfg.per_tree = true;
        } else if (arg == "--cache-limit-mb") {
            cfg.cache_limit_mb = std::stod(require_value(i, argc, argv));
        } else if (arg == "--fit-time") {
            cfg.fit_axes[0] = parse_doubles(require_value(i, argc, argv));
        } else if (arg == "--fit-prob") {
            cfg.fit_axes[1] = parse_doubles(require_value(i, argc, argv));
        } else if (arg == "--fit-N-base") {
            cfg.fit_axes[2] = parse_doubles(require_value(i, argc, argv));
        } else if (arg == "--fit-N-ghost") {
            cfg.fit_axes[3] = parse_doubles(require_value(i, argc, argv));
        } else if (arg == "--fit-maxeval") {
            cfg.fit_maxeval = std::stoul(require_value(i, argc, argv));
        } else if (arg == "--fit-starts") {
            cfg.fit_starts = std::stoul(require_value(i, argc, argv));
        } else if (arg == "--tree-cache") {
            cfg.tree_cache = require_value(i, argc, argv);
        } else if (arg == "--grid-time") {
            cfg.grid_time = parse_grid_spec(require_value(i, argc, argv));
            cfg.grid_time_specified = cfg.grid_requested = true;
        } else if (arg == "--grid-prob") {
            cfg.grid_prob = parse_grid_spec(require_value(i, argc, argv));
            cfg.grid_prob_specified = cfg.grid_requested = true;
        } else if (arg == "--grid-N-base") {
            cfg.grid_requested = true;
            if (cfg.grid_N.empty()) {
                cfg.grid_N = {{}, {}, {}};
            }
            cfg.grid_N[0] = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-N-ghost") {
            cfg.grid_requested = true;
            if (cfg.grid_N.empty()) {
                cfg.grid_N = {{}, {}, {}};
            }
            cfg.grid_N[1] = parse_grid_spec(require_value(i, argc, argv));
        } else if (arg == "--grid-N-outgroup") {
            cfg.grid_requested = true;
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
            cfg.tree_num_specified = true;
        } else if (arg == "--text-output") {
            cfg.binary_output = false;
        } else {
            throw std::runtime_error("unknown argument: " + arg);
        }
    }

    const bool fitting = cfg.mode == "fit-grid" || cfg.mode == "fit-opt" || cfg.mode == "fit-hybrid";
    if (cfg.mode != "llh" && cfg.mode != "grid" && cfg.mode != "ratio" && cfg.mode != "cached" && !fitting) {
        throw std::runtime_error("--mode must be llh, grid, ratio, cached, fit-grid, fit-opt, or fit-hybrid");
    }
    if ((cfg.mode == "cached" || fitting) && cfg.input.empty())
        throw std::runtime_error("cached and fit modes require --input FILE (one tree sequence ensemble)");
    if (!std::isfinite(cfg.cache_limit_mb) || cfg.cache_limit_mb < 0)
        throw std::runtime_error("--cache-limit-mb must be finite and nonnegative");
    if (cfg.mode != "cached" && (!cfg.queries.empty() || cfg.per_tree))
        throw std::runtime_error("--queries and --per-tree require --mode cached");
    if (cfg.mode != "cached" && !fitting && (cfg.estimate_only || cfg.cache_limit_mb != 0))
        throw std::runtime_error("cache options require cached or fit mode");
    if (fitting) {
        for (const auto& axis : cfg.fit_axes)
            if (axis.size() != 3 || !std::isfinite(axis[2]) || axis[2] < 1 || axis[2] > 1000000 || std::floor(axis[2]) != axis[2])
                throw std::runtime_error("each --fit axis requires LO,HI,positive integer COUNT");
    } else {
        for (const auto& axis : cfg.fit_axes)
            if (!axis.empty()) throw std::runtime_error("--fit axes require a fit mode");
        if (cfg.mode != "cached" && !cfg.tree_cache.empty())
            throw std::runtime_error("--tree-cache requires cached or fit mode");
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
    if (cfg.mode == "ratio" && cfg.grid_requested) {
        if (!cfg.grid_time_specified) cfg.grid_time = {cfg.parameters[0]};
        if (!cfg.grid_prob_specified) cfg.grid_prob = {cfg.parameters[3]};
        if (cfg.grid_time.empty() || cfg.grid_prob.empty())
            throw std::runtime_error("ratio grid axes must be nonempty");
        for (size_t i = 0; i < cfg.grid_time.size(); ++i)
            if (!std::isfinite(cfg.grid_time[i]) || cfg.grid_time[i] < 0
                || cfg.grid_time[i] >= cfg.parameters[1]
                || (i && cfg.grid_time[i] > cfg.grid_time[i - 1]))
                throw std::runtime_error("ratio --grid-time must descend and stay in [0, split2)");
        for (double p : cfg.grid_prob)
            if (!std::isfinite(p) || p < 0 || p > 1)
                throw std::runtime_error("ratio --grid-prob must be in [0,1]");
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
    if (cfg.mode == "ratio" && cfg.grid_requested) {
        size_t width = 1;
        for (const auto& axis : cfg.grid_N) width = std::max(width, axis.size());
        for (const auto& axis : cfg.grid_N) {
            if (axis.size() != 1 && axis.size() != width)
                throw std::runtime_error("ratio N axes must have length 1 or a common grid length");
            for (double n : axis)
                if (!std::isfinite(n) || n <= 0)
                    throw std::runtime_error("ratio N grid values must be finite and positive");
        }
    }
    return cfg;
}

std::string tree_cache_signature(const Config& cfg, const std::string& input)
{
    uint64_t hash = 14695981039346656037ULL;
    auto byte = [&](uint8_t x) { hash = (hash ^ x) * 1099511628211ULL; };
    auto word = [&](uint64_t x) {
        for (int i = 0; i < 8; ++i) { byte(static_cast<uint8_t>(x)); x >>= 8; }
    };
    for (char x : std::string("ARGLy-tree-statistics-v3")) byte(static_cast<uint8_t>(x));
    std::ifstream file(input, std::ios::binary);
    if (!file) throw std::runtime_error("cannot hash input for fit cache: " + input);
    std::array<char, 65536> buffer{};
    uint64_t length = 0;
    while (file) {
        file.read(buffer.data(), buffer.size());
        const auto count = file.gcount();
        length += static_cast<uint64_t>(count);
        for (std::streamsize i = 0; i < count; ++i)
            byte(static_cast<uint8_t>(buffer[static_cast<size_t>(i)]));
    }
    if (!file.eof()) throw std::runtime_error("failed hashing input for fit cache: " + input);
    word(length);
    for (size_t i : {size_t{1}, size_t{2}, size_t{4}, size_t{5}, size_t{8}})
        word(std::bit_cast<uint64_t>(cfg.parameters[i]));
    word(static_cast<uint64_t>(cfg.admixed_index));
    word(static_cast<uint64_t>(cfg.outgroup_index));
    word(cfg.sample_population.size());
    for (int population : cfg.sample_population) word(static_cast<uint64_t>(population));
    if (cfg.tree_num_specified) {
        for (char x : std::string("genomic-midpoints-v1")) byte(static_cast<uint8_t>(x));
        word(static_cast<uint64_t>(cfg.tree_num));
    }
    std::ostringstream signature;
    signature << std::hex << std::setw(16) << std::setfill('0') << hash;
    return signature.str();
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

    if (cfg.grid_requested)
        out << "tree_index\tleft\tright\ttime_index\tN_index\tprob_index"
            << "\tmigration_time\tNbase\tNghost\tNoutgroup\tmigration_prob"
            << "\tlog_base\tlog_ghost\tlog_ratio\n";
    else
        out << "tree_index\tleft\tright\tlog_base\tlog_ghost\tlog_ratio\n";
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    size_t tree_index = 0;
    size_t n_grid_size = 1;
    for (const auto& axis : cfg.grid_N) n_grid_size = std::max(n_grid_size, axis.size());
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)) {
        if (cfg.grid_requested) {
            auto [base, ghost] = computer.compute_grid_ratio_fast(
                tree, cfg.grid_time, cfg.grid_prob, cfg.grid_N, cfg.sample_id);
            for (size_t t = 0; t < cfg.grid_time.size(); ++t)
                for (size_t n = 0; n < n_grid_size; ++n)
                    for (size_t p = 0; p < cfg.grid_prob.size(); ++p) {
                        const auto N_value = [&](size_t population) {
                            const auto& values = cfg.grid_N[population];
                            return values[values.size() == 1 ? 0 : n];
                        };
                        const double log_base = base(t, n, p);
                        const double log_ghost = ghost(t, n, p);
                        const double log_ratio = log_base == -std::numeric_limits<double>::infinity()
                            && log_ghost == -std::numeric_limits<double>::infinity()
                            ? std::numeric_limits<double>::quiet_NaN() : log_base - log_ghost;
                        out << tree_index << '\t' << tree.interval.left << '\t' << tree.interval.right
                            << '\t' << t << '\t' << n << '\t' << p
                            << '\t' << cfg.grid_time[t] << '\t' << N_value(0) << '\t'
                            << N_value(1) << '\t' << N_value(2) << '\t' << cfg.grid_prob[p]
                            << '\t' << log_base << '\t' << log_ghost << '\t' << log_ratio << '\n';
                    }
            ++tree_index;
            continue;
        }
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

std::unique_ptr<treellh::TreeEnsemble> prepare_tree_ensemble(
    const Config& cfg, const std::string& input, tsk_treeseq_t& ts)
{
    using Clock = std::chrono::steady_clock;
    auto start = Clock::now();
    const bool persistent = cfg.mode != "cached" || !cfg.tree_cache.empty();
    const auto signature = persistent ? tree_cache_signature(cfg, input) : std::string{};
    const auto path = !persistent ? std::string{} : cfg.tree_cache.empty()
        ? input + ".tree-cache." + signature + ".bin" : cfg.tree_cache;
    auto ensemble = std::make_unique<treellh::TreeEnsemble>(ts, cfg.parameters,
        cfg.sample_population, cfg.admixed_index, cfg.outgroup_index, true,
        cfg.tree_num_specified ? static_cast<size_t>(cfg.tree_num) : 0);
    ProgressPrinter cache_progress;
    const treellh::CacheProgress report = [&](std::string_view stage, size_t done, size_t total) {
        cache_progress.print("Cache " + std::string(stage), done, total, done == total);
    };
    const bool saved = persistent && std::filesystem::exists(path);
    if (saved && !cfg.estimate_only) {
        const double max_bytes = cfg.cache_limit_mb * 1024.0 * 1024.0;
        ensemble->load(path, signature,
            max_bytes < static_cast<double>(std::numeric_limits<size_t>::max())
                ? static_cast<size_t>(max_bytes) : 0, report);
    }
    else ensemble->prepare(report);
    const auto estimate = ensemble->memory_estimate();
    std::cerr << "Cache estimate: trees=" << estimate.trees << " scenarios=" << estimate.scenarios
        << " records=" << estimate.records << " statistics_bytes=" << estimate.statistics_bytes
        << " resident_bytes=" << estimate.resident_bytes << " ("
        << estimate.resident_bytes / (1024.0 * 1024.0) << " MiB)\n"
        << "Tree statistics cache: " << (persistent ? path : "in memory")
        << (cfg.estimate_only ? " (estimate only)"
            : persistent ? (saved ? " (loaded)" : " (new)") : "") << '\n'
        << "Preparation/load seconds: " << std::chrono::duration<double>(Clock::now() - start).count() << '\n';
    std::cerr.flush();
    if (cfg.cache_limit_mb > 0 && estimate.resident_bytes / (1024.0 * 1024.0) > cfg.cache_limit_mb)
        throw std::runtime_error("estimated cache exceeds --cache-limit-mb");
    if (!cfg.estimate_only && !saved) {
        start = Clock::now();
        ensemble->build(report);
        std::cerr << "Cache build seconds: " << std::chrono::duration<double>(Clock::now() - start).count() << '\n';
        if (persistent) {
            start = Clock::now();
            ensemble->save(path, signature, report);
            std::cerr << "Cache save seconds: " << std::chrono::duration<double>(Clock::now() - start).count() << '\n';
        }
    }
    if (!cfg.estimate_only)
        std::cerr << "Cache ready: estimated resident " << estimate.resident_bytes / (1024.0 * 1024.0)
            << " MiB; statistics " << estimate.statistics_bytes / (1024.0 * 1024.0)
            << " MiB; peak process RSS " << peak_rss_mib() << " MiB\n";
    return ensemble;
}

void compute_cached_file(const Config& cfg, const std::string& input, const std::string& output)
{
    using Clock = std::chrono::steady_clock;
    tsk_treeseq_t ts = load_ts(input);
    struct Cleanup { tsk_treeseq_t* ts; ~Cleanup() { tsk_treeseq_free(ts); } } cleanup{&ts};
    auto ensemble = prepare_tree_ensemble(cfg, input, ts);
    if (cfg.estimate_only) return;

    std::ofstream file;
    if (!output.empty()) {
        file.open(output);
        if (!file) throw std::runtime_error("cannot open output file: " + output);
    }
    std::ostream& out = output.empty() ? std::cout : file;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << (cfg.per_tree ? "query_index\ttree_index\tleft\tright\tlog_llh\n" : "query_index\tlog_llh\n");
    out.flush();
    size_t query_index = 0;
    auto evaluate = [&](const std::vector<double>& parameters) {
        const auto query_start = Clock::now();
        if (cfg.per_tree) {
            const auto values = ensemble->evaluate_trees(parameters);
            const auto elapsed = std::chrono::duration<double>(Clock::now() - query_start).count();
            for (size_t i = 0; i < values.size(); ++i)
                out << query_index << '\t' << i << '\t' << ensemble->trees()[i].left << '\t'
                    << ensemble->trees()[i].right << '\t' << values[i] << '\n';
            std::cerr << "Query " << query_index << " seconds: " << elapsed << '\n';
        } else {
            const double value = ensemble->evaluate(parameters);
            const auto elapsed = std::chrono::duration<double>(Clock::now() - query_start).count();
            out << query_index << '\t' << value << '\n';
            std::cerr << "Query " << query_index << " seconds: " << elapsed << '\n';
        }
        out.flush();
        ++query_index;
    };
    if (cfg.queries.empty()) {
        evaluate(cfg.parameters);
    } else {
        std::ifstream queries_file;
        if (cfg.queries != "-") {
            queries_file.open(cfg.queries);
            if (!queries_file) throw std::runtime_error("cannot open query file: " + cfg.queries);
        }
        std::istream& queries = cfg.queries == "-" ? std::cin : queries_file;
        std::string line;
        size_t line_number = 0;
        while (std::getline(queries, line)) {
            ++line_number;
            line = line.substr(0, line.find('#'));
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream row(line);
            std::vector<double> parameters;
            double x;
            while (row >> x) parameters.push_back(x);
            if (!row.eof() || (!parameters.empty() && parameters.size() != 9))
                throw std::runtime_error("expected nine parameters on query line " + std::to_string(line_number));
            if (!parameters.empty()) evaluate(parameters);
        }
        if (queries.bad()) throw std::runtime_error("failed to read queries");
    }
}

void compute_fit_file(const Config& cfg, const std::string& input, const std::string& output)
{
    tsk_treeseq_t ts = load_ts(input);
    struct Cleanup { tsk_treeseq_t* ts; ~Cleanup() { tsk_treeseq_free(ts); } } cleanup{&ts};
    auto ensemble = prepare_tree_ensemble(cfg, input, ts);
    if (cfg.estimate_only) return;
    fit::Options options{};
    std::copy(cfg.parameters.begin(), cfg.parameters.end(), options.parameters.begin());
    options.max_evaluations = cfg.fit_maxeval;
    options.starts = cfg.fit_starts;
    for (size_t i = 0; i < 4; ++i)
        options.axes[i] = {cfg.fit_axes[i][0], cfg.fit_axes[i][1], static_cast<size_t>(cfg.fit_axes[i][2])};
    ProgressPrinter fit_progress;
    options.progress = [&](const fit::ProgressEvent& event) {
        const bool grid = event.phase == fit::ProgressEvent::Phase::grid;
        const std::string label = grid ? "Grid" : "Optimizer start "
            + std::to_string(event.start) + "/" + std::to_string(event.starts);
        fit_progress.print(label, event.done, event.total, event.complete, event.best_log_likelihood);
        if (grid && event.complete) {
            std::cerr << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "Grid optimum: log_llh=" << event.best_log_likelihood
                << " migration_time=" << event.best_parameters[0]
                << " migration_prob=" << event.best_parameters[3]
                << " Nbase=" << event.best_parameters[6]
                << " Nghost=" << event.best_parameters[7] << '\n';
        }
    };
    std::ofstream file;
    if (!output.empty()) {
        file.open(output);
        if (!file) throw std::runtime_error("cannot open output file: " + output);
    }
    std::ostream& out = output.empty() ? std::cout : file;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    const auto best = fit::search(*ensemble, options, cfg.mode.c_str(), out);
    std::cerr << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "Best log likelihood: " << best.log_likelihood << "; evaluations: " << best.evaluations
        << "; params: ";
    for (size_t i = 0; i < best.parameters.size(); ++i)
        std::cerr << (i ? "," : "") << best.parameters[i];
    std::cerr << '\n';
}

void run_single(const Config &cfg, const std::string &input, const std::string &output)
{
    if (cfg.mode == "cached") {
        compute_cached_file(cfg, input, output);
    } else if (cfg.mode == "fit-grid" || cfg.mode == "fit-opt" || cfg.mode == "fit-hybrid") {
        compute_fit_file(cfg, input, output);
    } else if (cfg.mode == "llh") {
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
