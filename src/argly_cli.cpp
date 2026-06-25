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

// Макрос проверки ошибок tskit
#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));                 \
    }

namespace {

// Структура конфигурации параметров командной строки
struct Config {
    std::string mode;           // Режим работы: "llh" (расчет в точке) или "grid" (на сетке параметров)
    std::string input;          // Путь к одиночному файлу дерева (*.trees или *.arg)
    std::string path;           // Путь к директории (для пакетного режима)
    std::string prefix = "arg";  // Префикс файлов симуляций
    std::string suffix = ".arg"; // Суффикс (расширение) файлов симуляций
    std::string output;         // Путь к выходному файлу результатов
    
    // Вектор из 9 параметров модели по умолчанию
    std::vector<double> parameters = {
        500, 1000, 2000,
        0.2, 0.5, 0.5,
        10000, 10000, 10000
    };
    
    std::vector<int> sample_population;        // Популяция каждого образца генома (0 - AFR, 2 - SAM)
    std::vector<double> grid_time = {500};      // Сетка времени миграции
    std::vector<double> grid_prob = {0.2};      // Сетка вероятности интрогрессии
    std::vector<std::vector<double>> grid_N;   // Сетка эффективных размеров Ne по популяциям
    int admixed_index = 0;                     // Индекс AFR в sample_population
    int outgroup_index = 1;                    // Индекс SAM в sample_population
    size_t sim_count = 0;                      // Количество симуляций в пакете
    long tree_num = 10000;                     // Число деревьев для прореживания (интервал сэмплирования)
    bool binary_output = true;                 // Флаг вывода результатов в бинарный (.bin), а не в текстовый формат
};

// Вывод справки по использованию CLI-интерфейса в поток std::cerr
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
        << "  --grid-N-base, --grid-N-ghost and --grid-N-outgroup add an N grid axis\n"
        << "  --tree-num K samples approximately K trees along the sequence\n"
        << "  --text-output writes grid as whitespace text instead of raw doubles\n";
}

// Разделение строки delim-символом (например, запятой)
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

// Преобразование строки с запятыми в вектор вещественных чисел double
std::vector<double> parse_doubles(const std::string &value)
{
    std::vector<double> result;
    for (const auto &item : split(value, ',')) {
        result.push_back(std::stod(item));
    }
    return result;
}

// Преобразование строки с запятыми в вектор целых чисел int
std::vector<int> parse_ints(const std::string &value)
{
    std::vector<int> result;
    for (const auto &item : split(value, ',')) {
        result.push_back(std::stoi(item));
    }
    return result;
}

// Генерация сетки параметров (linspace)
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

// Парсинг спецификации сетки. Принимает формат "start,end,count" для автоматического
// вызова linspace, либо обычный список значений через запятую (например: "500,600,700").
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

// Проверка наличия и извлечение значения следующего аргумента CLI
std::string require_value(int &i, int argc, char **argv)
{
    if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + argv[i]);
    }
    return argv[++i];
}

// Главный парсер аргументов CLI, возвращающий структуру Config
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

    // Проверка ограничений и валидности переданных аргументов
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
    
    // Если сетка эффективных размеров Ne не задана, инициализируем её значениями параметров по умолчанию
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

// Визуализация прогресса в консоли
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

// Формирует имя входного файла для i-й симуляции при пакетной обработке
std::string batch_input_name(const Config &cfg, size_t sim_i)
{
    std::filesystem::path path(cfg.path);
    path /= cfg.prefix + std::to_string(sim_i) + cfg.suffix;
    return path.string();
}

// Формирует имя выходного файла по умолчанию
std::string default_output_name(const std::string &input, const std::string &mode, bool binary)
{
    std::filesystem::path path(input);
    std::string ext = mode == "grid" ? (binary ? ".llh.bin" : ".llh.txt") : ".llh.txt";
    return path.replace_extension(ext).string();
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

// Сохраняет тензор результатов в текстовый формат
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

// Сохраняет скалярное значение логарифма правдоподобия в текстовый файл
void save_scalar_text(const std::string &filename, double value)
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("cannot open output file: " + filename);
    }
    file << value << '\n';
}

// Загружает tree sequence (последовательность деревьев) из файла при помощи tskit
tsk_treeseq_t load_ts(const std::string &filename)
{
    tsk_treeseq_t ts;
    int ret = tsk_treeseq_load(&ts, filename.c_str(), 0);
    if (ret != 0) {
        throw std::runtime_error("load error for " + filename + ": " + tsk_strerror(ret));
    }
    return ts;
}

// Расчет логарифма правдоподобия для одиночного файла при фиксированных параметрах (одна точка)
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
        result += computer.compute_llh(tree); // Суммируем логарифмы по всем локальным деревьям
    }
    check_tsk_error(ret);

    tsk_tree_free(&tree);
    tsk_treeseq_free(&ts);
    return result;
}

// Расчет тензора логарифмов правдоподобия для одиночного файла на сетке параметров
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
        // Выполняем сэмплирование/прореживание локальных деревьев
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

// Запуск расчета для одиночного файла (llh или grid) с выводом результатов
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
        
        // Режим одиночного запуска по файлу --input
        if (!cfg.input.empty()) {
            std::string output = cfg.output;
            if (output.empty() && cfg.mode == "grid") {
                output = default_output_name(cfg.input, cfg.mode, cfg.binary_output);
            }
            run_single(cfg, cfg.input, output);
            return EXIT_SUCCESS;
        }

        // Пакетный режим: итерируется по sim_count симуляциям
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
