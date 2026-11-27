#include "epoch_fit.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>
#ifdef __APPLE__
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

namespace {
std::vector<double> numbers(const std::string& s) {
    std::vector<double> r; std::stringstream in(s);std::string item;
    while(std::getline(in,item,',')) {size_t n=0;double x=std::stod(item,&n);if(n!=item.size() || !std::isfinite(x))throw std::invalid_argument("invalid numeric list");r.push_back(x);}
    return r;
}
size_t count(const std::string& s) {
    if(s.empty() || !std::all_of(s.begin(),s.end(),[](char c){return c>='0'&&c<='9';})) throw std::invalid_argument("expected a positive integer");
    size_t n=std::stoull(s);if(!n)throw std::invalid_argument("expected a positive integer");return n;
}
std::uint64_t seed(const std::string& s) {
    if(s.empty() || !std::all_of(s.begin(),s.end(),[](char c){return c>='0'&&c<='9';}))
        throw std::invalid_argument("expected an unsigned restart seed");
    return std::stoull(s);
}
double rss() {rusage u{};getrusage(RUSAGE_SELF,&u);
#ifdef __APPLE__
    return u.ru_maxrss/(1024.*1024.);
#else
    return u.ru_maxrss/1024.;
#endif
}
double current_rss() {
#ifdef __APPLE__
    mach_task_basic_info_data_t info{};
    mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,
        reinterpret_cast<task_info_t>(&info),&count)!=KERN_SUCCESS)
        throw std::runtime_error("cannot read current RSS");
    return info.resident_size/(1024.*1024.);
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    size_t total_pages=0,resident_pages=0;
    if(!(statm>>total_pages>>resident_pages))
        throw std::runtime_error("cannot read current RSS");
    const long page_size=sysconf(_SC_PAGESIZE);
    if(page_size<=0)throw std::runtime_error("cannot read page size");
    return static_cast<double>(resident_pages)*page_size/(1024.*1024.);
#else
    return rss();
#endif
}
void help() {
    std::cout << "Piecewise Nbase(t) and Nghost(t) ECME fit.\n"
        "argly-epoch-fit --input FILE --sample-population CSV --params CSV [options]\n"
        "argly-epoch-fit --prefix PREFIX --sim-count N --sample-population CSV --params CSV [options]\n"
        "  Batch files: PREFIX0.trees ... PREFIX(N-1).trees; --path DIR and --suffix are optional\n"
        "--params mt,split2,split1,mp,split2p,split1p,Nbase,Nghost,Noutgroup\n"
        "--base-epochs K --ghost-epochs K: separate grids (default: 1 each)\n"
        "--epoch-balance count|expected: adaptive boundaries by raw counts (default) or posterior expectations\n"
        "--base-boundaries CSV --ghost-boundaries CSV: explicit boundaries, including 0 and split1\n"
        "--epochs K: legacy shared grid of K equal time intervals from 0 to split1\n"
        "--epoch-boundaries CSV: K+1 increasing shared times; requires --epochs\n"
        "--ghost-sizes CSV: starting Nghost per epoch, youngest first (default: params Nghost)\n"
        "--base-sizes CSV: starting Nbase per base epoch, youngest first (default: params Nbase)\n"
        "--base-right-sizes CSV --ghost-right-sizes CSV: second N in linear/split modes\n"
        "--base-cuts CSV --ghost-cuts CSV: starting interior cut times for the selected cut mode\n"
        "--fit-base: accepted for compatibility; base fitting is now the default\n"
        "--fixed-base: keep supplied Nbase(t) fixed during fitting\n"
        "--fixed-ghost: keep supplied Nghost(t) fixed during fitting (default: estimate it)\n"
        "  Above split1, Nbase remains fixed at params Nbase.\n"
        "--time-bounds LO,HI --prob-bounds LO,HI --N-bounds LO,HI\n"
        "  defaults: [0,split2) with outgroup, [0,split1) without; [0.000001,0.999999], [100,10000000]\n"
        "--iterations K (30) --restarts K (10 random starts per migration search)\n"
        "--threads K (automatic by default; 1 runs tree loops serially)\n"
        "--linear-rates: fit 1/(2N(t)) linearly within each epoch (separate endpoint sizes)\n"
        "--movable-cuts: fit two constant sizes and one movable cut inside each epoch\n"
        "--base-movable-cuts: movable cuts in base epochs; constant N per ghost epoch\n"
        "--maxeval K (100 per restart) --restart-seed K (5489; reproducible starts)\n"
        "--tolerance X (1e-5 absolute log likelihood)\n"
        "--tree-num K: one file uses genomic probes; batch divides K tree indices evenly among files\n"
        "--admixed-index I (0) --outgroup-index I (1)\n"
        "--cache-limit-mb X (1024 MiB, 0 unlimited); checked before allocating each tree\n"
        "--evaluate-only: evaluate supplied parameters without fitting\n"
        "--output FILE: iteration TSV; otherwise stdout. Progress goes to stderr.\n"
        "--profile-output FILE: final epochs and information status (default: OUTPUT.epochs.tsv)\n"
        "--kernel-output FILE: posterior-weighted Gaussian N(t) after fitting\n"
        "--kernel-bandwidth X (default: split1/50) --kernel-points K (201)\n";
}
template<class Work> void with_tree_sequence(const std::string& path,Work work) {
    tsk_treeseq_t ts{};
    const int ret=tsk_treeseq_load(&ts,path.c_str(),0);
    if(ret<0)throw std::runtime_error("cannot load "+path+": "+tsk_strerror(ret));
    struct Cleanup {tsk_treeseq_t& ts;~Cleanup(){tsk_treeseq_free(&ts);}} cleanup{ts};
    work(ts);
}
}
int main(int argc,char** argv) {
    try {
        std::string input,path,prefix="arg",suffix=".trees",output,profile_output,kernel_output;
        size_t sim_count=0;
        std::vector<double> p,initial_sizes,initial_base,epoch_boundaries;
        std::vector<double> initial_base_right,initial_ghost_right;
        std::vector<double> initial_base_cuts,initial_ghost_cuts;
        std::vector<double> explicit_base_boundaries,explicit_ghost_boundaries;
        std::vector<int> populations;size_t epochs=1,base_epochs=1,shared_epochs=0,tree_limit=0;
        int admixed=0,outgroup=1;double memory_limit=1024;size_t threads=0;
        double kernel_bandwidth=0;size_t kernel_points=201;
        bool kernel_bandwidth_specified=false,kernel_points_specified=false;
        bool evaluate=false,has_time_bounds=false,fit_base=true,fit_base_requested=false,fixed_base_requested=false;
        bool separate_epochs_requested=false,epoch_boundaries_requested=false;
        bool base_epochs_requested=false,ghost_epochs_requested=false,balance_requested=false;
        bool base_boundaries_requested=false,ghost_boundaries_requested=false;
        std::string balance="count";
        epochfit::FitOptions options;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            auto value=[&](){if(++i==argc)throw std::invalid_argument("missing value for "+arg);return std::string(argv[i]);};
            auto bounds=[&](){auto v=numbers(value());if(v.size()!=2)throw std::invalid_argument("bounds need LO,HI");return std::array<double,2>{v[0],v[1]};};
            if(arg=="--help" || arg=="-h") {help();return 0;}
            else if(arg=="--input")input=value();
            else if(arg=="--path")path=value();
            else if(arg=="--prefix")prefix=value();
            else if(arg=="--suffix")suffix=value();
            else if(arg=="--sim-count")sim_count=count(value());
            else if(arg=="--output")output=value();
            else if(arg=="--profile-output")profile_output=value();
            else if(arg=="--kernel-output")kernel_output=value();
            else if(arg=="--kernel-bandwidth") {
                const auto values=numbers(value());
                if(values.size()!=1)throw std::invalid_argument("--kernel-bandwidth needs one number");
                kernel_bandwidth=values[0];kernel_bandwidth_specified=true;
            }
            else if(arg=="--kernel-points") {kernel_points=count(value());kernel_points_specified=true;}
            else if(arg=="--params")p=numbers(value());
            else if(arg=="--ghost-sizes")initial_sizes=numbers(value());
            else if(arg=="--base-sizes")initial_base=numbers(value());
            else if(arg=="--ghost-right-sizes")initial_ghost_right=numbers(value());
            else if(arg=="--base-right-sizes")initial_base_right=numbers(value());
            else if(arg=="--ghost-cuts")initial_ghost_cuts=numbers(value());
            else if(arg=="--base-cuts")initial_base_cuts=numbers(value());
            else if(arg=="--fit-base") {fit_base=true;fit_base_requested=true;}
            else if(arg=="--fixed-base") {fit_base=false;fixed_base_requested=true;}
            else if(arg=="--fixed-ghost")options.fit_ghost=false;
            else if(arg=="--sample-population") {
                for(double x:numbers(value())) {if(x!=std::floor(x)||x<0||x>2147483647)throw std::invalid_argument("invalid population ID");populations.push_back(int(x));}
            }
            else if(arg=="--epochs")shared_epochs=count(value());
            else if(arg=="--base-epochs") {base_epochs=count(value());separate_epochs_requested=base_epochs_requested=true;}
            else if(arg=="--ghost-epochs") {epochs=count(value());separate_epochs_requested=ghost_epochs_requested=true;}
            else if(arg=="--epoch-balance") {balance=value();balance_requested=true;
                if(balance!="count" && balance!="expected")throw std::invalid_argument("--epoch-balance needs count or expected");}
            else if(arg=="--base-boundaries") {explicit_base_boundaries=numbers(value());separate_epochs_requested=base_boundaries_requested=true;}
            else if(arg=="--ghost-boundaries") {explicit_ghost_boundaries=numbers(value());separate_epochs_requested=ghost_boundaries_requested=true;}
            else if(arg=="--epoch-boundaries") {epoch_boundaries=numbers(value());epoch_boundaries_requested=true;}
            else if(arg=="--tree-num")tree_limit=count(value());
            else if(arg=="--iterations")options.iterations=count(value());
            else if(arg=="--restarts")options.restarts=count(value());
            else if(arg=="--threads")threads=count(value());
            else if(arg=="--linear-rates")options.linear_rates=true;
            else if(arg=="--movable-cuts")options.movable_cuts=true;
            else if(arg=="--base-movable-cuts")options.base_movable_cuts=true;
            else if(arg=="--maxeval")options.max_evaluations=count(value());
            else if(arg=="--restart-seed")options.restart_seed=seed(value());
            else if(arg=="--tolerance")options.tolerance=std::stod(value());
            else if(arg=="--admixed-index")admixed=std::stoi(value());
            else if(arg=="--outgroup-index")outgroup=std::stoi(value());
            else if(arg=="--cache-limit-mb")memory_limit=std::stod(value());
            else if(arg=="--time-bounds") {options.time_bounds=bounds();has_time_bounds=true;}
            else if(arg=="--prob-bounds")options.probability_bounds=bounds();
            else if(arg=="--N-bounds") {auto b=bounds();options.min_size=b[0];options.max_size=b[1];}
            else if(arg=="--evaluate-only")evaluate=true;
            else throw std::invalid_argument("unknown option: "+arg);
        }
        if(p.size()!=9 || populations.empty())
            throw std::invalid_argument("require --params (9 values) and --sample-population");
        if((!input.empty())==(sim_count!=0) || (!input.empty() && !path.empty()))
            throw std::invalid_argument("provide --input FILE or --prefix PREFIX --sim-count N");
        std::vector<std::string> inputs;
        if(!input.empty())inputs.push_back(input);
        else {
            if(prefix.empty() || suffix.empty())throw std::invalid_argument("batch prefix and suffix must be nonempty");
            inputs.reserve(sim_count);
            for(size_t i=0;i<sim_count;++i)
                inputs.push_back((std::filesystem::path(path.empty()?".":path)
                    /(prefix+std::to_string(i)+suffix)).string());
        }
        for(const auto& source:inputs)
            if(!std::filesystem::is_regular_file(source))
                throw std::invalid_argument("input tree sequence is missing: "+source);
        if(tree_limit && tree_limit<inputs.size())
            throw std::invalid_argument("--tree-num must be at least the number of input files");
        auto file_tree_limit=[&](size_t i) {
            return tree_limit?tree_limit/inputs.size()+(i<tree_limit%inputs.size()):size_t(0);
        };
        const bool even_tree_indices=sim_count!=0;
        if(shared_epochs && separate_epochs_requested)
            throw std::invalid_argument("--epochs cannot be combined with separate population grids");
        if(shared_epochs && balance_requested)
            throw std::invalid_argument("--epoch-balance applies only to separate population grids");
        if(fit_base_requested && fixed_base_requested)
            throw std::invalid_argument("--fit-base conflicts with --fixed-base");
        if((options.linear_rates && (options.movable_cuts || options.base_movable_cuts))
            || (options.movable_cuts && options.base_movable_cuts))
            throw std::invalid_argument("choose one of --linear-rates, --movable-cuts, --base-movable-cuts");
        if(shared_epochs)epochs=base_epochs=shared_epochs;
        if(base_boundaries_requested) {
            if(explicit_base_boundaries.size()<2)throw std::invalid_argument("--base-boundaries needs at least two values");
            if(base_epochs_requested && explicit_base_boundaries.size()!=base_epochs+1)
                throw std::invalid_argument("--base-boundaries needs base epochs + 1 values");
            base_epochs=explicit_base_boundaries.size()-1;
        }
        if(ghost_boundaries_requested) {
            if(explicit_ghost_boundaries.size()<2)throw std::invalid_argument("--ghost-boundaries needs at least two values");
            if(ghost_epochs_requested && explicit_ghost_boundaries.size()!=epochs+1)
                throw std::invalid_argument("--ghost-boundaries needs ghost epochs + 1 values");
            epochs=explicit_ghost_boundaries.size()-1;
        }
        if(!base_epochs || !epochs)throw std::invalid_argument("boundaries need at least two values");
        if(epochs>100000 || base_epochs>100000)throw std::invalid_argument("epoch count must be at most 100000");
        if(epoch_boundaries_requested && !shared_epochs)
            throw std::invalid_argument("--epoch-boundaries requires --epochs shared-grid mode");
        if(epoch_boundaries_requested && epoch_boundaries.size()!=shared_epochs+1)
            throw std::invalid_argument("--epoch-boundaries needs K+1 values");
        if(profile_output.empty() && !output.empty())profile_output=output+".epochs.tsv";
        auto same_path=[](const std::string& a,const std::string& b) {
            if(a.empty() || b.empty())return false;
            if(std::filesystem::absolute(a).lexically_normal()
                ==std::filesystem::absolute(b).lexically_normal())return true;
            std::error_code error;
            return std::filesystem::equivalent(a,b,error) && !error;
        };
        for(const auto& source:inputs)
            if(same_path(output,source) || same_path(profile_output,source)
                || same_path(kernel_output,source))
                throw std::invalid_argument("output path must differ from every input file");
        if(same_path(output,profile_output))
            throw std::invalid_argument("iteration output and profile output must differ");
        if(!kernel_output.empty() && (same_path(kernel_output,output) || same_path(kernel_output,profile_output)))
            throw std::invalid_argument("kernel output must differ from input and other outputs");
        if(kernel_output.empty() && (kernel_bandwidth_specified || kernel_points_specified))
            throw std::invalid_argument("kernel settings require --kernel-output");
        if(!kernel_bandwidth_specified)kernel_bandwidth=p[2]/50.;
        if(!std::isfinite(kernel_bandwidth) || !(kernel_bandwidth>0) || kernel_points<2 || kernel_points>1000000)
            throw std::invalid_argument("invalid kernel bandwidth or point count");
        if(!std::isfinite(memory_limit)||memory_limit<0)throw std::invalid_argument("invalid memory limit");
        if(initial_sizes.empty())initial_sizes.assign(epochs,p[7]);
        if(initial_sizes.size()!=epochs)throw std::invalid_argument("--ghost-sizes needs one value per epoch");
        if(initial_base.empty())initial_base.assign(base_epochs,p[6]);
        if(initial_base.size()!=base_epochs)throw std::invalid_argument("--base-sizes needs one value per base epoch");
        if(!options.linear_rates && !options.movable_cuts && !options.base_movable_cuts
            && (!initial_base_right.empty() || !initial_ghost_right.empty()))
            throw std::invalid_argument("right-end sizes require --linear-rates or --movable-cuts");
        if(!options.movable_cuts && !options.base_movable_cuts
            && (!initial_base_cuts.empty() || !initial_ghost_cuts.empty()))
            throw std::invalid_argument("interior cuts require --movable-cuts");
        if(options.base_movable_cuts && (!initial_ghost_cuts.empty() || !initial_ghost_right.empty()))
            throw std::invalid_argument("--base-movable-cuts has no ghost cut or right size");
        if(initial_base_right.empty())initial_base_right=initial_base;
        if(initial_ghost_right.empty() && !options.base_movable_cuts)initial_ghost_right=initial_sizes;
        if(initial_base_right.size()!=base_epochs
            || (!options.base_movable_cuts && initial_ghost_right.size()!=epochs))
            throw std::invalid_argument("right-end sizes need one value per epoch");
        options.fit_base=fit_base;
        std::ofstream file;
        if(!output.empty()) {file.open(output);if(!file)throw std::runtime_error("cannot open output file");}
        std::ostream& out=output.empty()?std::cout:file;
        out<<std::setprecision(17);std::cerr<<std::setprecision(12);
        std::cerr<<"Input tree sequences ("<<inputs.size()<<"):\n";
        for(size_t i=0;i<inputs.size();++i)
            std::cerr<<"  ["<<i+1<<'/'<<inputs.size()<<"] "<<inputs[i]<<'\n';
        if(!has_time_bounds) {
            with_tree_sequence(inputs.front(),[&](const tsk_treeseq_t& source) {
                bool has_outgroup=false;
                const auto* samples=tsk_treeseq_get_samples(&source);
                for(size_t i=0;i<tsk_treeseq_get_num_samples(&source);++i)
                    if(size_t(samples[i])<populations.size() && populations[samples[i]]==outgroup)
                        {has_outgroup=true;break;}
                options.time_bounds={0,std::nextafter(p[has_outgroup?1:2],0.)};
            });
        }
        using Clock=std::chrono::steady_clock;
        auto start=Clock::now();
        std::vector<double> ghost_boundaries=shared_epochs?epoch_boundaries:explicit_ghost_boundaries;
        std::vector<double> base_boundaries=explicit_base_boundaries;
        const bool choose_boundaries=!shared_epochs
            && ((epochs>1 && ghost_boundaries.empty()) || (base_epochs>1 && base_boundaries.empty()));
        if(choose_boundaries) {
            if(balance=="expected") {
                std::optional<epochfit::Ensemble> pilot;
                for(size_t i=0;i<inputs.size();++i)
                    with_tree_sequence(inputs[i],[&](const tsk_treeseq_t& next) {
                        if(i==0)pilot.emplace(next,p,populations,admixed,outgroup,1,
                            file_tree_limit(i),memory_limit,
                            std::function<void(size_t,size_t,size_t)>{},
                            std::vector<double>{},0,std::vector<double>{},
                            std::array<double,2>{p[0],p[0]},even_tree_indices);
                        else pilot->append(next,p,populations,admixed,outgroup,
                            file_tree_limit(i),memory_limit,{},even_tree_indices);
                    });
                if(inputs.size()>1)pilot->set_sizes(std::vector<double>(1,p[7]));
                const auto masses=pilot->coalescence_masses(p[0],p[3]);
                if(base_boundaries.empty() && base_epochs>1)
                    base_boundaries=epochfit::balanced_epoch_boundaries(masses[0],base_epochs,p[2]);
                if(ghost_boundaries.empty() && epochs>1)
                    ghost_boundaries=epochfit::balanced_epoch_boundaries(masses[1],epochs,p[2]);
            } else {
                std::map<double,double> combined;
                size_t scanned=0;
                for(size_t i=0;i<inputs.size();++i) {
                    auto collect=[&](const tsk_treeseq_t& source) {
                        scanned+=source.num_trees;
                        for(const auto& [time,mass]:epochfit::raw_coalescence_counts(source,p[2]))
                            combined[time]+=mass;
                    };
                    with_tree_sequence(inputs[i],collect);
                }
                const std::vector<std::pair<double,double>> raw_counts(combined.begin(),combined.end());
                std::cerr<<"Counted coalescences in all "<<scanned<<" input trees at "
                    <<raw_counts.size()<<" distinct times\n";
                if(base_boundaries.empty() && base_epochs>1)
                    base_boundaries=epochfit::balanced_epoch_boundaries(raw_counts,base_epochs,p[2]);
                if(ghost_boundaries.empty() && epochs>1)
                    ghost_boundaries=epochfit::balanced_epoch_boundaries(raw_counts,epochs,p[2]);
            }
        }
        auto progress=[&](size_t file_index) {
            return [&,file_index,last=-1](size_t done,size_t total,size_t bytes) mutable {
                const int bucket=total?int(10.*done/total)*10:100;
                if(bucket>last) {
                    last=bucket;
                    std::cerr<<"Epoch cache file "<<file_index+1<<'/'<<inputs.size()<<": "<<bucket
                        <<"% (scanned "<<done<<'/'<<total<<" trees), "
                        <<bytes/(1024.*1024.)<<" MiB total\n";
                }
            };
        };
        std::optional<epochfit::Ensemble> pooled;
        std::array<std::vector<epochfit::RawEpochStatistics>,2> raw_information;
        size_t raw_tree_count=0;
        for(size_t i=0;i<inputs.size();++i)
            with_tree_sequence(inputs[i],[&](const tsk_treeseq_t& next) {
                if(i==0)pooled.emplace(next,p,populations,admixed,outgroup,epochs,
                    file_tree_limit(i),memory_limit,progress(i),ghost_boundaries,
                    base_epochs,base_boundaries,options.time_bounds,even_tree_indices);
                else pooled->append(next,p,populations,admixed,outgroup,
                    file_tree_limit(i),memory_limit,progress(i),even_tree_indices);
                const auto part=epochfit::raw_epoch_statistics(
                    next,pooled->base_boundaries(),pooled->boundaries());
                raw_tree_count+=next.num_trees;
                if(i==0)raw_information=part;
                else for(size_t q=0;q<2;++q)for(size_t e=0;e<raw_information[q].size();++e) {
                        raw_information[q][e].coalescences+=part[q][e].coalescences;
                        raw_information[q][e].lineage_time+=part[q][e].lineage_time;
                        raw_information[q][e].pair_time+=part[q][e].pair_time;
                    }
            });
        auto& ensemble=*pooled;
        std::cerr<<"Cache built: selected trees="<<ensemble.tree_count()
            <<" scenarios="<<ensemble.scenario_count()
            <<" resident arrays="<<ensemble.bytes()/(1024.*1024.)
            <<" MiB; current RSS="<<current_rss()<<" MiB; peak RSS="<<rss()<<" MiB\n";
        ensemble.set_threads(threads);
        if(options.movable_cuts || options.base_movable_cuts) {
            auto defaults=[](const std::vector<double>& bounds) {
                std::vector<double> cuts;
                for(size_t e=0;e+1<bounds.size();++e)
                    cuts.push_back(bounds[e]+(bounds[e+1]-bounds[e])*0.5);
                return cuts;
            };
            if(initial_base_cuts.empty())initial_base_cuts=defaults(ensemble.base_boundaries());
            if(options.movable_cuts && initial_ghost_cuts.empty())
                initial_ghost_cuts=defaults(ensemble.boundaries());
        }
        if(options.linear_rates)ensemble.set_linear_sizes(initial_base,initial_base_right,
            initial_sizes,initial_ghost_right);
        else if(options.movable_cuts)ensemble.set_split_sizes(initial_base,initial_base_right,
            initial_base_cuts,initial_sizes,initial_ghost_right,initial_ghost_cuts);
        else if(options.base_movable_cuts)ensemble.set_split_base_sizes(
            initial_base,initial_base_right,initial_base_cuts,initial_sizes);
        else ensemble.set_sizes(initial_base,initial_sizes);
        std::cerr<<"Cache ready: selected trees="<<ensemble.tree_count()
            <<" scenarios="<<ensemble.scenario_count()
            <<" resident arrays="<<ensemble.bytes()/(1024.*1024.)
            <<" MiB; current RSS="<<current_rss()<<" MiB; peak RSS="<<rss()<<" MiB; seconds="
            <<std::chrono::duration<double>(Clock::now()-start).count()<<'\n';
        std::cerr<<"Base boundaries (youngest first):";for(double b:ensemble.base_boundaries())std::cerr<<' '<<b;std::cerr<<'\n';
        std::cerr<<"Ghost boundaries (youngest first):";for(double b:ensemble.boundaries())std::cerr<<' '<<b;std::cerr<<'\n';
        const auto initial_expected=ensemble.expectation(p[0],p[3]);
        std::cerr<<"Initial epoch information (raw columns: all "<<raw_tree_count
            <<" input trees, uncoloured; expected C: "<<ensemble.tree_count()
            <<" fitted trees at initial parameters):\n"
            <<"population\tepoch\t(left,right]\traw_C\traw_mean_K\traw_lineage_time"
            <<"\traw_pair_time\texpected_C\tN0\tideal_var_N\tideal_sd_N\n";
        auto report_epochs=[&](const char* name,const std::vector<double>& boundaries,
            const std::vector<double>& sizes,const std::vector<double>& expected_counts,
            const std::vector<epochfit::RawEpochStatistics>& raw) {
            double total_raw=0,total_expected=0;
            double min_raw=std::numeric_limits<double>::infinity(),max_raw=0;
            for(size_t e=0;e<sizes.size();++e) {
                const double c=expected_counts[e];
                const double sd=c>0?sizes[e]/std::sqrt(c):std::numeric_limits<double>::infinity();
                const double mean_k=raw[e].lineage_time/
                    (double(raw_tree_count)*(boundaries[e+1]-boundaries[e]));
                std::cerr<<name<<'\t'<<e<<"\t("<<boundaries[e]<<','<<boundaries[e+1]
                    <<"]\t"<<raw[e].coalescences<<'\t'<<mean_k<<'\t'
                    <<raw[e].lineage_time<<'\t'<<raw[e].pair_time<<'\t'
                    <<c<<'\t'<<sizes[e]<<'\t'<<sd*sd<<'\t'<<sd<<'\n';
                total_raw+=raw[e].coalescences;total_expected+=c;
                min_raw=std::min(min_raw,raw[e].coalescences);
                max_raw=std::max(max_raw,raw[e].coalescences);
            }
            std::cerr<<name<<" mean per epoch: raw_C="<<total_raw/sizes.size()
                <<" (min="<<min_raw<<", max="<<max_raw<<"), expected_C="
                <<total_expected/sizes.size()<<'\n';
        };
        report_epochs("base",ensemble.base_boundaries(),initial_base,
            initial_expected.base_coalescences,raw_information[0]);
        report_epochs("ghost",ensemble.boundaries(),initial_sizes,
            initial_expected.coalescences,raw_information[1]);
        std::cerr<<"Raw mean K is the time- and tree-averaged active lineage count; "
            <<"raw pair time integrates K(K-1)/2 over each window. "
            <<"Raw columns are computed once before ECME and ignore population colouring.\n";
        std::cerr<<"Idealized variance uses N0^2/expected_C (SD=N0/sqrt(expected_C)); "
            <<"it is not a confidence interval and ignores latent-colouring and tree dependence.\n";
        out<<"iteration\tlog_llh\tT_adm\tp\timprovement";
        for(size_t e=0;e<base_epochs;++e) {
            out<<((options.linear_rates || options.movable_cuts || options.base_movable_cuts)
                ?"\tNbase_left_":"\tNbase_")<<e;
            if(options.linear_rates || options.movable_cuts || options.base_movable_cuts)
                out<<"\tNbase_right_"<<e;
            if(options.movable_cuts || options.base_movable_cuts)out<<"\tbase_cut_"<<e;
        }
        for(size_t e=0;e<epochs;++e) {
            out<<((options.linear_rates || options.movable_cuts)?"\tNghost_left_":"\tNghost_")<<e;
            if(options.linear_rates || options.movable_cuts)out<<"\tNghost_right_"<<e;
            if(options.movable_cuts)out<<"\tghost_cut_"<<e;
        }
        out<<'\n';
        options.iteration=[&](size_t it,double llh,const epochfit::Parameters& state,double gain){
            out<<it<<'\t'<<llh<<'\t'<<state.time<<'\t'<<state.probability<<'\t'<<gain;
            for(size_t e=0;e<state.base_sizes.size();++e) {
                out<<'\t'<<state.base_sizes[e];
                if(options.linear_rates || options.movable_cuts || options.base_movable_cuts)
                    out<<'\t'<<state.base_right_sizes[e];
                if(options.movable_cuts || options.base_movable_cuts)out<<'\t'<<state.base_cuts[e];
            }
            for(size_t e=0;e<state.ghost_sizes.size();++e) {
                out<<'\t'<<state.ghost_sizes[e];
                if(options.linear_rates || options.movable_cuts)out<<'\t'<<state.ghost_right_sizes[e];
                if(options.movable_cuts)out<<'\t'<<state.ghost_cuts[e];
            }
            out<<'\n';out.flush();
            std::cerr<<"ECME "<<it<<": log_llh="<<llh<<" improvement="<<gain<<" T_adm="<<state.time<<" p="<<state.probability;
            std::cerr<<" Nbase=";for(double n:state.base_sizes)std::cerr<<n<<',';
            if(options.linear_rates || options.movable_cuts || options.base_movable_cuts) {std::cerr<<" Nbase_right=";
                for(double n:state.base_right_sizes)std::cerr<<n<<',';}
            std::cerr<<" Nghost=";
            for(double n:state.ghost_sizes)std::cerr<<n<<',';
            if(options.linear_rates || options.movable_cuts) {std::cerr<<" Nghost_right=";
                for(double n:state.ghost_right_sizes)std::cerr<<n<<',';}
            std::cerr<<'\n';
        };
        options.optimizer=[&](size_t it,size_t calls,double best){
            const size_t budget=options.max_evaluations*options.restarts;
            if(calls==1 || calls%10==0 || calls==budget)
                std::cerr<<"Migration search "<<it<<": "<<calls<<'/'<<budget<<" evaluations across "
                    <<options.restarts<<" restarts, best="<<best<<'\n';
        };
        epochfit::Parameters state{p[0],p[3],initial_sizes,initial_base};
        if(options.linear_rates || options.movable_cuts || options.base_movable_cuts) {
            if(!options.base_movable_cuts)state.ghost_right_sizes=initial_ghost_right;
            state.base_right_sizes=initial_base_right;
        }
        if(options.movable_cuts || options.base_movable_cuts) {
            if(options.movable_cuts)state.ghost_cuts=initial_ghost_cuts;
            state.base_cuts=initial_base_cuts;
        }
        if(evaluate) options.iteration(0,ensemble.evaluate(state.time,state.probability),state,0);
        else {
            const auto result=epochfit::fit(ensemble,state,options);
            state=result.parameters;
            std::cerr<<"Finished: "<<result.status<<"; iterations="<<result.iterations<<"; log_llh="<<result.log_likelihood<<'\n';
        }
        if(!profile_output.empty()) {
            std::ofstream profile(profile_output);if(!profile)throw std::runtime_error("cannot open profile output");
            profile<<std::setprecision(17);
            const auto expected=ensemble.expectation(state.time,state.probability);
            auto ghost_status=[&](size_t e) {
                std::string status="active";
                if(ensemble.boundaries()[e+1]<=state.time)status="below_migration";
                else if(!options.fit_ghost)status="fixed";
                else if(expected.exposure[e]==0 && expected.coalescences[e]==0)status="no_information";
                else if(state.ghost_sizes[e]==options.max_size
                    || ((options.linear_rates || options.movable_cuts)
                        && state.ghost_right_sizes[e]==options.max_size))status="upper_bound";
                else if(state.ghost_sizes[e]==options.min_size
                    || ((options.linear_rates || options.movable_cuts)
                        && state.ghost_right_sizes[e]==options.min_size))status="lower_bound";
                return status;
            };
            auto base_status=[&](size_t e) {
                std::string status="active";
                if(!options.fit_base)status="fixed";
                else if(expected.base_exposure[e]==0 && expected.base_coalescences[e]==0)status="no_information";
                else if(state.base_sizes[e]==options.max_size
                    || ((options.linear_rates || options.movable_cuts || options.base_movable_cuts)
                        && state.base_right_sizes[e]==options.max_size))status="upper_bound";
                else if(state.base_sizes[e]==options.min_size
                    || ((options.linear_rates || options.movable_cuts || options.base_movable_cuts)
                        && state.base_right_sizes[e]==options.min_size))status="lower_bound";
                return status;
            };
            if(options.base_movable_cuts) {
                profile<<"population\tepoch\tleft_time\tright_time\tN_left\tN_right\tcut_time"
                    <<"\texpected_coalescences\texposure\tstatus\n";
                for(size_t e=0;e<base_epochs;++e)
                    profile<<"base\t"<<e<<'\t'<<ensemble.base_boundaries()[e]
                        <<'\t'<<ensemble.base_boundaries()[e+1]<<'\t'<<state.base_sizes[e]
                        <<'\t'<<state.base_right_sizes[e]<<'\t'<<state.base_cuts[e]
                        <<'\t'<<expected.base_coalescences[e]
                        <<'\t'<<expected.base_exposure[e]<<'\t'<<base_status(e)<<'\n';
                for(size_t e=0;e<epochs;++e)
                    profile<<"ghost\t"<<e<<'\t'<<ensemble.boundaries()[e]
                        <<'\t'<<ensemble.boundaries()[e+1]<<'\t'<<state.ghost_sizes[e]
                        <<"\t\t\t"<<expected.coalescences[e]
                        <<'\t'<<expected.exposure[e]<<'\t'<<ghost_status(e)<<'\n';
            } else if(options.linear_rates || options.movable_cuts) {
                profile<<"population\tepoch\tleft_time\tright_time\tN_left\tN_right";
                if(options.movable_cuts)profile<<"\tcut_time";
                profile
                    <<"\texpected_coalescences\texposure\tstatus\n";
                for(size_t e=0;e<base_epochs;++e) {
                    profile<<"base\t"<<e<<'\t'<<ensemble.base_boundaries()[e]
                        <<'\t'<<ensemble.base_boundaries()[e+1]<<'\t'<<state.base_sizes[e]
                        <<'\t'<<state.base_right_sizes[e];
                    if(options.movable_cuts)profile<<'\t'<<state.base_cuts[e];
                    profile<<'\t'<<expected.base_coalescences[e]
                        <<'\t'<<expected.base_exposure[e]<<'\t'<<base_status(e)<<'\n';
                }
                for(size_t e=0;e<epochs;++e) {
                    profile<<"ghost\t"<<e<<'\t'<<ensemble.boundaries()[e]
                        <<'\t'<<ensemble.boundaries()[e+1]<<'\t'<<state.ghost_sizes[e]
                        <<'\t'<<state.ghost_right_sizes[e];
                    if(options.movable_cuts)profile<<'\t'<<state.ghost_cuts[e];
                    profile<<'\t'<<expected.coalescences[e]
                        <<'\t'<<expected.exposure[e]<<'\t'<<ghost_status(e)<<'\n';
                }
            } else if(shared_epochs) {
                profile<<"epoch\tleft_time\tright_time\tNbase\tbase_expected_coalescences\tbase_exposure\tbase_status"
                    <<"\tNghost\texpected_coalescences\texposure\tstatus\n";
                for(size_t e=0;e<epochs;++e) {
                    profile<<e<<'\t'<<ensemble.boundaries()[e]<<'\t'<<ensemble.boundaries()[e+1]
                        <<'\t'<<state.base_sizes[e]<<'\t'<<expected.base_coalescences[e]
                        <<'\t'<<expected.base_exposure[e]<<'\t'<<base_status(e)
                        <<'\t'<<state.ghost_sizes[e]<<'\t'<<expected.coalescences[e]
                        <<'\t'<<expected.exposure[e]<<'\t'<<ghost_status(e)<<'\n';
                }
            } else {
                profile<<"population\tepoch\tleft_time\tright_time\tN\texpected_coalescences\texposure\tstatus\n";
                for(size_t e=0;e<base_epochs;++e)
                    profile<<"base\t"<<e<<'\t'<<ensemble.base_boundaries()[e]
                        <<'\t'<<ensemble.base_boundaries()[e+1]<<'\t'<<state.base_sizes[e]
                        <<'\t'<<expected.base_coalescences[e]<<'\t'<<expected.base_exposure[e]
                        <<'\t'<<base_status(e)<<'\n';
                for(size_t e=0;e<epochs;++e)
                    profile<<"ghost\t"<<e<<'\t'<<ensemble.boundaries()[e]
                        <<'\t'<<ensemble.boundaries()[e+1]<<'\t'<<state.ghost_sizes[e]
                        <<'\t'<<expected.coalescences[e]<<'\t'<<expected.exposure[e]
                        <<'\t'<<ghost_status(e)<<'\n';
            }
            if(!profile)throw std::runtime_error("failed writing profile output");
            std::cerr<<"Epoch profile: "<<profile_output<<'\n';
        }
        if(!kernel_output.empty()) {
            const auto curve=ensemble.kernel_profile(
                state.time,state.probability,kernel_bandwidth,kernel_points);
            std::ofstream kernel(kernel_output);
            if(!kernel)throw std::runtime_error("cannot open kernel output file");
            kernel<<std::setprecision(17)
                <<"population\ttime\tN_kernel\tkernel_C\tkernel_E\tbandwidth\tstatus\n";
            for(size_t q=0;q<2;++q)for(const auto& row:curve[q]) {
                const char* status=q==1 && row.time<=state.time?"below_migration"
                    :row.coalescences<=0?"no_events"
                    :row.exposure<=0?"no_exposure"
                    :row.coalescences<1?"low_events":"ok";
                kernel<<(q==0?"base":"ghost")<<'\t'<<row.time<<'\t';
                if(std::isfinite(row.size))kernel<<row.size;
                else kernel<<"NA";
                kernel<<'\t'<<row.coalescences<<'\t'<<row.exposure<<'\t'
                    <<kernel_bandwidth<<'\t'<<status<<'\n';
            }
            if(!kernel)throw std::runtime_error("failed writing kernel output");
            std::cerr<<"Kernel profile: "<<kernel_output<<'\n';
        }
        std::cerr<<"Total seconds="<<std::chrono::duration<double>(Clock::now()-start).count()<<"; peak RSS="<<rss()<<" MiB\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<"epoch-fit: "<<e.what()<<'\n';return 1;}
}
