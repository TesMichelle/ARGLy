#include "epoch_fit.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>

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
double rss() {rusage u{};getrusage(RUSAGE_SELF,&u);
#ifdef __APPLE__
    return u.ru_maxrss/(1024.*1024.);
#else
    return u.ru_maxrss/1024.;
#endif
}
void help() {
    std::cout << "Piecewise Nghost(t) and optional Nbase(t) ECME fit.\n"
        "argly-epoch-fit --input FILE --sample-population CSV --params CSV --epochs K [options]\n"
        "--params mt,split2,split1,mp,split2p,split1p,Nbase,Nghost,Noutgroup\n"
        "--epochs K: equal time intervals from 0 to split1\n"
        "--epoch-boundaries CSV: K+1 increasing times from 0 to split1; overrides equal intervals\n"
        "--ghost-sizes CSV: starting Nghost per epoch, youngest first (default: params Nghost)\n"
        "--fit-base: also estimate Nbase(t) on the same epochs (starts at params Nbase)\n"
        "--base-sizes CSV: starting Nbase per epoch; enables --fit-base\n"
        "--fixed-ghost: keep supplied Nghost(t) fixed during fitting (default: estimate it)\n"
        "  Above split1, Nbase remains fixed at params Nbase.\n"
        "--time-bounds LO,HI --prob-bounds LO,HI --N-bounds LO,HI\n"
        "  defaults: [0,split2) with outgroup, [0,split1) without; [0.000001,0.999999], [100,10000000]\n"
        "--iterations K (30) --maxeval K (100 per migration search) --tolerance X (1e-5 absolute log likelihood)\n"
        "--tree-num K: equally spaced genomic probes, deduplicated (default: all trees)\n"
        "--admixed-index I (0) --outgroup-index I (1)\n"
        "--cache-limit-mb X (1024 MiB, 0 unlimited); checked before allocating each tree\n"
        "--evaluate-only: evaluate supplied parameters without fitting\n"
        "--output FILE: iteration TSV; otherwise stdout. Progress goes to stderr.\n"
        "--profile-output FILE: final epochs and information status (default: OUTPUT.epochs.tsv)\n";
}
}
int main(int argc,char** argv) {
    try {
        std::string input,output,profile_output;std::vector<double> p,initial_sizes,initial_base,epoch_boundaries;
        std::vector<int> populations;size_t epochs=0,tree_limit=0;
        int admixed=0,outgroup=1;double memory_limit=1024;
        bool evaluate=false,has_time_bounds=false,fit_base=false;
        epochfit::FitOptions options;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            auto value=[&](){if(++i==argc)throw std::invalid_argument("missing value for "+arg);return std::string(argv[i]);};
            auto bounds=[&](){auto v=numbers(value());if(v.size()!=2)throw std::invalid_argument("bounds need LO,HI");return std::array<double,2>{v[0],v[1]};};
            if(arg=="--help" || arg=="-h") {help();return 0;}
            else if(arg=="--input")input=value();
            else if(arg=="--output")output=value();
            else if(arg=="--profile-output")profile_output=value();
            else if(arg=="--params")p=numbers(value());
            else if(arg=="--ghost-sizes")initial_sizes=numbers(value());
            else if(arg=="--base-sizes") {initial_base=numbers(value());fit_base=true;}
            else if(arg=="--fit-base")fit_base=true;
            else if(arg=="--fixed-ghost")options.fit_ghost=false;
            else if(arg=="--sample-population") {
                for(double x:numbers(value())) {if(x!=std::floor(x)||x<0||x>2147483647)throw std::invalid_argument("invalid population ID");populations.push_back(int(x));}
            }
            else if(arg=="--epochs")epochs=count(value());
            else if(arg=="--epoch-boundaries")epoch_boundaries=numbers(value());
            else if(arg=="--tree-num")tree_limit=count(value());
            else if(arg=="--iterations")options.iterations=count(value());
            else if(arg=="--maxeval")options.max_evaluations=count(value());
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
        if(input.empty() || p.size()!=9 || populations.empty() || !epochs)throw std::invalid_argument("require --input, --params (9 values), --sample-population and --epochs");
        if(epochs>100000)throw std::invalid_argument("--epochs must be at most 100000");
        if(!epoch_boundaries.empty() && epoch_boundaries.size()!=epochs+1)
            throw std::invalid_argument("--epoch-boundaries needs K+1 values");
        if(profile_output.empty() && !output.empty())profile_output=output+".epochs.tsv";
        if(!output.empty() && (output==input || output==profile_output))throw std::invalid_argument("input, iteration output and profile output must have different paths");
        if(!profile_output.empty() && profile_output==input)throw std::invalid_argument("profile output must differ from input");
        if(!std::isfinite(memory_limit)||memory_limit<0)throw std::invalid_argument("invalid memory limit");
        if(initial_sizes.empty())initial_sizes.assign(epochs,p[7]);
        if(initial_sizes.size()!=epochs)throw std::invalid_argument("--ghost-sizes needs one value per epoch");
        if(fit_base && initial_base.empty())initial_base.assign(epochs,p[6]);
        if(fit_base && initial_base.size()!=epochs)throw std::invalid_argument("--base-sizes needs one value per epoch");
        std::ofstream file;
        if(!output.empty()) {file.open(output);if(!file)throw std::runtime_error("cannot open output file");}
        std::ostream& out=output.empty()?std::cout:file;
        out<<std::setprecision(17);std::cerr<<std::setprecision(12);
        tsk_treeseq_t ts{};
        int ret=tsk_treeseq_load(&ts,input.c_str(),0);
        if(ret<0)throw std::runtime_error(tsk_strerror(ret));
        struct Cleanup {tsk_treeseq_t& ts;~Cleanup(){tsk_treeseq_free(&ts);}} cleanup{ts};
        if(!has_time_bounds) {
            bool has_outgroup=false;
            const auto* samples=tsk_treeseq_get_samples(&ts);
            for(size_t i=0;i<tsk_treeseq_get_num_samples(&ts);++i)
                if(size_t(samples[i])<populations.size() && populations[samples[i]]==outgroup) {has_outgroup=true;break;}
            options.time_bounds={0,std::nextafter(p[has_outgroup?1:2],0.)};
        }
        using Clock=std::chrono::steady_clock;
        auto start=Clock::now();int last=-1;
        epochfit::Ensemble ensemble(ts,p,populations,admixed,outgroup,epochs,tree_limit,memory_limit,
            [&](size_t done,size_t total,size_t bytes){
                int bucket=total?int(10.*done/total)*10:100;
                if(bucket>last) {last=bucket;std::cerr<<"Epoch cache: "<<bucket<<"% (scanned "<<done<<"/"<<total<<" trees), "<<bytes/(1024.*1024.)<<" MiB\n";}
            },epoch_boundaries);
        if(fit_base)ensemble.set_sizes(initial_base,initial_sizes);
        else ensemble.set_sizes(initial_sizes);
        std::cerr<<"Cache ready: selected trees="<<ensemble.tree_count()<<" scenarios="<<ensemble.scenario_count()
            <<" resident arrays="<<ensemble.bytes()/(1024.*1024.)<<" MiB; peak RSS="<<rss()<<" MiB; seconds="
            <<std::chrono::duration<double>(Clock::now()-start).count()<<'\n';
        std::cerr<<"Epoch boundaries (youngest first):";for(double b:ensemble.boundaries())std::cerr<<' '<<b;std::cerr<<'\n';
        out<<"iteration\tlog_llh\tT_adm\tp\timprovement";
        if(fit_base)for(size_t e=0;e<epochs;++e)out<<"\tNbase_"<<e;
        for(size_t e=0;e<epochs;++e)out<<"\tNghost_"<<e;
        out<<'\n';
        options.iteration=[&](size_t it,double llh,const epochfit::Parameters& state,double gain){
            out<<it<<'\t'<<llh<<'\t'<<state.time<<'\t'<<state.probability<<'\t'<<gain;
            for(double n:state.base_sizes)out<<'\t'<<n;
            for(double n:state.ghost_sizes)out<<'\t'<<n;out<<'\n';out.flush();
            std::cerr<<"ECME "<<it<<": log_llh="<<llh<<" improvement="<<gain<<" T_adm="<<state.time<<" p="<<state.probability;
            if(fit_base) {std::cerr<<" Nbase=";for(double n:state.base_sizes)std::cerr<<n<<',';}
            std::cerr<<" Nghost=";
            for(double n:state.ghost_sizes)std::cerr<<n<<',';std::cerr<<'\n';
        };
        options.optimizer=[&](size_t it,size_t calls,double best){
            if(calls==1 || calls%10==0 || calls==options.max_evaluations)
                std::cerr<<"Migration search "<<it<<": "<<calls<<'/'<<options.max_evaluations<<" evaluations, best="<<best<<'\n';
        };
        epochfit::Parameters state{p[0],p[3],initial_sizes,initial_base};
        if(evaluate) options.iteration(0,ensemble.evaluate(state.time,state.probability),state,0);
        else {
            const auto result=epochfit::fit(ensemble,state,options);
            state=result.parameters;
            std::cerr<<"Finished: "<<result.status<<"; iterations="<<result.iterations<<"; log_llh="<<result.log_likelihood<<'\n';
        }
        if(!profile_output.empty()) {
            std::ofstream profile(profile_output);if(!profile)throw std::runtime_error("cannot open profile output");
            profile<<std::setprecision(17)<<"epoch\tleft_time\tright_time";
            if(fit_base)profile<<"\tNbase\tbase_expected_coalescences\tbase_exposure\tbase_status";
            profile<<"\tNghost\texpected_coalescences\texposure\tstatus\n";
            const auto expected=ensemble.expectation(state.time,state.probability);
            for(size_t e=0;e<epochs;++e) {
                std::string status="active";
                if(ensemble.boundaries()[e+1]<=state.time)status="below_migration";
                else if(!options.fit_ghost)status="fixed";
                else if(expected.exposure[e]==0 && expected.coalescences[e]==0)status="no_information";
                else if(state.ghost_sizes[e]==options.max_size)status="upper_bound";
                else if(state.ghost_sizes[e]==options.min_size)status="lower_bound";
                profile<<e<<'\t'<<ensemble.boundaries()[e]<<'\t'<<ensemble.boundaries()[e+1];
                if(fit_base) {
                    std::string base_status="active";
                    if(expected.base_exposure[e]==0 && expected.base_coalescences[e]==0)base_status="no_information";
                    else if(state.base_sizes[e]==options.max_size)base_status="upper_bound";
                    else if(state.base_sizes[e]==options.min_size)base_status="lower_bound";
                    profile<<'\t'<<state.base_sizes[e]<<'\t'<<expected.base_coalescences[e]
                        <<'\t'<<expected.base_exposure[e]<<'\t'<<base_status;
                }
                profile<<'\t'<<state.ghost_sizes[e]
                    <<'\t'<<expected.coalescences[e]<<'\t'<<expected.exposure[e]<<'\t'<<status<<'\n';
            }
            if(!profile)throw std::runtime_error("failed writing profile output");
            std::cerr<<"Epoch profile: "<<profile_output<<'\n';
        }
        std::cerr<<"Total seconds="<<std::chrono::duration<double>(Clock::now()-start).count()<<"; peak RSS="<<rss()<<" MiB\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<"epoch-fit: "<<e.what()<<'\n';return 1;}
}
