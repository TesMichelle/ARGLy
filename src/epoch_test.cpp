#include "epoch_fit.h"
#include <likelihood_cache.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
void check(int r){if(r<0)throw std::runtime_error(tsk_strerror(r));}
void close(double a,double b,double tolerance=1e-8){if(a==b)return;if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance){std::cerr<<a<<" != "<<b<<'\n';throw std::runtime_error("epoch test mismatch");}}
void require(bool condition,const char* msg){if(!condition)throw std::runtime_error(msg);}
}
int main(){try{
    tsk_table_collection_t tables{};check(tsk_table_collection_init(&tables,0));tables.sequence_length=3;
    for(int i=0;i<6;++i)check(tsk_node_table_add_row(&tables.nodes,TSK_NODE_IS_SAMPLE,0,TSK_NULL,TSK_NULL,nullptr,0));
    for(double t:{300.,700.,400.,2200.,3000.,1300.,1500.})check(tsk_node_table_add_row(&tables.nodes,0,t,TSK_NULL,TSK_NULL,nullptr,0));
    std::vector<std::vector<std::pair<int,int>>> edges={
        {{6,0},{6,1},{7,2},{7,3},{8,4},{8,5},{9,6},{9,7},{10,9},{10,8}},
        {{6,0},{6,2},{7,1},{7,3},{8,4},{8,5},{11,6},{11,8},{10,11},{10,7}},
        {{6,0},{6,1},{7,2},{7,3},{8,4},{8,5},{12,6},{12,7},{10,12},{10,8}}};
    for(size_t i=0;i<edges.size();++i)for(auto [parent,child]:edges[i])check(tsk_edge_table_add_row(&tables.edges,i,i+1,parent,child,nullptr,0));
    check(tsk_table_collection_sort(&tables,nullptr,0));tsk_treeseq_t ts{};check(tsk_treeseq_init(&ts,&tables,TSK_TS_INIT_BUILD_INDEXES));
    std::vector<double> p={500,1000,2000,0.2,0.5,0.5,10000,12000,9000};std::vector<int> pops={0,0,0,0,1,1};
    treellh::TreeEnsemble legacy(ts,p,pops,0,1);legacy.build();
    size_t comparisons=0;
    for(size_t bins:{size_t(1),size_t(4),size_t(7)}){
        epochfit::Ensemble ensemble(ts,p,pops,0,1,bins);
        for(double n:{120.,12000.,30000.}){
            ensemble.set_sizes(std::vector<double>(bins,n));p[7]=n;
            for(double time:{0.,100.,std::nextafter(300.,0.),300.,std::nextafter(300.,400.),
                    301.,500.,699.,700.,701.,900.,999.})
                for(double probability:{0.,1e-10,0.2,0.8,1.-1e-10,1.}){
                p[0]=time;p[3]=probability;
                close(ensemble.evaluate(time,probability),legacy.evaluate(p));++comparisons;
            }
        }
    }
    for(double base:{5000.,10000.,20000.}) {
        p[6]=base;
        epochfit::Ensemble constant(ts,p,pops,0,1,4);
        constant.set_sizes(std::vector<double>(4,base),std::vector<double>(4,12000.));
        p[7]=12000.;
        for(double time:{0.,300.,500.,700.,999.})for(double probability:{0.05,0.2,0.8}) {
            p[0]=time;p[3]=probability;
            close(constant.evaluate(time,probability),legacy.evaluate(p));++comparisons;
        }
    }
    p[6]=10000.;
    epochfit::Ensemble ensemble(ts,p,pops,0,1,4);
    const std::vector<double> sizes={5000,8000,13000,20000};ensemble.set_sizes(sizes);
    epochfit::Ensemble threaded(ts,p,pops,0,1,4);
    threaded.set_threads(2);
    threaded.set_sizes(sizes);
    for(double time:{0.,std::nextafter(300.,0.),300.,std::nextafter(300.,400.),
            301.,500.,700.,701.,999.})
        for(double probability:{0.,1e-10,0.2,0.8,1.-1e-10,1.}) {
            close(threaded.evaluate(time,probability),ensemble.evaluate(time,probability));
            if(probability==0.2 || probability==0.8) {
                const auto serial=ensemble.expectation(time,probability);
                const auto parallel=threaded.expectation(time,probability);
                close(parallel.log_likelihood,serial.log_likelihood);
                for(size_t e=0;e<sizes.size();++e) {
                    close(parallel.exposure[e],serial.exposure[e]);
                    close(parallel.coalescences[e],serial.coalescences[e]);
                    close(parallel.base_exposure[e],serial.base_exposure[e]);
                    close(parallel.base_coalescences[e],serial.base_coalescences[e]);
                }
            }
        }
    const std::vector<double> changed_base={7000,9000,11000,15000};
    const std::vector<double> changed_ghost={9000,10000,11000,12000};
    ensemble.set_sizes(changed_base,changed_ghost);
    threaded.set_sizes(changed_base,changed_ghost);
    close(threaded.evaluate(550,0.3),ensemble.evaluate(550,0.3));
    ensemble.set_sizes(sizes);
    for(double time:{0.,300.,500.,701.,999.}){
        const auto expected=ensemble.expectation(time,0.2);
        close(expected.log_likelihood,ensemble.evaluate(time,0.2));
        for(size_t e=0;e<sizes.size();++e){
            auto perturbed=sizes;const double eps=1e-5;
            perturbed[e]=sizes[e]*std::exp(eps);ensemble.set_sizes(perturbed);double plus=ensemble.evaluate(time,0.2);
            perturbed[e]=sizes[e]*std::exp(-eps);ensemble.set_sizes(perturbed);double minus=ensemble.evaluate(time,0.2);
            close((plus-minus)/(2*eps),expected.exposure[e]/(2*sizes[e])-expected.coalescences[e],1e-7);
        }
        ensemble.set_sizes(sizes);
        auto updated=sizes;
        for(size_t e=0;e<sizes.size();++e)if(expected.exposure[e]>0 || expected.coalescences[e]>0)
            updated[e]=expected.coalescences[e]>0?std::clamp(expected.exposure[e]/(2*expected.coalescences[e]),100.,1e7):1e7;
        ensemble.set_sizes(updated);
        require(ensemble.evaluate(time,0.2)>=expected.log_likelihood-1e-8,"N M-step decreased likelihood");
        ensemble.set_sizes(sizes);
    }
    const std::vector<double> base_sizes={7000,9000,11000,15000};
    ensemble.set_sizes(base_sizes,sizes);
    {
        epochfit::Ensemble linear(ts,p,pops,0,1,4);
        linear.set_linear_sizes(base_sizes,base_sizes,sizes,sizes);
        for(double time:{300.,500.,701.,999.})
            close(linear.evaluate(time,0.2),ensemble.evaluate(time,0.2),1e-7);
        const std::vector<double> base_right={8000,8500,12000,13000};
        const std::vector<double> ghost_right={7000,9000,11000,18000};
        linear.set_linear_sizes(base_sizes,base_right,sizes,ghost_right);
        const double initial=linear.evaluate(550,0.2);
        const auto evidence=linear.expectation(550,0.2);
        close(evidence.log_likelihood,initial,1e-8);
        epochfit::Ensemble linear_parallel(ts,p,pops,0,1,4);
        linear_parallel.set_threads(2);
        linear_parallel.set_linear_sizes(base_sizes,base_right,sizes,ghost_right);
        close(linear_parallel.evaluate(550,0.2),initial,1e-8);
        const auto parallel_evidence=linear_parallel.expectation(550,0.2);
        for(size_t e=0;e<4;++e) {
            close(parallel_evidence.base_rates[e].left_exposure,
                evidence.base_rates[e].left_exposure,1e-7);
            close(parallel_evidence.ghost_rates[e].right_exposure,
                evidence.ghost_rates[e].right_exposure,1e-7);
        }
        const auto derivatives=linear.derivatives(550,0.2);
        close(derivatives.time,
            (linear.evaluate(550.01,0.2)-linear.evaluate(549.99,0.2))/0.02,1e-6);
        close(derivatives.probability,
            (linear.evaluate(550,0.20001)-linear.evaluate(550,0.19999))/0.00002,1e-6);
        for(size_t e=0;e<4;++e) {
            close(evidence.base_rates[e].left_exposure+evidence.base_rates[e].right_exposure,
                evidence.base_exposure[e],1e-7);
            close(evidence.ghost_rates[e].left_exposure+evidence.ghost_rates[e].right_exposure,
                evidence.exposure[e],1e-7);
        }
        auto trial=base_right;
        const double eps=1e-5;
        trial[1]=base_right[1]*std::exp(eps);
        linear.set_linear_sizes(base_sizes,trial,sizes,ghost_right);
        const double plus=linear.evaluate(550,0.2);
        trial[1]=base_right[1]*std::exp(-eps);
        linear.set_linear_sizes(base_sizes,trial,sizes,ghost_right);
        const double minus=linear.evaluate(550,0.2);
        const double lambda=0.5/base_right[1];
        double event_gradient=0;
        for(const auto& [z,mass]:evidence.base_rates[1].events)
            event_gradient+=mass*z/((1-z)*0.5/base_sizes[1]+z*lambda);
        close((plus-minus)/(2*eps),
            -lambda*(event_gradient-evidence.base_rates[1].right_exposure),1e-6);
        epochfit::FitOptions linear_options;
        linear_options.linear_rates=true;linear_options.time_bounds={550,550};
        linear_options.probability_bounds={0.2,0.2};linear_options.iterations=4;
        epochfit::Parameters start{550,0.2,sizes,base_sizes,ghost_right,base_right};
        double previous_linear=-std::numeric_limits<double>::infinity();
        linear_options.iteration=[&](size_t,double value,const epochfit::Parameters&,double) {
            require(value>=previous_linear-1e-8,"linear-rate ECME decreased likelihood");
            previous_linear=value;
        };
        const auto fitted=epochfit::fit(linear,start,linear_options);
        require(fitted.log_likelihood>=initial-1e-8,"linear-rate fit lost initial likelihood");
        close(fitted.log_likelihood,linear.evaluate(550,0.2),1e-7);
    }
    {
        epochfit::Ensemble split(ts,p,pops,0,1,4);
        const std::vector<double> cuts={250,750,1250,1750};
        split.set_split_sizes(base_sizes,base_sizes,cuts,sizes,sizes,cuts);
        for(double time:{300.,500.,701.,999.})
            close(split.evaluate(time,0.2),ensemble.evaluate(time,0.2),1e-7);
        const std::vector<double> base_right={8000,8500,12000,13000};
        const std::vector<double> ghost_right={7000,9000,11000,18000};
        split.set_split_sizes(base_sizes,base_right,cuts,sizes,ghost_right,cuts);
        const double initial=split.evaluate(550,0.2);
        const auto evidence=split.expectation(550,0.2);
        close(initial,evidence.log_likelihood,1e-8);
        epochfit::Ensemble split_parallel(ts,p,pops,0,1,4);
        split_parallel.set_threads(2);
        split_parallel.set_split_sizes(base_sizes,base_right,cuts,sizes,ghost_right,cuts);
        close(split_parallel.evaluate(550,0.2),initial,1e-8);
        const auto parallel_evidence=split_parallel.expectation(550,0.2);
        for(size_t e=0;e<4;++e) {
            close(parallel_evidence.base_exposure[e],evidence.base_exposure[e],1e-7);
            close(parallel_evidence.base_coalescences[e],evidence.base_coalescences[e],1e-7);
        }
        const auto split_derivatives=split.derivatives(550,0.2);
        close(split_derivatives.time,
            (split.evaluate(550.01,0.2)-split.evaluate(549.99,0.2))/0.02,1e-6);
        bool cut_derivative=false;
        try{split.derivatives(750,0.2);}
        catch(const std::invalid_argument&){cut_derivative=true;}
        require(cut_derivative,"time derivative accepted at movable cut");
        for(size_t e=0;e<4;++e) {
            close(evidence.base_rates[e].left_exposure+evidence.base_rates[e].right_exposure,
                evidence.base_exposure[e],1e-7);
            close(evidence.ghost_rates[e].left_exposure+evidence.ghost_rates[e].right_exposure,
                evidence.exposure[e],1e-7);
        }
        epochfit::FitOptions split_options;
        split_options.movable_cuts=true;split_options.time_bounds={550,550};
        split_options.probability_bounds={0.2,0.2};split_options.iterations=4;
        epochfit::Parameters start{550,0.2,sizes,base_sizes,ghost_right,base_right,cuts,cuts};
        double previous_split=-std::numeric_limits<double>::infinity();
        split_options.iteration=[&](size_t,double value,const epochfit::Parameters& state,double) {
            require(value>=previous_split-1e-8,"movable-cut ECME decreased likelihood");
            for(size_t e=0;e<4;++e)
                require(split.base_boundaries()[e]<state.base_cuts[e]
                    && state.base_cuts[e]<split.base_boundaries()[e+1],"cut left base epoch");
            previous_split=value;
        };
        const auto fitted=epochfit::fit(split,start,split_options);
        require(fitted.log_likelihood>=initial-1e-8,"movable-cut fit lost initial likelihood");
        close(fitted.log_likelihood,split.evaluate(550,0.2),1e-7);
    }
    {
        epochfit::Ensemble mixed(ts,p,pops,0,1,4);
        const std::vector<double> cuts={250,750,1250,1750};
        const std::vector<double> base_right={8000,8500,12000,13000};
        mixed.set_split_base_sizes(base_sizes,base_right,cuts,sizes);
        epochfit::Ensemble equivalent(ts,p,pops,0,1,4);
        equivalent.set_split_sizes(base_sizes,base_right,cuts,sizes,sizes,cuts);
        for(double time:{300.,550.,701.,999.})
            close(mixed.evaluate(time,0.2),equivalent.evaluate(time,0.2),1e-7);
        const double initial=mixed.evaluate(550,0.2);
        const auto evidence=mixed.expectation(550,0.2);
        close(initial,evidence.log_likelihood,1e-8);
        require(evidence.ghost_rates.empty(),"base-only cuts retained ghost rate evidence");
        require(evidence.base_rates.size()==base_sizes.size(),"missing base cut evidence");
        epochfit::FitOptions options;
        options.base_movable_cuts=true;options.time_bounds={550,550};
        options.probability_bounds={0.2,0.2};options.iterations=4;
        epochfit::Parameters start{550,0.2,sizes,base_sizes,{},base_right,{},cuts};
        double previous=-std::numeric_limits<double>::infinity();
        options.iteration=[&](size_t,double value,const epochfit::Parameters& state,double) {
            require(value>=previous-1e-8,"base-only movable-cut ECME decreased likelihood");
            require(state.ghost_cuts.empty() && state.ghost_right_sizes.empty(),
                "base-only movable-cut fit added a ghost cut");
            previous=value;
        };
        const auto fitted=epochfit::fit(mixed,start,options);
        require(fitted.log_likelihood>=initial-1e-8,"base-only movable-cut fit lost likelihood");
        close(fitted.log_likelihood,mixed.evaluate(550,0.2),1e-7);
    }
    for(double time:{0.,300.,500.,701.,999.}) {
        const auto expected=ensemble.expectation(time,0.2);
        close(expected.log_likelihood,ensemble.evaluate(time,0.2));
        for(size_t e=0;e<base_sizes.size();++e) {
            auto perturbed=base_sizes;const double eps=1e-5;
            perturbed[e]*=std::exp(eps);ensemble.set_sizes(perturbed,sizes);
            const double plus=ensemble.evaluate(time,0.2);
            perturbed[e]=base_sizes[e]*std::exp(-eps);ensemble.set_sizes(perturbed,sizes);
            const double minus=ensemble.evaluate(time,0.2);
            close((plus-minus)/(2*eps),expected.base_exposure[e]/(2*base_sizes[e])-expected.base_coalescences[e],1e-7);
        }
        ensemble.set_sizes(base_sizes,sizes);
        auto updated_base=base_sizes,updated_ghost=sizes;
        for(size_t e=0;e<sizes.size();++e) {
            if(expected.base_exposure[e]>0 || expected.base_coalescences[e]>0)
                updated_base[e]=expected.base_coalescences[e]>0
                    ?std::clamp(expected.base_exposure[e]/(2*expected.base_coalescences[e]),100.,1e7):1e7;
            if(expected.exposure[e]>0 || expected.coalescences[e]>0)
                updated_ghost[e]=expected.coalescences[e]>0
                    ?std::clamp(expected.exposure[e]/(2*expected.coalescences[e]),100.,1e7):1e7;
        }
        ensemble.set_sizes(updated_base,updated_ghost);
        require(ensemble.evaluate(time,0.2)>=expected.log_likelihood-1e-8,"two-population N M-step decreased likelihood");
        ensemble.set_sizes(base_sizes,sizes);
    }
    epochfit::FitOptions options;options.time_bounds={10,999};options.iterations=5;options.max_evaluations=30;
    require(options.restarts==10,"wrong default ECME restart count");
    options.restarts=2;
    double previous=-std::numeric_limits<double>::infinity();
    options.iteration=[&](size_t,double llh,const epochfit::Parameters&,double){require(llh>=previous-1e-8,"ECME decreased likelihood");previous=llh;};
    auto result=epochfit::fit(ensemble,{500,0.2,sizes},options);
    close(result.log_likelihood,ensemble.evaluate(result.parameters.time,result.parameters.probability));
    auto threaded_options=options;
    threaded_options.iteration={};
    const auto threaded_result=epochfit::fit(threaded,{500,0.2,sizes},threaded_options);
    close(threaded_result.log_likelihood,result.log_likelihood,1e-7);
    close(threaded_result.parameters.time,result.parameters.time,1e-7);
    close(threaded_result.parameters.probability,result.parameters.probability,1e-7);
    previous=-std::numeric_limits<double>::infinity();
    auto dual=epochfit::fit(ensemble,{500,0.2,sizes,base_sizes},options);
    close(dual.log_likelihood,ensemble.evaluate(dual.parameters.time,dual.parameters.probability));
    options.fit_ghost=false; previous=-std::numeric_limits<double>::infinity();
    auto fixed=epochfit::fit(ensemble,{500,0.2,sizes,base_sizes},options);
    require(fixed.parameters.ghost_sizes==sizes,"fixed ghost profile changed");
    close(fixed.log_likelihood,ensemble.evaluate(fixed.parameters.time,fixed.parameters.probability));
    epochfit::FitOptions multi;
    multi.time_bounds={10,999};multi.iterations=1;multi.max_evaluations=12;
    multi.fit_ghost=false;multi.restart_seed=42;multi.restarts=1;
    size_t calls=0;
    multi.optimizer=[&](size_t,size_t n,double){calls=n;};
    const auto one_start=epochfit::fit(ensemble,{500,0.2,sizes},multi);
    const size_t one_calls=calls;
    multi.restarts=10;calls=0;
    const auto ten_starts=epochfit::fit(ensemble,{500,0.2,sizes},multi);
    require(calls>one_calls,"extra ECME restarts did not run");
    require(ten_starts.log_likelihood>=one_start.log_likelihood-1e-10,
        "extra ECME restarts worsened the conditional maximum");
    const auto repeated=epochfit::fit(ensemble,{500,0.2,sizes},multi);
    close(repeated.log_likelihood,ten_starts.log_likelihood,1e-10);
    close(repeated.parameters.time,ten_starts.parameters.time,1e-10);
    close(repeated.parameters.probability,ten_starts.parameters.probability,1e-10);
    multi.probability_bounds={0.2,0.2};multi.restarts=3;
    const auto fixed_probability=epochfit::fit(ensemble,{500,0.2,sizes},multi);
    close(fixed_probability.parameters.probability,0.2);
    require(fixed_probability.parameters.time>=10 && fixed_probability.parameters.time<=999,
        "random ECME time start left the bounds");
    multi.restarts=0;
    bool bad_restarts=false;
    try{epochfit::fit(ensemble,{500,0.2,sizes},multi);}
    catch(const std::invalid_argument&){bad_restarts=true;}
    require(bad_restarts,"zero ECME restarts accepted");
    const auto balanced=epochfit::balanced_epoch_boundaries(
        {{100,1},{200,1},{300,1},{400,1}},2,500);
    require(balanced==std::vector<double>({0,250,500}),
        "balanced epoch boundary did not fall between coalescences");
    bool too_many_epochs=false;
    try{epochfit::balanced_epoch_boundaries({{100,2},{200,2}},3,500);}
    catch(const std::invalid_argument&){too_many_epochs=true;}
    require(too_many_epochs,"balanced grid split tied coalescence times");
    epochfit::Ensemble boundary_case(ts,p,pops,0,1,2,0,0,{},
        {0,700,2000},2,{0,300,2000});
    const auto boundary_expected=boundary_case.expectation(500,0.2);
    const auto boundary_masses=boundary_case.coalescence_masses(500,0.2);
    double low_base=0,low_ghost=0;
    for(const auto& [time,mass]:boundary_masses[0])if(time<=300)low_base+=mass;
    for(const auto& [time,mass]:boundary_masses[1])if(time<=700)low_ghost+=mass;
    close(boundary_expected.base_coalescences[0],low_base);
    close(boundary_expected.coalescences[0],low_ghost);
    const auto raw_counts=epochfit::raw_coalescence_counts(ts,2000);
    double raw_total=0;
    for(const auto& [time,count]:raw_counts)raw_total+=count;
    close(raw_total,11);
    close(raw_counts.front().first,300);
    close(raw_counts[1].second,3); // Three simultaneous outgroup coalescences.
    const auto raw_information=epochfit::raw_epoch_statistics(
        ts,{0,600,2000},{0,400,1000,2000});
    close(raw_information[1][0].coalescences,6); // Events at 400 belong below it.
    close(raw_information[1][0].lineage_time,6900); // 3 × (6 × 300 + 5 × 100).
    close(raw_information[1][0].pair_time,16500); // 3 × (15 × 300 + 10 × 100).
    close(raw_information[0][0].coalescences,6);
    close(raw_information[0][0].lineage_time,9300);
    close(raw_information[0][0].pair_time,20100);
    double through_700=0;
    for(const auto& [time,count]:epochfit::raw_coalescence_counts(ts,700))through_700+=count;
    close(through_700,9); // Events at the upper split boundary are included.
    const auto window_counts=epochfit::coalescences_per_epoch(
        {{100,1},{200,2},{300,3}}, {0,200,400});
    require(window_counts==std::vector<double>({3,3}),
        "coalescence on an epoch boundary was not assigned to the younger epoch");
    epochfit::Ensemble two_grids(ts,p,pops,0,1,3,0,0,{},
        {0,400,1000,2000},2,{0,600,2000});
    const std::vector<double> separate_base={7000,11000};
    const std::vector<double> separate_ghost={5000,9000,13000};
    two_grids.set_sizes({10000,10000},{12000,12000,12000});
    p[6]=10000;p[7]=12000;p[3]=0.2;
    for(double time:{300.,400.,500.,600.,701.}) {
        p[0]=time;
        close(two_grids.evaluate(time,0.2),legacy.evaluate(p));
    }
    two_grids.set_sizes(separate_base,separate_ghost);
    {
        const auto expected=two_grids.expectation(500,0.3);
        const auto kernel=two_grids.kernel_profile(500,0.3,1e6,5);
        auto totals=[](const std::vector<double>& values){
            double sum=0;for(double value:values)sum+=value;return sum;
        };
        const auto near=[](double actual,double target){
            return std::abs(actual-target)<=1e-5*std::max(1.,std::abs(target));
        };
        require(kernel[0].size()==5 && kernel[1].size()==5,
            "kernel profile omitted a requested time point");
        require(near(kernel[0][2].coalescences,totals(expected.base_coalescences))
            && near(kernel[0][2].exposure,totals(expected.base_exposure)),
            "kernel base statistics do not recover E-step totals for a wide kernel");
        require(near(kernel[1][2].coalescences,totals(expected.coalescences))
            && near(kernel[1][2].exposure,totals(expected.exposure)),
            "kernel ghost statistics do not recover E-step totals for a wide kernel");
        require(std::isnan(kernel[1][0].size) && kernel[1][0].coalescences==0,
            "ghost kernel profile was extrapolated below migration");
    }
    for(double time:{350.,500.,650.,850.})for(double probability:{0.1,0.3,0.8}) {
        const auto analytic=two_grids.derivatives(time,probability);
        close(analytic.log_likelihood,two_grids.evaluate(time,probability));
        const double dt=0.001,dp=1e-5;
        close(analytic.time,
            (two_grids.evaluate(time+dt,probability)-two_grids.evaluate(time-dt,probability))/(2*dt),1e-7);
        close(analytic.probability,
            (two_grids.evaluate(time,probability+dp)-two_grids.evaluate(time,probability-dp))/(2*dp),1e-6);
        double tree_llh=0,tree_dt=0,tree_dp=0;
        for(const auto& tree:two_grids.derivatives_by_tree(time,probability)) {
            tree_llh+=tree.log_likelihood;tree_dt+=tree.time;tree_dp+=tree.probability;
        }
        close(tree_llh,analytic.log_likelihood);close(tree_dt,analytic.time);close(tree_dp,analytic.probability);
    }
    bool derivative_boundary_rejected=false;
    try{two_grids.derivatives(400,0.3);}catch(const std::invalid_argument&){derivative_boundary_rejected=true;}
    require(derivative_boundary_rejected,"time derivative reported at a cache boundary");
    bool derivative_endpoint_rejected=false;
    try{two_grids.derivatives(500,0);}catch(const std::invalid_argument&){derivative_endpoint_rejected=true;}
    require(derivative_endpoint_rejected,"probability derivative reported at an endpoint");
    treellh::Scenario_Computer profiled_ratio(ts,p,pops,0,1);
    profiled_ratio.set_population_profiles(
        {{0,600,2000},separate_base},
        {{0,400,1000,2000},separate_ghost});
    tsk_tree_t ratio_tree{};check(tsk_tree_init(&ratio_tree,&ts,0));
    for(double time:{0.,300.,400.,500.,600.,700.,999.})
        for(double probability:{0.,0.2,1.}) {
            p[0]=time;p[3]=probability;
            double sum=0;
            int next;
            for(next=tsk_tree_first(&ratio_tree);next==TSK_TREE_OK;next=tsk_tree_next(&ratio_tree)) {
                const auto [base,ghost]=profiled_ratio.compute_llh_ratio(ratio_tree,p,0);
                const double total=base==-std::numeric_limits<double>::infinity()?ghost
                    :ghost==-std::numeric_limits<double>::infinity()?base
                    :std::max(base,ghost)+std::log1p(std::exp(-std::abs(base-ghost)));
                sum+=total;
            }
            check(next);
            close(sum,two_grids.evaluate(time,probability),1e-7);
        }
    tsk_tree_free(&ratio_tree);
    require(two_grids.base_boundaries().size()==3 && two_grids.boundaries().size()==4,
        "population epoch grids were merged");
    for(double time:{300.,400.,500.,600.,701.}) {
        const auto expected=two_grids.expectation(time,0.2);
        const auto masses=two_grids.coalescence_masses(time,0.2);
        double total_base=0,total_ghost=0;
        for(const auto& [event,weight]:masses[0])total_base+=weight;
        for(const auto& [event,weight]:masses[1])total_ghost+=weight;
        double count_base=0,count_ghost=0;
        for(double count:expected.base_coalescences)count_base+=count;
        for(double count:expected.coalescences)count_ghost+=count;
        close(total_base,count_base);close(total_ghost,count_ghost);
        for(size_t e=0;e<separate_base.size();++e) {
            auto perturbed=separate_base;const double eps=1e-5;
            perturbed[e]*=std::exp(eps);two_grids.set_sizes(perturbed,separate_ghost);
            const double plus=two_grids.evaluate(time,0.2);
            perturbed[e]=separate_base[e]*std::exp(-eps);two_grids.set_sizes(perturbed,separate_ghost);
            const double minus=two_grids.evaluate(time,0.2);
            close((plus-minus)/(2*eps),expected.base_exposure[e]/(2*separate_base[e])-expected.base_coalescences[e],1e-7);
        }
        for(size_t e=0;e<separate_ghost.size();++e) {
            auto perturbed=separate_ghost;const double eps=1e-5;
            perturbed[e]*=std::exp(eps);two_grids.set_sizes(separate_base,perturbed);
            const double plus=two_grids.evaluate(time,0.2);
            perturbed[e]=separate_ghost[e]*std::exp(-eps);two_grids.set_sizes(separate_base,perturbed);
            const double minus=two_grids.evaluate(time,0.2);
            close((plus-minus)/(2*eps),expected.exposure[e]/(2*separate_ghost[e])-expected.coalescences[e],1e-7);
        }
        two_grids.set_sizes(separate_base,separate_ghost);
    }
    epochfit::FitOptions two_grid_fit;
    two_grid_fit.time_bounds={10,999};two_grid_fit.iterations=2;
    two_grid_fit.max_evaluations=12;two_grid_fit.restarts=2;
    const double two_grid_initial=two_grids.evaluate(500,0.2);
    const auto two_grid_result=epochfit::fit(two_grids,
        {500,0.2,separate_ghost,separate_base},two_grid_fit);
    require(two_grid_result.parameters.base_sizes.size()==2
        && two_grid_result.parameters.ghost_sizes.size()==3
        && two_grid_result.log_likelihood>=two_grid_initial-1e-8,
        "two-grid ECME changed dimensions or decreased likelihood");
    ensemble.set_sizes(sizes);
    epochfit::Ensemble pooled(ts,p,pops,0,1,4,1);
    pooled.append(ts,p,pops,0,1,1);
    pooled.set_sizes(sizes);
    epochfit::Ensemble one_selected(ts,p,pops,0,1,4,1);
    one_selected.set_sizes(sizes);
    require(pooled.tree_count()==2*one_selected.tree_count(),"appended tree count is wrong");
    close(pooled.evaluate(500,0.2),2*one_selected.evaluate(500,0.2));
    const auto pooled_expected=pooled.expectation(500,0.2);
    const auto single_expected=one_selected.expectation(500,0.2);
    for(size_t e=0;e<sizes.size();++e) {
        close(pooled_expected.exposure[e],2*single_expected.exposure[e]);
        close(pooled_expected.coalescences[e],2*single_expected.coalescences[e]);
    }
    epochfit::Ensemble bounded(ts,p,pops,0,1,4,0,0,{}, {},0,{}, {300,700});
    bounded.set_sizes(sizes);
    require(bounded.bytes()<ensemble.bytes(),"bounded prefixes did not save memory");
    for(double time:{300.,400.,500.,699.,700.}) {
        close(bounded.evaluate(time,0.2),ensemble.evaluate(time,0.2));
        close(bounded.expectation(time,0.2).log_likelihood,
            ensemble.expectation(time,0.2).log_likelihood);
    }
    bool outside_bounds=false;
    try {bounded.evaluate(701,0.2);}catch(const std::invalid_argument&){outside_bounds=true;}
    require(outside_bounds,"bounded prefix cache accepted a time outside its range");
    bool rejected=false;try{epochfit::Ensemble limited(ts,p,pops,0,1,4,0,0.0001);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"memory limit ignored");
    size_t last_reported=0;
    epochfit::Ensemble monitored(ts,p,pops,0,1,4,0,ensemble.bytes()/(1024.*1024.)+0.001,
        [&](size_t,size_t,size_t bytes){require(bytes>=last_reported,"memory counter decreased during construction");last_reported=bytes;});
    require(monitored.bytes()>=last_reported,"final memory counter excludes constructed data");
    const size_t stable_bytes=monitored.bytes();
    monitored.set_sizes(sizes);
    require(monitored.bytes()==stable_bytes,"memory counter changed without a new allocation");
    tsk_treeseq_free(&ts);tsk_table_collection_free(&tables);
    // Independent closed form: two tips, root older than split1, four colourings.
    check(tsk_table_collection_init(&tables,0));tables.sequence_length=1;
    for(int i=0;i<2;++i)check(tsk_node_table_add_row(&tables.nodes,TSK_NODE_IS_SAMPLE,0,TSK_NULL,TSK_NULL,nullptr,0));
    check(tsk_node_table_add_row(&tables.nodes,0,3000,TSK_NULL,TSK_NULL,nullptr,0));
    for(int i=0;i<2;++i)check(tsk_edge_table_add_row(&tables.edges,0,1,2,i,nullptr,0));
    check(tsk_table_collection_sort(&tables,nullptr,0));check(tsk_treeseq_init(&ts,&tables,TSK_TS_INIT_BUILD_INDEXES));
    p={750,1000,2000,0.3,0.5,0.5,10000,12000,9000};
    epochfit::Ensemble simple(ts,p,{0,0},0,1,4);simple.set_sizes(sizes);
    const double ghost_integral=250/(2*sizes[1])+500/(2*sizes[2])+500/(2*sizes[3]);
    const double expected=std::log(0.5/10000)-(1000+750)/(2.*10000)
        +std::log(0.7*0.7*std::exp(-1250/(2.*10000))+2*0.3*0.7+0.3*0.3*std::exp(-ghost_integral));
    close(simple.evaluate(750,0.3),expected);
    const double late=std::log(0.5/10000)-(1000+1500)/(2.*10000)
        +std::log(0.7*0.7*std::exp(-500/(2.*10000))+2*0.3*0.7
            +0.3*0.3*std::exp(-500/(2*sizes[3])));
    close(simple.evaluate(1500,0.3),late);
    p[1]=500;
    epochfit::Ensemble placeholder(ts,p,{0,0},0,1,4);placeholder.set_sizes(sizes);
    close(placeholder.evaluate(1500,0.3),late);
    require(placeholder.migration_limit()==2000 && !placeholder.has_outgroup(),"no-outgroup migration limit depends on split2");
    bool time_rejected=false;try{placeholder.evaluate(2000,0.3);}catch(const std::invalid_argument&){time_rejected=true;}
    require(time_rejected,"migration at split1 must be rejected");
    epochfit::Ensemble custom(ts,p,{0,0},0,1,4,0,0,{}, {0,200,400,600,2000});
    custom.set_sizes(sizes);
    const double custom_expected=std::log(0.5/10000)-(1000+750)/(2.*10000)
        +std::log(0.7*0.7*std::exp(-1250/(2.*10000))+2*0.3*0.7
            +0.3*0.3*std::exp(-1250/(2*sizes[3])));
    close(custom.evaluate(750,0.3),custom_expected);
    bool boundaries_rejected=false;
    try {epochfit::Ensemble invalid(ts,p,{0,0},0,1,4,0,0,{}, {0,200,200,600,2000});}
    catch(const std::invalid_argument&){boundaries_rejected=true;}
    require(boundaries_rejected,"duplicate epoch boundary accepted");
    tsk_treeseq_free(&ts);tsk_table_collection_free(&tables);
    // A 256-sample tree requires the two-byte count representation.
    check(tsk_table_collection_init(&tables,0));tables.sequence_length=1;
    for(int i=0;i<256;++i)
        check(tsk_node_table_add_row(&tables.nodes,TSK_NODE_IS_SAMPLE,0,TSK_NULL,TSK_NULL,nullptr,0));
    for(int i=1;i<256;++i)
        check(tsk_node_table_add_row(&tables.nodes,0,i,TSK_NULL,TSK_NULL,nullptr,0));
    for(int i=1;i<256;++i) {
        const int parent=255+i;
        check(tsk_edge_table_add_row(&tables.edges,0,1,parent,i==1?0:parent-1,nullptr,0));
        check(tsk_edge_table_add_row(&tables.edges,0,1,parent,i,nullptr,0));
    }
    check(tsk_table_collection_sort(&tables,nullptr,0));
    check(tsk_treeseq_init(&ts,&tables,TSK_TS_INIT_BUILD_INDEXES));
    p={50,100,300,0.2,0.5,0.5,10000,12000,9000};
    epochfit::Ensemble wide_counts(ts,p,std::vector<int>(256,0),0,1,2,
        0,0,{}, {},0,{}, {50,50});
    require(std::isfinite(wide_counts.evaluate(50,0.2)),"two-byte count cache gave nonfinite likelihood");
    tsk_treeseq_free(&ts);tsk_table_collection_free(&tables);
    std::cout<<"Passed "<<comparisons<<" constant-N comparisons, two-population score derivatives, M-step and ECME monotonicity, memory guard\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
