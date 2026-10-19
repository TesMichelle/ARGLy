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
            for(double time:{0.,100.,300.,301.,500.,699.,700.,701.,900.,999.})for(double probability:{0.,0.2,0.8,1.}){
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
    double previous=-std::numeric_limits<double>::infinity();
    options.iteration=[&](size_t,double llh,const epochfit::Parameters&,double){require(llh>=previous-1e-8,"ECME decreased likelihood");previous=llh;};
    auto result=epochfit::fit(ensemble,{500,0.2,sizes},options);
    close(result.log_likelihood,ensemble.evaluate(result.parameters.time,result.parameters.probability));
    previous=-std::numeric_limits<double>::infinity();
    auto dual=epochfit::fit(ensemble,{500,0.2,sizes,base_sizes},options);
    close(dual.log_likelihood,ensemble.evaluate(dual.parameters.time,dual.parameters.probability));
    options.fit_ghost=false; previous=-std::numeric_limits<double>::infinity();
    auto fixed=epochfit::fit(ensemble,{500,0.2,sizes,base_sizes},options);
    require(fixed.parameters.ghost_sizes==sizes,"fixed ghost profile changed");
    close(fixed.log_likelihood,ensemble.evaluate(fixed.parameters.time,fixed.parameters.probability));
    ensemble.set_sizes(sizes);
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
    std::cout<<"Passed "<<comparisons<<" constant-N comparisons, two-population score derivatives, M-step and ECME monotonicity, memory guard\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
