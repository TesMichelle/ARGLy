#include "epoch_fit.h"
#include <nlopt.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace epochfit {
FitResult fit(Ensemble& ensemble, Parameters p, const FitOptions& options) {
    if(!options.iterations || !options.max_evaluations || options.max_evaluations>2147483647
        || !(options.min_size>0 && options.max_size>=options.min_size)
        || !std::isfinite(options.max_size) || !(options.tolerance>0)
        || !std::isfinite(options.tolerance)) throw std::invalid_argument("invalid ECME options");
    auto tb=options.time_bounds, pb=options.probability_bounds;
    if(!(tb[0]>=0 && tb[0]<=tb[1] && tb[1]<ensemble.migration_limit())
        || !(pb[0]>=0 && pb[0]<=pb[1] && pb[1]<=1)) throw std::invalid_argument("invalid time/probability bounds");
    if(p.time<tb[0] || p.time>tb[1] || p.probability<pb[0] || p.probability>pb[1])
        throw std::invalid_argument("starting time/probability outside fit bounds");
    for(double n:p.ghost_sizes) if(n<options.min_size || n>options.max_size) throw std::invalid_argument("starting N outside fit bounds");
    for(double n:p.base_sizes) if(n<options.min_size || n>options.max_size) throw std::invalid_argument("starting base N outside fit bounds");
    if(!p.base_sizes.empty() && p.base_sizes.size()!=p.ghost_sizes.size())
        throw std::invalid_argument("base and ghost profiles need the same number of epochs");
    auto apply_sizes=[&]{if(p.base_sizes.empty()) ensemble.set_sizes(p.ghost_sizes);
        else ensemble.set_sizes(p.base_sizes,p.ghost_sizes);};
    apply_sizes();
    double current=ensemble.evaluate(p.time,p.probability);
    if(!std::isfinite(current)) throw std::runtime_error("nonfinite starting likelihood");
    if(options.iteration) options.iteration(0,current,p,0);
    FitResult result{p,current,0,"iteration_limit"};
    for(size_t it=1;it<=options.iterations;++it) {
        const Parameters previous=p;
        const double old=current;
        const auto expected=ensemble.expectation(p.time,p.probability);
        if(options.fit_ghost) for(size_t e=0;e<p.ghost_sizes.size();++e) {
            if(expected.exposure[e]==0 && expected.coalescences[e]==0) continue;
            p.ghost_sizes[e]=expected.coalescences[e]>0
                ? std::clamp(expected.exposure[e]/(2*expected.coalescences[e]),options.min_size,options.max_size)
                : options.max_size;
        }
        for(size_t e=0;e<p.base_sizes.size();++e) {
            if(expected.base_exposure[e]==0 && expected.base_coalescences[e]==0) continue;
            p.base_sizes[e]=expected.base_coalescences[e]>0
                ? std::clamp(expected.base_exposure[e]/(2*expected.base_coalescences[e]),options.min_size,options.max_size)
                : options.max_size;
        }
        apply_sizes();
        current=ensemble.evaluate(p.time,p.probability);
        struct Context {
            Ensemble& ensemble;
            const FitOptions& options;
            size_t iteration, calls=0;
            std::array<double,2> lo,hi,point;
            std::vector<size_t> active;
            double best;
            static double objective(const std::vector<double>& x,std::vector<double>&,void* ptr) {
                auto& c=*static_cast<Context*>(ptr);
                auto point=c.point;
                for(size_t j=0;j<x.size();++j) {size_t i=c.active[j];point[i]=c.lo[i]+x[j]*(c.hi[i]-c.lo[i]);}
                double value=c.ensemble.evaluate(point[0],point[1]);
                if(value>c.best) {c.best=value;c.point=point;}
                ++c.calls;
                if(c.options.optimizer) c.options.optimizer(c.iteration,c.calls,c.best);
                return std::isfinite(value)?value:-1e300;
            }
        } context{ensemble,options,it,0,{tb[0],pb[0]},{tb[1],pb[1]}, {p.time,p.probability},{},current};
        std::vector<double> x;
        for(size_t i=0;i<2;++i) if(context.hi[i]>context.lo[i]) {
            context.active.push_back(i);x.push_back((context.point[i]-context.lo[i])/(context.hi[i]-context.lo[i]));
        }
        if(!x.empty()) {
            nlopt::opt optimizer(x.size()==1?nlopt::LN_NELDERMEAD:nlopt::LN_BOBYQA,x.size());
            optimizer.set_lower_bounds(0.);optimizer.set_upper_bounds(1.);
            optimizer.set_max_objective(Context::objective,&context);
            optimizer.set_maxeval(options.max_evaluations);optimizer.set_xtol_abs(1e-6);
            double min_epoch=std::numeric_limits<double>::infinity();
            const auto& boundaries=ensemble.boundaries();
            for(size_t e=1;e<boundaries.size();++e)
                min_epoch=std::min(min_epoch,boundaries[e]-boundaries[e-1]);
            std::vector<double> initial_step;
            for(size_t i:context.active) {
                const double natural=i==0?min_epoch/4:0.01;
                initial_step.push_back(std::clamp(natural/(context.hi[i]-context.lo[i]),1e-5,0.05));
            }
            optimizer.set_initial_step(initial_step);
            double optimum=current;
            try {optimizer.optimize(x,optimum);} catch(const nlopt::roundoff_limited&) { /* Use the best evaluated point. */ }
        }
        p.time=context.point[0];p.probability=context.point[1];current=context.best;
        if(!std::isfinite(current) || current<old-1e-8*(1+std::abs(old)))
            throw std::runtime_error("ECME likelihood decreased: numerical/model inconsistency");
        if(current<old) {p=previous;current=old;apply_sizes();}
        if(options.iteration) options.iteration(it,current,p,current-old);
        result={p,current,it,"iteration_limit"};
        if(current-old<=options.tolerance) {result.status="converged";break;}
    }
    return result;
}
}
