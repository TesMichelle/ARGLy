#include "epoch_fit.h"
#include <nlopt.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
namespace epochfit {
FitResult fit(Ensemble& ensemble, Parameters p, const FitOptions& options) {
    if(!options.iterations || !options.max_evaluations || !options.restarts
        || options.max_evaluations>2147483647
        || options.restarts>std::numeric_limits<size_t>::max()/options.max_evaluations
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
    if((options.linear_rates && (options.movable_cuts || options.base_movable_cuts))
        || (options.movable_cuts && options.base_movable_cuts))
        throw std::invalid_argument("choose one epoch rate model");
    if(options.linear_rates || options.movable_cuts || options.base_movable_cuts) {
        if(p.base_sizes.empty() || p.base_right_sizes.size()!=p.base_sizes.size()
            || ((options.linear_rates || options.movable_cuts)
                && p.ghost_right_sizes.size()!=p.ghost_sizes.size()))
            throw std::invalid_argument("right epoch sizes do not match left sizes");
        for(double n:p.base_right_sizes)if(n<options.min_size || n>options.max_size)
            throw std::invalid_argument("starting base right N outside fit bounds");
        if(options.linear_rates || options.movable_cuts)
            for(double n:p.ghost_right_sizes)if(n<options.min_size || n>options.max_size)
                throw std::invalid_argument("starting ghost right N outside fit bounds");
    }
    if(options.movable_cuts && (p.base_cuts.size()!=p.base_sizes.size()
        || p.ghost_cuts.size()!=p.ghost_sizes.size()))
        throw std::invalid_argument("movable cuts need one point per epoch");
    if(options.base_movable_cuts && (p.base_cuts.size()!=p.base_sizes.size()
        || !p.ghost_cuts.empty() || !p.ghost_right_sizes.empty()))
        throw std::invalid_argument("base movable cuts need base points and constant ghost sizes");
    auto apply_sizes=[&]{
        if(options.linear_rates)ensemble.set_linear_sizes(p.base_sizes,p.base_right_sizes,
            p.ghost_sizes,p.ghost_right_sizes);
        else if(options.movable_cuts)ensemble.set_split_sizes(p.base_sizes,p.base_right_sizes,
            p.base_cuts,p.ghost_sizes,p.ghost_right_sizes,p.ghost_cuts);
        else if(options.base_movable_cuts)ensemble.set_split_base_sizes(
            p.base_sizes,p.base_right_sizes,p.base_cuts,p.ghost_sizes);
        else if(p.base_sizes.empty()) ensemble.set_sizes(p.ghost_sizes);
        else ensemble.set_sizes(p.base_sizes,p.ghost_sizes);
    };
    auto update_rate=[&](const RateEvidence& data,double& left_size,double& right_size) {
        if(data.events.empty() && data.left_exposure==0 && data.right_exposure==0)return;
        const double lo=0.5/options.max_size,hi=0.5/options.min_size;
        double rate[2]={0.5/left_size,0.5/right_size};
        for(size_t pass=0;pass<30;++pass) {
            const double old[2]={rate[0],rate[1]};
            for(size_t coordinate=0;coordinate<2;++coordinate) {
                auto derivative=[&](double trial) {
                    double value=-(coordinate?data.right_exposure:data.left_exposure);
                    for(const auto& [z,mass]:data.events) {
                        const double weight=coordinate?z:1-z;
                        value+=mass*weight/(weight*trial+(1-weight)*rate[1-coordinate]);
                    }
                    return value;
                };
                if(derivative(lo)<=0)rate[coordinate]=lo;
                else if(derivative(hi)>=0)rate[coordinate]=hi;
                else {
                    double low=lo,high=hi;
                    for(size_t k=0;k<55;++k) {
                        const double mid=std::midpoint(low,high);
                        if(derivative(mid)>0)low=mid;else high=mid;
                    }
                    rate[coordinate]=std::midpoint(low,high);
                }
            }
            if(std::abs(rate[0]-old[0])+std::abs(rate[1]-old[1])
                <1e-10*(rate[0]+rate[1]))break;
        }
        left_size=0.5/rate[0];right_size=0.5/rate[1];
    };
    auto update_cut=[&](const RateEvidence& data,double left,double right,
        double& cut,double& left_size,double& right_size) {
        if(data.cut_events.empty() && data.spans.empty())return;
        auto events=data.cut_events;
        std::sort(events.begin(),events.end());
        std::vector<std::pair<double,double>> changes;
        changes.reserve(data.spans.size()*2);
        double total_exposure=0,total_count=0;
        for(const auto& [time,mass]:events)total_count+=mass;
        for(const auto& span:data.spans) {
            changes.emplace_back(span[0],span[2]);
            changes.emplace_back(span[1],-span[2]);
            total_exposure+=(span[1]-span[0])*span[2];
        }
        std::sort(changes.begin(),changes.end());
        if(total_count==0) {
            if(total_exposure>0)left_size=right_size=options.max_size;
            return;
        }
        std::vector<double> candidates{std::nextafter(left,right),cut,std::nextafter(right,left)};
        candidates.reserve(events.size()*2+3);
        for(const auto& [time,mass]:events) {
            (void)mass;
            if(time>left && time<right) {
                candidates.push_back(time);
                const double before=std::nextafter(time,left);
                if(before>left)candidates.push_back(before);
            }
        }
        std::sort(candidates.begin(),candidates.end());
        candidates.erase(std::unique(candidates.begin(),candidates.end()),candidates.end());
        const double min_rate=0.5/options.max_size,max_rate=0.5/options.min_size;
        double best=-std::numeric_limits<double>::infinity();
        double best_left=left_size,best_right=right_size,best_cut=cut;
        size_t event_i=0,change_i=0;
        double cursor=left,active=0,exposure=0,count=0;
        for(double candidate:candidates) {
            if(!(candidate>left && candidate<right))continue;
            while(change_i<changes.size() && changes[change_i].first<=candidate) {
                const double next=changes[change_i].first;
                exposure+=active*(next-cursor);cursor=next;
                while(change_i<changes.size() && changes[change_i].first==next)
                    active+=changes[change_i++].second;
            }
            exposure+=active*(candidate-cursor);cursor=candidate;
            while(event_i<events.size() && events[event_i].first<=candidate)
                count+=events[event_i++].second;
            const double e0=std::clamp(exposure,0.,total_exposure);
            const double e1=std::max(0.,total_exposure-e0);
            const double c0=std::clamp(count,0.,total_count),c1=std::max(0.,total_count-c0);
            const double r0=e0>0?std::clamp(c0/e0,min_rate,max_rate)
                :c0>0?max_rate:0.5/left_size;
            const double r1=e1>0?std::clamp(c1/e1,min_rate,max_rate)
                :c1>0?max_rate:0.5/right_size;
            const double value=c0*std::log(r0)-e0*r0+c1*std::log(r1)-e1*r1;
            if(value>best) {
                best=value;best_cut=candidate;
                best_left=0.5/r0;best_right=0.5/r1;
            }
        }
        cut=best_cut;left_size=best_left;right_size=best_right;
    };
    apply_sizes();
    std::mt19937_64 rng(options.restart_seed);
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    double current=ensemble.evaluate(p.time,p.probability);
    if(!std::isfinite(current)) throw std::runtime_error("nonfinite starting likelihood");
    if(options.iteration) options.iteration(0,current,p,0);
    FitResult result{p,current,0,"iteration_limit"};
    for(size_t it=1;it<=options.iterations;++it) {
        const Parameters previous=p;
        const double old=current;
        const auto expected=ensemble.expectation(p.time,p.probability);
        if(options.movable_cuts || options.base_movable_cuts) {
            if(options.movable_cuts && options.fit_ghost)
                for(size_t e=0;e<p.ghost_sizes.size();++e)
                update_cut(expected.ghost_rates[e],ensemble.boundaries()[e],ensemble.boundaries()[e+1],
                    p.ghost_cuts[e],p.ghost_sizes[e],p.ghost_right_sizes[e]);
            if(options.fit_base)for(size_t e=0;e<p.base_sizes.size();++e)
                update_cut(expected.base_rates[e],ensemble.base_boundaries()[e],ensemble.base_boundaries()[e+1],
                    p.base_cuts[e],p.base_sizes[e],p.base_right_sizes[e]);
            if(options.base_movable_cuts && options.fit_ghost)
                for(size_t e=0;e<p.ghost_sizes.size();++e) {
                    if(expected.exposure[e]==0 && expected.coalescences[e]==0)continue;
                    p.ghost_sizes[e]=expected.coalescences[e]>0
                        ?std::clamp(expected.exposure[e]/(2*expected.coalescences[e]),
                            options.min_size,options.max_size)
                        :options.max_size;
                }
        } else if(options.linear_rates) {
            if(options.fit_ghost)for(size_t e=0;e<p.ghost_sizes.size();++e)
                update_rate(expected.ghost_rates[e],p.ghost_sizes[e],p.ghost_right_sizes[e]);
            if(options.fit_base)for(size_t e=0;e<p.base_sizes.size();++e)
                update_rate(expected.base_rates[e],p.base_sizes[e],p.base_right_sizes[e]);
        } else {
        if(options.fit_ghost) for(size_t e=0;e<p.ghost_sizes.size();++e) {
            if(expected.exposure[e]==0 && expected.coalescences[e]==0) continue;
            p.ghost_sizes[e]=expected.coalescences[e]>0
                ? std::clamp(expected.exposure[e]/(2*expected.coalescences[e]),options.min_size,options.max_size)
                : options.max_size;
        }
        if(options.fit_base)for(size_t e=0;e<p.base_sizes.size();++e) {
            if(expected.base_exposure[e]==0 && expected.base_coalescences[e]==0) continue;
            p.base_sizes[e]=expected.base_coalescences[e]>0
                ? std::clamp(expected.base_exposure[e]/(2*expected.base_coalescences[e]),options.min_size,options.max_size)
                : options.max_size;
        }
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
        std::vector<double> initial_step;
        double min_epoch=std::numeric_limits<double>::infinity();
        for(const auto* boundaries:{&ensemble.base_boundaries(),&ensemble.boundaries()})
            for(size_t e=1;e<boundaries->size();++e)
                min_epoch=std::min(min_epoch,(*boundaries)[e]-(*boundaries)[e-1]);
        for(size_t i=0;i<2;++i) if(context.hi[i]>context.lo[i]) {
            context.active.push_back(i);
            const double natural=i==0?min_epoch/4:0.01;
            initial_step.push_back(std::clamp(natural/(context.hi[i]-context.lo[i]),1e-5,0.05));
        }
        if(!context.active.empty()) {
            for(size_t restart=0;restart<options.restarts;++restart) {
                std::vector<double> x;
                x.reserve(context.active.size());
                for(size_t j=0;j<context.active.size();++j) x.push_back(uniform(rng));
                nlopt::opt optimizer(x.size()==1?nlopt::LN_NELDERMEAD:nlopt::LN_BOBYQA,x.size());
                optimizer.set_lower_bounds(0.);optimizer.set_upper_bounds(1.);
                optimizer.set_max_objective(Context::objective,&context);
                optimizer.set_maxeval(options.max_evaluations);optimizer.set_xtol_abs(1e-6);
                optimizer.set_initial_step(initial_step);
                double optimum=current;
                try {optimizer.optimize(x,optimum);}
                catch(const nlopt::roundoff_limited&) { /* Keep the best evaluated point. */ }
            }
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
