#include "epoch_model.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace epochfit {
namespace {
constexpr double neginf = -std::numeric_limits<double>::infinity();
double pairs(int k) { return double(k) * (k-1) / 2; }
double add(double a, double b) {
    if (a == neginf) return b;
    if (b == neginf) return a;
    return std::max(a,b) + std::log1p(std::exp(-std::abs(a-b)));
}
void check(int r) { if (r < 0) throw std::runtime_error(tsk_strerror(r)); }
struct Guard {
    tsk_tree_t tree{};
    Guard(const tsk_treeseq_t& ts) { check(tsk_tree_init(&tree,&ts,0)); }
    ~Guard() { tsk_tree_free(&tree); }
};
size_t product(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max()/b)
        throw std::overflow_error("epoch cache size overflow");
    return a*b;
}
size_t tree_storage_bytes(const TreeCache& tree) {
    return tree.times.capacity()*sizeof(double)
        + tree.epochs.capacity()*sizeof(std::array<size_t,2>)
        + tree.steps.capacity()*sizeof(Step)
        + tree.prefix.capacity()*sizeof(double);
}
template<class Cost, class Work>
void for_trees(const std::vector<TreeCache>& trees,size_t requested,
                 Cost cost,Work work) {
    size_t threads=requested;
    if(!threads) threads=std::min<size_t>(8,std::max(1u,std::thread::hardware_concurrency()));
    threads=std::min(threads,trees.size());
    if(threads<2) {
        for(size_t i=0;i<trees.size();++i)work(0,i);
        return;
    }
    std::vector<size_t> cumulative(trees.size()+1,0);
    for(size_t i=0;i<trees.size();++i) {
        const size_t weight=cost(trees[i]);
        if(weight>std::numeric_limits<size_t>::max()-cumulative[i])
            throw std::overflow_error("epoch parallel work size overflow");
        cumulative[i+1]=cumulative[i]+weight;
    }
    // For small caches, launching workers costs more than the tree loop.
    if(!requested && cumulative.back()<3000) {
        for(size_t i=0;i<trees.size();++i)work(0,i);
        return;
    }
    std::vector<size_t> boundaries(threads+1);
    boundaries.front()=0;boundaries.back()=trees.size();
    for(size_t w=1;w<threads;++w) {
        const size_t target=cumulative.back()/threads*w
            +(cumulative.back()%threads)*w/threads;
        const auto first=cumulative.begin()+boundaries[w-1]+1;
        const auto last=cumulative.end()-(threads-w);
        boundaries[w]=std::clamp<size_t>(
            std::lower_bound(first,last,target)-cumulative.begin(),
            boundaries[w-1]+1,trees.size()-(threads-w));
    }
    std::exception_ptr failure;
    std::mutex failure_mutex;
    auto run=[&](size_t w) {
        try {
            for(size_t i=boundaries[w];i<boundaries[w+1];++i)work(w,i);
        } catch(...) {
            std::lock_guard lock(failure_mutex);
            if(!failure)failure=std::current_exception();
        }
    };
    {
        std::vector<std::jthread> workers;
        workers.reserve(threads-1);
        for(size_t w=1;w<threads;++w)workers.emplace_back(run,w);
        run(0);
    }
    if(failure)std::rethrow_exception(failure);
}
}

std::vector<double> balanced_epoch_boundaries(
    const std::vector<std::pair<double,double>>& events,
    size_t epochs, double split_time) {
    if(!epochs || !std::isfinite(split_time) || !(split_time>0))
        throw std::invalid_argument("invalid balanced epoch request");
    if(epochs==1)return {0,split_time};
    if(events.size()<epochs)
        throw std::invalid_argument("more epochs requested than distinct informative coalescence times");
    std::vector<double> cumulative(events.size()+1,0);
    for(size_t i=0;i<events.size();++i) {
        const auto [time,mass]=events[i];
        if(!std::isfinite(time) || !(time>0) || time>split_time
            || (i && !(time>events[i-1].first))
            || !std::isfinite(mass) || !(mass>0))
            throw std::invalid_argument("invalid coalescence times or weights for balanced epochs");
        cumulative[i+1]=cumulative[i]+mass;
    }
    if(!std::isfinite(cumulative.back()))
        throw std::invalid_argument("coalescence weights overflow");
    std::vector<double> boundaries{0};
    size_t previous=0;
    for(size_t e=1;e<epochs;++e) {
        const size_t first=previous+1, last=events.size()-(epochs-e);
        const double target=cumulative.back()*double(e)/epochs;
        const auto found=std::lower_bound(cumulative.begin()+first,cumulative.begin()+last+1,target);
        size_t best=found==cumulative.begin()+last+1?last:size_t(found-cumulative.begin());
        if(best>first && std::abs(cumulative[best-1]-target)<=std::abs(cumulative[best]-target))--best;
        const double boundary=std::midpoint(events[best-1].first,events[best].first);
        if(!(boundary>boundaries.back()) || !(boundary<split_time)
            || !(boundary>events[best-1].first) || !(boundary<events[best].first))
            throw std::invalid_argument("no representable boundary between coalescence times");
        boundaries.push_back(boundary);
        previous=best;
    }
    boundaries.push_back(split_time);
    return boundaries;
}

std::vector<std::pair<double,double>> raw_coalescence_counts(
    const tsk_treeseq_t& ts, double split_time) {
    if(!std::isfinite(split_time) || !(split_time>0))
        throw std::invalid_argument("invalid split time for coalescence counts");
    Guard guard(ts);
    std::map<double,double> counts;
    int ret;
    for(ret=tsk_tree_first(&guard.tree);ret==TSK_TREE_OK;ret=tsk_tree_next(&guard.tree)) {
        if(tsk_tree_get_num_roots(&guard.tree)!=1)
            throw std::invalid_argument("balanced epochs require single-root trees");
        std::vector<tsk_id_t> stack{guard.tree.left_child[guard.tree.virtual_root]};
        while(!stack.empty()) {
            const auto id=stack.back();stack.pop_back();
            const auto child=guard.tree.left_child[id];
            if(child==TSK_NULL)continue;
            if(guard.tree.right_sib[child]==TSK_NULL || guard.tree.right_sib[guard.tree.right_sib[child]]!=TSK_NULL)
                throw std::invalid_argument("balanced epochs require binary trees");
            const double time=ts.tables->nodes.time[id];
            if(time>0 && time<=split_time)++counts[time];
            for(auto ch=child;ch!=TSK_NULL;ch=guard.tree.right_sib[ch])stack.push_back(ch);
        }
    }
    check(ret);
    return {counts.begin(),counts.end()};
}

std::vector<double> coalescences_per_epoch(
    const std::vector<std::pair<double,double>>& events,
    const std::vector<double>& boundaries) {
    if(boundaries.size()<2)throw std::invalid_argument("epoch boundaries are empty");
    std::vector<double> counts(boundaries.size()-1,0);
    for(const auto& [time,count]:events) {
        const auto found=std::lower_bound(boundaries.begin()+1,boundaries.end(),time);
        if(time<=boundaries.front() || found==boundaries.end())
            throw std::invalid_argument("coalescence time outside epoch boundaries");
        counts[size_t(found-boundaries.begin()-1)]+=count;
    }
    return counts;
}

std::array<std::vector<RawEpochStatistics>,2> raw_epoch_statistics(
    const tsk_treeseq_t& ts, const std::vector<double>& base_boundaries,
    const std::vector<double>& ghost_boundaries) {
    const std::array<const std::vector<double>*,2> grids{&base_boundaries,&ghost_boundaries};
    std::array<std::vector<RawEpochStatistics>,2> result;
    for(size_t q=0;q<2;++q) {
        const auto& bounds=*grids[q];
        if(bounds.size()<2 || bounds.front()!=0)
            throw std::invalid_argument("raw epoch statistics need boundaries starting at zero");
        for(size_t e=1;e<bounds.size();++e)
            if(!std::isfinite(bounds[e]) || !(bounds[e]>bounds[e-1]))
                throw std::invalid_argument("raw epoch statistics need increasing finite boundaries");
        result[q].resize(bounds.size()-1);
    }
    const double upper=std::max(base_boundaries.back(),ghost_boundaries.back());
    auto add_span=[&](double left,double right,int lineages) {
        if(lineages<0)throw std::invalid_argument("negative raw lineage count");
        for(size_t q=0;q<2;++q) {
            const auto& bounds=*grids[q];
            double start=std::max(0.,left),end=std::min(right,bounds.back());
            while(start<end) {
                const size_t e=std::upper_bound(bounds.begin(),bounds.end(),start)-bounds.begin()-1;
                const double stop=std::min(end,bounds[e+1]);
                const double dt=stop-start;
                result[q][e].lineage_time+=lineages*dt;
                result[q][e].pair_time+=pairs(lineages)*dt;
                start=stop;
            }
        }
    };
    struct Event {double time;int delta;bool coalescence;};
    Guard guard(ts);
    int ret;
    for(ret=tsk_tree_first(&guard.tree);ret==TSK_TREE_OK;ret=tsk_tree_next(&guard.tree)) {
        if(tsk_tree_get_num_roots(&guard.tree)!=1)
            throw std::invalid_argument("raw epoch statistics require single-root trees");
        std::vector<Event> events;
        std::vector<tsk_id_t> stack{guard.tree.left_child[guard.tree.virtual_root]};
        while(!stack.empty()) {
            const auto id=stack.back();stack.pop_back();
            const auto child=guard.tree.left_child[id];
            const double time=ts.tables->nodes.time[id];
            if(child==TSK_NULL) {
                events.push_back({time,1,false});
            } else {
                if(guard.tree.right_sib[child]==TSK_NULL
                    || guard.tree.right_sib[guard.tree.right_sib[child]]!=TSK_NULL)
                    throw std::invalid_argument("raw epoch statistics require binary trees");
                events.push_back({time,-1,true});
                for(auto ch=child;ch!=TSK_NULL;ch=guard.tree.right_sib[ch])stack.push_back(ch);
            }
        }
        std::sort(events.begin(),events.end(),[](const Event& a,const Event& b){return a.time<b.time;});
        int lineages=0;
        double previous=0;
        for(const auto& event:events) {
            if(event.time>upper)break;
            add_span(previous,event.time,lineages);
            if(event.coalescence && event.time>0)for(size_t q=0;q<2;++q) {
                const auto& bounds=*grids[q];
                const auto found=std::lower_bound(bounds.begin()+1,bounds.end(),event.time);
                if(found!=bounds.end())
                    ++result[q][size_t(found-bounds.begin()-1)].coalescences;
            }
            lineages+=event.delta;
            previous=event.time;
        }
        add_span(previous,upper,lineages);
    }
    check(ret);
    return result;
}

Ensemble::Ensemble(const tsk_treeseq_t& ts, const std::vector<double>& p,
    const std::vector<int>& pops, int admixed, int outgroup, size_t epoch_count,
    size_t limit, double memory_limit_mb,
    const std::function<void(size_t,size_t,size_t)>& progress,
    const std::vector<double>& epoch_boundaries,
    size_t base_epoch_count, const std::vector<double>& base_epoch_boundaries)
    : base_size_(p.at(6)), split2_(p.at(1)), split1_(p.at(2)) {
    if (p.size()!=9 || !epoch_count || epoch_count>100000 || admixed==outgroup
        || !(base_size_>0) || !(p[7]>0) || !(split1_>0))
        throw std::invalid_argument("invalid epoch model parameters");
    if(base_epoch_count==0)base_epoch_count=epoch_count;
    if(base_epoch_count>100000)throw std::invalid_argument("too many base epochs");
    for (double x:p) if (!std::isfinite(x)) throw std::invalid_argument("parameters must be finite");
    auto make_boundaries=[&](size_t count,const std::vector<double>& supplied) {
        std::vector<double> result;
        if(supplied.empty()) {
            result.resize(count+1);
            for(size_t e=0;e<=count;++e)result[e]=split1_*double(e)/count;
        } else {
            if(supplied.size()!=count+1 || supplied.front()!=0 || supplied.back()!=split1_)
                throw std::invalid_argument("epoch boundaries must have K+1 values from 0 to split1");
            result=supplied;
        }
        for(size_t e=1;e<result.size();++e)
            if(!std::isfinite(result[e]) || !(result[e]>result[e-1]))
                throw std::invalid_argument("epoch boundaries must be finite and strictly increasing");
        return result;
    };
    boundaries_=make_boundaries(epoch_count,epoch_boundaries);
    base_boundaries_=base_epoch_boundaries.empty() && base_epoch_count==epoch_count
        ?boundaries_:make_boundaries(base_epoch_count,base_epoch_boundaries);
    resident_bytes_=sizeof(*this)+(boundaries_.capacity()+base_boundaries_.capacity())*sizeof(double);
    const auto* samples=tsk_treeseq_get_samples(&ts);
    const size_t ns=tsk_treeseq_get_num_samples(&ts), nn=tsk_treeseq_get_num_nodes(&ts);
    for (size_t i=0;i<ns;++i)
        if (size_t(samples[i])>=pops.size() || (pops[samples[i]]!=admixed && pops[samples[i]]!=outgroup))
            throw std::invalid_argument("sample populations must identify every sample as admixed or outgroup");
    for(size_t i=0;i<ns;++i) if(pops[samples[i]]==outgroup) {has_outgroup_=true;break;}
    if(has_outgroup_ && !(0<split2_ && split2_<split1_))
        throw std::invalid_argument("outgroup split2 must lie between 0 and split1");
    treellh::Scenario_Computer helper(ts,p,pops,admixed,outgroup);
    Guard guard(ts);
    size_t probe=0, visited=0;
    const size_t target=limit ? std::min(limit,size_t(ts.num_trees)) : size_t(ts.num_trees);
    if(memory_limit_mb>0 && (static_cast<long double>(bytes())+product(target,sizeof(TreeCache)))/(1024*1024)>memory_limit_mb)
        throw std::runtime_error("epoch cache metadata exceeds --cache-limit-mb");
    trees_.reserve(target);
    resident_bytes_+=trees_.capacity()*sizeof(TreeCache);
    if (progress) progress(0,ts.num_trees,bytes());
    int ret;
    for (ret=tsk_tree_first(&guard.tree);ret==TSK_TREE_OK;ret=tsk_tree_next(&guard.tree)) {
        const auto& tree=guard.tree;
        bool selected=!limit || limit>=size_t(ts.num_trees);
        if (!selected) {
            while (probe<target && (static_cast<long double>(probe)+0.5L)*ts.tables->sequence_length/target < tree.interval.right) {
                selected=true;
                ++probe;
            }
        }
        ++visited;
        if (!selected) { if(progress) progress(visited,ts.num_trees,bytes()); continue; }
        if (tsk_tree_get_num_roots(&tree)!=1) throw std::invalid_argument("epoch model requires single-root trees");
        std::vector<tsk_id_t> nodes, stack;
        stack.push_back(tree.left_child[tree.virtual_root]);
        while (!stack.empty()) {
            auto id=stack.back(); stack.pop_back(); nodes.push_back(id);
            size_t children=0;
            for(auto ch=tree.left_child[id];ch!=TSK_NULL;ch=tree.right_sib[ch]) { stack.push_back(ch); ++children; }
            if ((children!=0 && children!=2) || (children==0 && (!(ts.tables->nodes.flags[id]&TSK_NODE_IS_SAMPLE) || ts.tables->nodes.time[id]!=0))
                || (children && (ts.tables->nodes.flags[id]&TSK_NODE_IS_SAMPLE)))
                throw std::invalid_argument("epoch model requires binary trees and present-day sample leaves");
        }
        std::sort(nodes.begin(),nodes.end(),[&](auto a,auto b){return ts.tables->nodes.time[a]>ts.tables->nodes.time[b];});
        // These helpers use tree topology and node times directly; the local
        // event order has already been sorted into nodes above.
        helper.compute_lower_nodes_split_1(tree); helper.compute_fixed_nodes(tree);
        TreeCache cache;
        cache.impossible=has_outgroup_ && helper.check_imposible_split_2(tree)!=0;
        if (cache.impossible) { trees_.push_back(std::move(cache)); if(progress)progress(visited,ts.num_trees,bytes()); continue; }
        const auto frontier=helper.get_lower_nodes_split_1();
        std::vector<bool> fixed(nn,false);
        int out_lines=0;
        for(size_t i=0;i<ns;++i) if(pops[samples[i]]==outgroup) {
            auto id=samples[i];
            while(id!=TSK_NULL && !fixed[id]) { fixed[id]=true; id=tree.parent[id]; }
            if(id==TSK_NULL || ts.tables->nodes.time[id]>=split2_) ++out_lines;
        }
        size_t bits=0; for(auto id:frontier) if(!fixed[id]) ++bits;
        if(bits>=63 || bits>=std::numeric_limits<size_t>::digits) throw std::overflow_error("too many epoch colouring scenarios");
        cache.scenarios=size_t{1}<<bits;
        cache.times=boundaries_;
        cache.times.insert(cache.times.end(),base_boundaries_.begin(),base_boundaries_.end());
        if(has_outgroup_)cache.times.push_back(split2_);
        for(auto id:nodes) if(ts.tables->nodes.time[id]>0 && ts.tables->nodes.time[id]<=split1_) cache.times.push_back(ts.tables->nodes.time[id]);
        std::sort(cache.times.begin(),cache.times.end(),std::greater<double>());
        cache.times.erase(std::unique(cache.times.begin(),cache.times.end()),cache.times.end());
        cache.epochs.resize(cache.times.size());
        for(size_t j=0;j+1<cache.times.size();++j) {
            const double mid=(cache.times[j]+cache.times[j+1])/2;
            cache.epochs[j][0]=std::upper_bound(base_boundaries_.begin(),base_boundaries_.end(),mid)-base_boundaries_.begin()-1;
            cache.epochs[j][1]=std::upper_bound(boundaries_.begin(),boundaries_.end(),mid)-boundaries_.begin()-1;
        }
        const size_t records=product(cache.scenarios,cache.times.size());
        const size_t allocation=product(records,sizeof(Step)+sizeof(double));
        if(memory_limit_mb>0 && (static_cast<long double>(bytes())+allocation
            +cache.times.capacity()*sizeof(double)+cache.epochs.capacity()*sizeof(std::array<size_t,2>)
            +(epoch_count+base_epoch_count)*sizeof(double))/(1024*1024)>memory_limit_mb)
            throw std::runtime_error("epoch cache exceeds --cache-limit-mb before allocation");
        cache.steps.resize(records);
        cache.prefix.resize(records);
        double previous=std::max(split1_,ts.tables->nodes.time[nodes.front()]);
        int root_k=1;
        for(auto id:nodes) {
            double t=ts.tables->nodes.time[id]; if(t<=split1_) break;
            cache.root_log_likelihood-=pairs(root_k)*(previous-t)/(2*base_size_);
            cache.root_log_likelihood+=std::log(0.5/base_size_); ++root_k; previous=t;
        }
        cache.root_log_likelihood-=pairs(root_k)*(previous-split1_)/(2*base_size_);
        std::vector<int> colour(nn,0);
        for(size_t s=0;s<cache.scenarios;++s) {
            std::fill(colour.begin(),colour.end(),0);
            for(auto id:nodes) if(fixed[id] && ts.tables->nodes.time[id]<split2_) colour[id]=2;
            std::array<int,2> k{};
            size_t bit=0;
            for(auto id:frontier) {
                int pop=fixed[id]?0:int((s>>bit++)&1); ++k[pop];
                stack={id};
                while(!stack.empty()) {
                    auto v=stack.back(); stack.pop_back(); if(colour[v]==2) continue;
                    colour[v]=pop;
                    for(auto ch=tree.left_child[v];ch!=TSK_NULL;ch=tree.right_sib[ch]) stack.push_back(ch);
                }
            }
            size_t node_i=0;
            while(node_i<nodes.size() && ts.tables->nodes.time[nodes[node_i]]>split1_) ++node_i;
            for(size_t j=0;j<cache.times.size();++j) {
                auto& step=cache.steps[j*cache.scenarios+s];
                if(has_outgroup_ && cache.times[j]==split2_) k[0]-=out_lines;
                while(node_i<nodes.size() && ts.tables->nodes.time[nodes[node_i]]==cache.times[j]) {
                    auto id=nodes[node_i++];
                    if(cache.times[j]>0 && colour[id]!=2) { ++step.coalescences[colour[id]]; ++k[colour[id]]; }
                }
                step.lineages=k;
                if(k[0]<0 || k[1]<0) throw std::logic_error("negative lineage count in epoch cache");
            }
        }
        const size_t cache_bytes=tree_storage_bytes(cache);
        trees_.push_back(std::move(cache));
        resident_bytes_+=cache_bytes;
        if(progress) progress(visited,ts.num_trees,bytes());
    }
    check(ret);
    if(trees_.empty()) throw std::invalid_argument("no trees selected");
    set_sizes(std::vector<double>(epoch_count,p[7]));
}

void Ensemble::set_sizes(const std::vector<double>& sizes) {
    set_sizes(std::vector<double>(base_boundaries_.size()-1,base_size_),sizes);
}
void Ensemble::set_threads(size_t threads) {
    if(threads>64)throw std::invalid_argument("epoch worker count must be at most 64");
    threads_=threads;
}
void Ensemble::set_sizes(const std::vector<double>& base, const std::vector<double>& ghost) {
    set_profile(base,base,ghost,ghost,{},{},false,false,false);
}
void Ensemble::set_linear_sizes(const std::vector<double>& base_left,
    const std::vector<double>& base_right,const std::vector<double>& ghost_left,
    const std::vector<double>& ghost_right) {
    set_profile(base_left,base_right,ghost_left,ghost_right,{},{},true,false,false);
}
void Ensemble::set_split_sizes(const std::vector<double>& base_left,
    const std::vector<double>& base_right,const std::vector<double>& base_cuts,
    const std::vector<double>& ghost_left,const std::vector<double>& ghost_right,
    const std::vector<double>& ghost_cuts) {
    set_profile(base_left,base_right,ghost_left,ghost_right,base_cuts,ghost_cuts,false,true,true);
}
void Ensemble::set_split_base_sizes(const std::vector<double>& base_left,
    const std::vector<double>& base_right,const std::vector<double>& base_cuts,
    const std::vector<double>& ghost) {
    set_profile(base_left,base_right,ghost,ghost,base_cuts,{},false,true,false);
}
double Ensemble::rate(size_t pop,size_t epoch,double time) const {
    const auto& left=pop? sizes_:base_sizes_;
    const auto& right=pop? right_sizes_:base_right_sizes_;
    if(pop?split_ghost_:split_base_)return time<=(pop?cuts_[epoch]:base_cuts_[epoch])
        ?0.5/left[epoch]:0.5/right[epoch];
    if(!linear_rates_)return 0.5/left[epoch];
    const auto& bounds=pop?boundaries_:base_boundaries_;
    const double z=std::clamp((time-bounds[epoch])/(bounds[epoch+1]-bounds[epoch]),0.,1.);
    return (1-z)*(0.5/left[epoch])+z*(0.5/right[epoch]);
}
double Ensemble::rate_integral(size_t pop,size_t epoch,double low,double high) const {
    if(pop?split_ghost_:split_base_) {
        const double cut=pop?cuts_[epoch]:base_cuts_[epoch];
        const double left=pop?0.5/sizes_[epoch]:0.5/base_sizes_[epoch];
        const double right=pop?0.5/right_sizes_[epoch]:0.5/base_right_sizes_[epoch];
        return (std::min(high,cut)-std::min(low,cut))*left
            +(std::max(high,cut)-std::max(low,cut))*right;
    }
    if(!linear_rates_)
        return (high-low)*(pop?0.5/sizes_[epoch]:0.5/base_sizes_[epoch]);
    return (high-low)*0.5*(rate(pop,epoch,low)+rate(pop,epoch,high));
}
void Ensemble::set_profile(const std::vector<double>& base,const std::vector<double>& base_right,
    const std::vector<double>& ghost,const std::vector<double>& ghost_right,
    const std::vector<double>& base_cuts,const std::vector<double>& ghost_cuts,
    bool linear,bool split_base,bool split_ghost) {
    if(base.size()+1!=base_boundaries_.size() || ghost.size()+1!=boundaries_.size()
        || base_right.size()!=base.size() || ghost_right.size()!=ghost.size())
        throw std::invalid_argument("wrong number of epoch sizes");
    for(double n:base) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("base epoch sizes must be positive and finite");
    for(double n:ghost) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("ghost epoch sizes must be positive and finite");
    for(double n:base_right) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("base right epoch sizes must be positive and finite");
    for(double n:ghost_right) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("ghost right epoch sizes must be positive and finite");
    if(linear && (split_base || split_ghost))
        throw std::invalid_argument("choose linear rates or movable cuts");
    if(base_cuts.size()!=(split_base?base.size():0)
        || ghost_cuts.size()!=(split_ghost?ghost.size():0))
        throw std::invalid_argument("movable cuts need one point per split epoch");
    if(split_base) {
        for(size_t e=0;e<base.size();++e)
            if(!std::isfinite(base_cuts[e]) || !(base_boundaries_[e]<base_cuts[e]
                && base_cuts[e]<base_boundaries_[e+1]))
                throw std::invalid_argument("base cut must lie strictly inside its epoch");
    }
    if(split_ghost) {
        for(size_t e=0;e<ghost.size();++e)
            if(!std::isfinite(ghost_cuts[e]) || !(boundaries_[e]<ghost_cuts[e]
                && ghost_cuts[e]<boundaries_[e+1]))
                throw std::invalid_argument("ghost cut must lie strictly inside its epoch");
    }
    const size_t old_capacity=base_sizes_.capacity()+sizes_.capacity()
        +base_right_sizes_.capacity()+right_sizes_.capacity()
        +base_cuts_.capacity()+cuts_.capacity();
    base_sizes_=base; sizes_=ghost;linear_rates_=linear;
    split_base_=split_base;split_ghost_=split_ghost;
    base_right_sizes_=base_right;right_sizes_=ghost_right;
    base_cuts_=base_cuts;cuts_=ghost_cuts;
    resident_bytes_+=(base_sizes_.capacity()+sizes_.capacity()
        +base_right_sizes_.capacity()+right_sizes_.capacity()
        +base_cuts_.capacity()+cuts_.capacity()-old_capacity)*sizeof(double);
    for_trees(trees_,threads_,[](const TreeCache& tree){return tree.prefix.size();},[&](size_t,size_t i) {
        auto& tree=trees_[i];
        if(tree.impossible) return;
        std::fill(tree.prefix.begin(),tree.prefix.end(),0);
        for(size_t j=0;j+1<tree.times.size();++j) {
            const double inv[2]={rate(0,tree.epochs[j][0],tree.times[j]),rate(1,tree.epochs[j][1],tree.times[j])};
            const double log_inv[2]={std::log(inv[0]),std::log(inv[1])};
            const double integrated[2]={
                rate_integral(0,tree.epochs[j][0],tree.times[j+1],tree.times[j]),
                rate_integral(1,tree.epochs[j][1],tree.times[j+1],tree.times[j])};
            for(size_t s=0;s<tree.scenarios;++s) {
                const auto& step=tree.steps[j*tree.scenarios+s];
                double v=tree.prefix[j*tree.scenarios+s];
                for(size_t pop=0;pop<2;++pop)
                    v+=step.coalescences[pop]*log_inv[pop]
                        -pairs(step.lineages[pop])*integrated[pop];
                tree.prefix[(j+1)*tree.scenarios+s]=v;
            }
        }
    });
}
size_t Ensemble::locate(TreeCache& tree,double t) {
    auto contains=[&](size_t i){return i<tree.times.size() && tree.times[i]<=t && (!i || t<tree.times[i-1]);};
    if(contains(tree.cursor)) return tree.cursor;
    if(tree.cursor && contains(tree.cursor-1)) return --tree.cursor;
    if(contains(tree.cursor+1)) return ++tree.cursor;
    return tree.cursor=std::lower_bound(tree.times.begin(),tree.times.end(),t,std::greater<double>())-tree.times.begin();
}
double Ensemble::partial(const TreeCache& tree,size_t s,size_t boundary,double t) const {
    if(tree.times[boundary]==t) return tree.prefix[boundary*tree.scenarios+s];
    const size_t j=boundary-1;
    const auto& step=tree.steps[j*tree.scenarios+s];
    const double inv[2]={rate(0,tree.epochs[j][0],tree.times[j]),rate(1,tree.epochs[j][1],tree.times[j])};
    const double integrated[2]={
        rate_integral(0,tree.epochs[j][0],t,tree.times[j]),
        rate_integral(1,tree.epochs[j][1],t,tree.times[j])};
    double v=tree.prefix[j*tree.scenarios+s];
    for(size_t pop=0;pop<2;++pop)
        v+=step.coalescences[pop]*std::log(inv[pop])
            -pairs(step.lineages[pop])*integrated[pop];
    return v;
}
double Ensemble::score(TreeCache& tree,size_t s,size_t boundary,double t,double p) const {
    const auto& k=tree.steps[(boundary-1)*tree.scenarios+s].lineages;
    if((p==0 && k[1]) || (p==1 && k[0])) return neginf;
    double v=partial(tree,s,boundary,t);
    if(p>0 && p<1) v+=k[0]*std::log1p(-p)+k[1]*std::log(p);
    return v;
}
void Ensemble::validate_query(double t,double p) const {
    if(!std::isfinite(t)||!std::isfinite(p)||t<0||t>=migration_limit()||p<0||p>1)
        throw std::invalid_argument("require 0 <= time < migration limit and 0 <= probability <= 1");
}
double Ensemble::evaluate(double t,double p) {
    validate_query(t,p);
    const bool interior=p>0 && p<1;
    const double log_base=interior?std::log1p(-p):0;
    const double log_ghost=interior?std::log(p):0;
    std::vector<double> values(trees_.size());
    for_trees(trees_,threads_,[](const TreeCache& tree){return tree.scenarios;},[&](size_t,size_t i) {
        auto& tree=trees_[i];
        if(tree.impossible) {values[i]=neginf;return;}
        size_t b=locate(tree,t);
        const size_t j=b-1;
        const bool exact=tree.times[b]==t;
        double inv[2]{},log_inv[2]{},integrated[2]{};
        if(!exact) {
            inv[0]=rate(0,tree.epochs[j][0],tree.times[j]);
            inv[1]=rate(1,tree.epochs[j][1],tree.times[j]);
            log_inv[0]=std::log(inv[0]);
            log_inv[1]=std::log(inv[1]);
            integrated[0]=rate_integral(0,tree.epochs[j][0],t,tree.times[j]);
            integrated[1]=rate_integral(1,tree.epochs[j][1],t,tree.times[j]);
        }
        double sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) {
            const auto& step=tree.steps[j*tree.scenarios+s];
            const auto& k=step.lineages;
            if((p==0 && k[1]) || (p==1 && k[0]))continue;
            double v=tree.prefix[(exact?b:j)*tree.scenarios+s];
            if(!exact) for(size_t pop=0;pop<2;++pop)
                v+=step.coalescences[pop]*log_inv[pop]
                    -pairs(k[pop])*integrated[pop];
            if(interior)v+=k[0]*log_base+k[1]*log_ghost;
            sum=add(sum,v);
        }
        values[i]=tree.root_log_likelihood+tree.prefix[(tree.times.size()-1)*tree.scenarios]-partial(tree,0,b,t)+sum;
    });
    double value=0;
    for(double v:values) {
        if(v==neginf)return neginf;
        value+=v;
    }
    return value;
}
std::vector<LikelihoodDerivatives> Ensemble::derivatives_by_tree(double t,double p) {
    validate_query(t,p);
    if(!(t>0) || !(p>0 && p<1))
        throw std::invalid_argument("derivatives require 0 < time < migration limit and 0 < probability < 1");
    std::vector<LikelihoodDerivatives> result;
    result.reserve(trees_.size());
    for(auto& tree:trees_) {
        if(tree.impossible)throw std::runtime_error("derivatives undefined for an impossible tree");
        const size_t b=locate(tree,t);
        if(tree.times[b]==t)
            throw std::invalid_argument("time derivative undefined at a cached event or epoch boundary");
        const size_t j=b-1;
        if((split_base_ && t==base_cuts_[tree.epochs[j][0]])
            || (split_ghost_ && t==cuts_[tree.epochs[j][1]]))
            throw std::invalid_argument("time derivative undefined at a movable cut");
        const double inv_base=rate(0,tree.epochs[j][0],t);
        const double inv_ghost=rate(1,tree.epochs[j][1],t);
        const auto& all=tree.steps[j*tree.scenarios].lineages;
        const double tail_rate=pairs(all[0]+all[1])*inv_base;
        std::vector<double> log_scores(tree.scenarios);
        double log_sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) {
            log_scores[s]=score(tree,s,b,t,p);
            log_sum=add(log_sum,log_scores[s]);
        }
        if(!std::isfinite(log_sum))
            throw std::runtime_error("derivatives undefined at zero tree likelihood");
        LikelihoodDerivatives entry;
        entry.log_likelihood=tree.root_log_likelihood
            +tree.prefix[(tree.times.size()-1)*tree.scenarios]-partial(tree,0,b,t)+log_sum;
        entry.time=-tail_rate;
        for(size_t s=0;s<tree.scenarios;++s) {
            const double weight=std::exp(log_scores[s]-log_sum);
            const auto& k=tree.steps[j*tree.scenarios+s].lineages;
            entry.time+=weight*(pairs(k[0])*inv_base+pairs(k[1])*inv_ghost);
            entry.probability+=weight*(k[1]/p-k[0]/(1-p));
        }
        result.push_back(entry);
    }
    return result;
}
LikelihoodDerivatives Ensemble::derivatives(double t,double p) {
    LikelihoodDerivatives result;
    for(const auto& tree:derivatives_by_tree(t,p)) {
        result.log_likelihood+=tree.log_likelihood;
        result.time+=tree.time;
        result.probability+=tree.probability;
    }
    return result;
}
Sufficient Ensemble::expectation(double t,double p) {
    validate_query(t,p);
    const size_t slots=std::min(trees_.size(),threads_?threads_:size_t(8));
    std::vector<Sufficient> partials(slots);
    for(auto& part:partials) {
        part.exposure.resize(sizes_.size(),0);part.coalescences.resize(sizes_.size(),0);
        part.base_exposure.resize(base_sizes_.size(),0);part.base_coalescences.resize(base_sizes_.size(),0);
        if(linear_rates_ || split_ghost_)part.ghost_rates.resize(sizes_.size());
        if(linear_rates_ || split_base_)part.base_rates.resize(base_sizes_.size());
    }
    for_trees(trees_,threads_,[](const TreeCache& tree){return tree.steps.size();},[&](size_t worker,size_t i) {
        auto& result=partials[worker];
        auto& tree=trees_[i];
        if(tree.impossible) throw std::runtime_error("cannot fit an impossible tree");
        size_t b=locate(tree,t);
        std::vector<double> weights(tree.scenarios);
        double sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) { weights[s]=score(tree,s,b,t,p); sum=add(sum,weights[s]); }
        if(!std::isfinite(sum)) throw std::runtime_error("zero likelihood at current parameters; use an interior probability start");
        result.log_likelihood+=tree.root_log_likelihood+tree.prefix[(tree.times.size()-1)*tree.scenarios]-partial(tree,0,b,t)+sum;
        for(auto& w:weights) w=std::exp(w-sum);
        auto add_rate=[&](size_t pop,size_t e,double low,double high,double pair_mass,double event_mass) {
            const bool split=pop?split_ghost_:split_base_;
            if(!linear_rates_ && !split)return;
            auto& evidence=pop?result.ghost_rates[e]:result.base_rates[e];
            const auto& bounds=pop?boundaries_:base_boundaries_;
            const double width=bounds[e+1]-bounds[e];
            const double zlo=std::clamp((low-bounds[e])/width,0.,1.);
            const double zhi=std::clamp((high-bounds[e])/width,0.,1.);
            const double right=pair_mass*(high-low)*(zlo+zhi)*0.5;
            evidence.right_exposure+=right;
            evidence.left_exposure+=pair_mass*(high-low)-right;
            if(event_mass>0) {
                if(linear_rates_)evidence.events.emplace_back(zhi,event_mass);
                if(split)evidence.cut_events.emplace_back(high,event_mass);
            }
            if(split && pair_mass>0 && high>low)
                evidence.spans.push_back({low,high,pair_mass});
        };
        for(size_t j=0;j<b;++j) {
            const double dt=tree.times[j]-std::max(t,tree.times[j+1]);
            const auto e=tree.epochs[j];
            std::array<double,2> pair_mass{},event_mass{};
            for(size_t s=0;s<tree.scenarios;++s) {
                const auto& step=tree.steps[j*tree.scenarios+s];
                result.exposure[e[1]]+=weights[s]*pairs(step.lineages[1])*dt;
                result.coalescences[e[1]]+=weights[s]*step.coalescences[1];
                result.base_exposure[e[0]]+=weights[s]*pairs(step.lineages[0])*dt;
                result.base_coalescences[e[0]]+=weights[s]*step.coalescences[0];
                if(linear_rates_ || split_base_ || split_ghost_)for(size_t pop=0;pop<2;++pop) {
                    pair_mass[pop]+=weights[s]*pairs(step.lineages[pop]);
                    event_mass[pop]+=weights[s]*step.coalescences[pop];
                }
            }
            if(linear_rates_ || split_base_ || split_ghost_)for(size_t pop=0;pop<2;++pop)
                add_rate(pop,e[pop],std::max(t,tree.times[j+1]),tree.times[j],
                    pair_mass[pop],event_mass[pop]);
        }
        // Below migration all colourings have joined the base population.
        // These intervals are the scenario-independent tail used by evaluate().
        if(tree.times[b]<t) {
            const size_t j=b-1, e=tree.epochs[j][0];
            const auto& step=tree.steps[j*tree.scenarios];
            result.base_exposure[e]+=pairs(step.lineages[0])*(t-tree.times[b]);
            add_rate(0,e,tree.times[b],t,pairs(step.lineages[0]),0);
        }
        for(size_t j=b;j+1<tree.times.size();++j) {
            const size_t e=tree.epochs[j][0];
            const auto& step=tree.steps[j*tree.scenarios];
            result.base_exposure[e]+=pairs(step.lineages[0])*(tree.times[j]-tree.times[j+1]);
            result.base_coalescences[e]+=step.coalescences[0];
            add_rate(0,e,tree.times[j+1],tree.times[j],pairs(step.lineages[0]),step.coalescences[0]);
        }
    });
    Sufficient result=std::move(partials[0]);
    for(size_t w=1;w<partials.size();++w) {
        const auto& part=partials[w];
        result.log_likelihood+=part.log_likelihood;
        for(size_t e=0;e<result.exposure.size();++e) {
            result.exposure[e]+=part.exposure[e];
            result.coalescences[e]+=part.coalescences[e];
        }
        for(size_t e=0;e<result.base_exposure.size();++e) {
            result.base_exposure[e]+=part.base_exposure[e];
            result.base_coalescences[e]+=part.base_coalescences[e];
        }
        if(linear_rates_ || split_base_ || split_ghost_) {
            auto merge=[&](std::vector<RateEvidence>& to,const std::vector<RateEvidence>& from) {
                for(size_t e=0;e<to.size();++e) {
                    to[e].left_exposure+=from[e].left_exposure;
                    to[e].right_exposure+=from[e].right_exposure;
                    to[e].events.insert(to[e].events.end(),from[e].events.begin(),from[e].events.end());
                    to[e].cut_events.insert(to[e].cut_events.end(),from[e].cut_events.begin(),from[e].cut_events.end());
                    to[e].spans.insert(to[e].spans.end(),from[e].spans.begin(),from[e].spans.end());
                }
            };
            merge(result.base_rates,part.base_rates);
            merge(result.ghost_rates,part.ghost_rates);
        }
    }
    return result;
}
std::array<std::vector<KernelEstimate>,2> Ensemble::kernel_profile(
    double t,double p,double bandwidth,size_t points) {
    validate_query(t,p);
    if(!std::isfinite(bandwidth) || !(bandwidth>0) || points<2 || points>1000000)
        throw std::invalid_argument("kernel profile needs positive bandwidth and 2..1000000 points");
    // Posterior event masses and pair-exposure rates are accumulated once.
    std::array<std::map<double,double>,2> events,changes;
    auto interval=[&](size_t pop,double low,double high,double rate) {
        if(high>low && rate>0) {
            changes[pop][low]+=rate;
            changes[pop][high]-=rate;
        }
    };
    for(auto& tree:trees_) {
        if(tree.impossible)throw std::runtime_error("cannot smooth an impossible tree");
        const size_t b=locate(tree,t);
        std::vector<double> weights(tree.scenarios);
        double total=neginf;
        for(size_t s=0;s<tree.scenarios;++s) {
            weights[s]=score(tree,s,b,t,p);
            total=add(total,weights[s]);
        }
        if(!std::isfinite(total))
            throw std::runtime_error("zero likelihood while smoothing population sizes");
        for(auto& weight:weights)weight=std::exp(weight-total);
        for(size_t j=0;j<b;++j) {
            const double high=tree.times[j],low=std::max(t,tree.times[j+1]);
            std::array<double,2> pair_rate{},coal{};
            for(size_t s=0;s<tree.scenarios;++s) {
                const auto& step=tree.steps[j*tree.scenarios+s];
                for(size_t q=0;q<2;++q) {
                    pair_rate[q]+=weights[s]*pairs(step.lineages[q]);
                    coal[q]+=weights[s]*step.coalescences[q];
                }
            }
            for(size_t q=0;q<2;++q) {
                interval(q,low,high,pair_rate[q]);
                if(coal[q]>0)events[q][high]+=coal[q];
            }
        }
        // Below migration the base lineage count does not depend on colouring.
        if(tree.times[b]<t) {
            const auto& step=tree.steps[(b-1)*tree.scenarios];
            interval(0,tree.times[b],t,pairs(step.lineages[0]));
        }
        for(size_t j=b;j+1<tree.times.size();++j) {
            const auto& step=tree.steps[j*tree.scenarios];
            interval(0,tree.times[j+1],tree.times[j],pairs(step.lineages[0]));
            if(step.coalescences[0])events[0][tree.times[j]]+=step.coalescences[0];
        }
    }
    struct Segment {double low,high,rate;};
    std::array<std::vector<Segment>,2> segments;
    for(size_t q=0;q<2;++q) {
        double rate=0;
        for(auto it=changes[q].begin();it!=changes[q].end();++it) {
            rate+=it->second;
            const auto next=std::next(it);
            if(next!=changes[q].end() && next->first>it->first && rate>0)
                segments[q].push_back({it->first,next->first,rate});
        }
    }
    std::vector<double> grid;
    grid.reserve(points+1);
    for(size_t i=0;i<points;++i)
        grid.push_back(split1_*double(i)/double(points-1));
    grid.push_back(t);
    std::sort(grid.begin(),grid.end());
    grid.erase(std::unique(grid.begin(),grid.end()),grid.end());
    std::array<std::vector<KernelEstimate>,2> result;
    const double scale=bandwidth*std::sqrt(std::numbers::pi/2.),
        divisor=std::sqrt(2.)*bandwidth;
    for(size_t q=0;q<2;++q) {
        result[q].reserve(grid.size());
        for(double x:grid) {
            KernelEstimate row;
            row.time=x;
            row.size=std::numeric_limits<double>::quiet_NaN();
            if(q==1 && x<=t) {result[q].push_back(row);continue;}
            for(const auto& [event_time,mass]:events[q]) {
                const double z=(event_time-x)/bandwidth;
                row.coalescences+=mass*std::exp(-0.5*z*z);
            }
            for(const auto& segment:segments[q]) {
                const double kernel_integral=scale*(
                    std::erf((segment.high-x)/divisor)
                    -std::erf((segment.low-x)/divisor));
                row.exposure+=segment.rate*std::max(0.,kernel_integral);
            }
            if(row.coalescences>0 && row.exposure>0)
                row.size=row.exposure/(2*row.coalescences);
            result[q].push_back(row);
        }
    }
    return result;
}
std::array<std::vector<std::pair<double,double>>,2>
Ensemble::coalescence_masses(double t,double p) {
    validate_query(t,p);
    std::map<double,std::array<double,2>> mass;
    for(auto& tree:trees_) {
        if(tree.impossible)throw std::runtime_error("cannot build balanced epochs from an impossible tree");
        const size_t b=locate(tree,t);
        std::vector<double> weights(tree.scenarios);
        double sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) {
            weights[s]=score(tree,s,b,t,p);
            sum=add(sum,weights[s]);
        }
        if(!std::isfinite(sum))
            throw std::runtime_error("zero likelihood while choosing balanced epoch boundaries");
        for(auto& w:weights)w=std::exp(w-sum);
        for(size_t j=0;j+1<tree.times.size();++j) {
            const double event_time=tree.times[j];
            if(event_time>t) {
                for(size_t s=0;s<tree.scenarios;++s) {
                    const auto& c=tree.steps[j*tree.scenarios+s].coalescences;
                    for(size_t q=0;q<2;++q)mass[event_time][q]+=weights[s]*c[q];
                }
            } else {
                const auto& c=tree.steps[j*tree.scenarios].coalescences;
                mass[event_time][0]+=c[0];
            }
        }
    }
    std::array<std::vector<std::pair<double,double>>,2> result;
    for(const auto& [time,pair]:mass)
        for(size_t q=0;q<2;++q)if(pair[q]>0)result[q].emplace_back(time,pair[q]);
    return result;
}
size_t Ensemble::bytes() const { return resident_bytes_; }
size_t Ensemble::scenario_count() const { size_t n=0;for(const auto& t:trees_)n+=t.scenarios;return n; }
}
