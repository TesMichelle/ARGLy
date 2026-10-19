#include "epoch_model.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

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
        + tree.epochs.capacity()*sizeof(size_t)
        + tree.steps.capacity()*sizeof(Step)
        + tree.prefix.capacity()*sizeof(double);
}
}

Ensemble::Ensemble(const tsk_treeseq_t& ts, const std::vector<double>& p,
    const std::vector<int>& pops, int admixed, int outgroup, size_t epoch_count,
    size_t limit, double memory_limit_mb,
    const std::function<void(size_t,size_t,size_t)>& progress,
    const std::vector<double>& epoch_boundaries)
    : base_size_(p.at(6)), split2_(p.at(1)), split1_(p.at(2)) {
    if (p.size()!=9 || !epoch_count || epoch_count>100000 || admixed==outgroup
        || !(base_size_>0) || !(p[7]>0) || !(split1_>0))
        throw std::invalid_argument("invalid epoch model parameters");
    for (double x:p) if (!std::isfinite(x)) throw std::invalid_argument("parameters must be finite");
    if(epoch_boundaries.empty()) {
        boundaries_.resize(epoch_count+1);
        for (size_t e=0;e<=epoch_count;++e) boundaries_[e]=split1_*double(e)/epoch_count;
    } else {
        if(epoch_boundaries.size()!=epoch_count+1 || epoch_boundaries.front()!=0
            || epoch_boundaries.back()!=split1_)
            throw std::invalid_argument("epoch boundaries must have K+1 values from 0 to split1");
        boundaries_=epoch_boundaries;
    }
    for (size_t e=1;e<boundaries_.size();++e)
        if (!std::isfinite(boundaries_[e]) || !(boundaries_[e]>boundaries_[e-1]))
            throw std::invalid_argument("epoch boundaries must be finite and strictly increasing");
    resident_bytes_=sizeof(*this)+boundaries_.capacity()*sizeof(double);
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
        helper.sort_by_time(tree); helper.compute_lower_nodes_split_1(tree); helper.compute_fixed_nodes(tree);
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
        if(has_outgroup_)cache.times.push_back(split2_);
        for(auto id:nodes) if(ts.tables->nodes.time[id]>0 && ts.tables->nodes.time[id]<=split1_) cache.times.push_back(ts.tables->nodes.time[id]);
        std::sort(cache.times.begin(),cache.times.end(),std::greater<double>());
        cache.times.erase(std::unique(cache.times.begin(),cache.times.end()),cache.times.end());
        cache.epochs.resize(cache.times.size(),0);
        for(size_t j=0;j+1<cache.times.size();++j) {
            const double mid=(cache.times[j]+cache.times[j+1])/2;
            cache.epochs[j]=std::upper_bound(boundaries_.begin(),boundaries_.end(),mid)-boundaries_.begin()-1;
        }
        const size_t records=product(cache.scenarios,cache.times.size());
        const size_t allocation=product(records,sizeof(Step)+sizeof(double));
        if(memory_limit_mb>0 && (static_cast<long double>(bytes())+allocation
            +cache.times.capacity()*sizeof(double)+cache.epochs.capacity()*sizeof(size_t)
            +2*epoch_count*sizeof(double))/(1024*1024)>memory_limit_mb)
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
    set_sizes(std::vector<double>(sizes.size(),base_size_),sizes);
}
void Ensemble::set_sizes(const std::vector<double>& base, const std::vector<double>& ghost) {
    if(base.size()+1!=boundaries_.size() || ghost.size()!=base.size())
        throw std::invalid_argument("wrong number of epoch sizes");
    for(double n:base) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("base epoch sizes must be positive and finite");
    for(double n:ghost) if(!std::isfinite(n)||n<=0) throw std::invalid_argument("ghost epoch sizes must be positive and finite");
    const size_t old_base=base_sizes_.capacity(), old_ghost=sizes_.capacity();
    base_sizes_=base; sizes_=ghost;
    resident_bytes_+=(base_sizes_.capacity()-old_base+sizes_.capacity()-old_ghost)*sizeof(double);
    for(auto& tree:trees_) {
        if(tree.impossible) continue;
        std::fill(tree.prefix.begin(),tree.prefix.end(),0);
        for(size_t j=0;j+1<tree.times.size();++j) {
            const double inv[2]={0.5/base_sizes_[tree.epochs[j]],0.5/sizes_[tree.epochs[j]]};
            const double log_inv[2]={std::log(inv[0]),std::log(inv[1])};
            for(size_t s=0;s<tree.scenarios;++s) {
                const auto& step=tree.steps[j*tree.scenarios+s];
                double v=tree.prefix[j*tree.scenarios+s];
                for(size_t pop=0;pop<2;++pop)
                    v+=step.coalescences[pop]*log_inv[pop]-pairs(step.lineages[pop])*(tree.times[j]-tree.times[j+1])*inv[pop];
                tree.prefix[(j+1)*tree.scenarios+s]=v;
            }
        }
    }
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
    const double inv[2]={0.5/base_sizes_[tree.epochs[j]],0.5/sizes_[tree.epochs[j]]};
    double v=tree.prefix[j*tree.scenarios+s];
    for(size_t pop=0;pop<2;++pop)
        v+=step.coalescences[pop]*std::log(inv[pop])-pairs(step.lineages[pop])*(tree.times[j]-t)*inv[pop];
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
    double value=0;
    for(auto& tree:trees_) {
        if(tree.impossible) return neginf;
        size_t b=locate(tree,t);
        double sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) sum=add(sum,score(tree,s,b,t,p));
        value+=tree.root_log_likelihood+tree.prefix[(tree.times.size()-1)*tree.scenarios]-partial(tree,0,b,t)+sum;
    }
    return value;
}
Sufficient Ensemble::expectation(double t,double p) {
    validate_query(t,p);
    Sufficient result;
    result.exposure.resize(sizes_.size(),0); result.coalescences.resize(sizes_.size(),0);
    result.base_exposure.resize(sizes_.size(),0); result.base_coalescences.resize(sizes_.size(),0);
    for(auto& tree:trees_) {
        if(tree.impossible) throw std::runtime_error("cannot fit an impossible tree");
        size_t b=locate(tree,t);
        std::vector<double> weights(tree.scenarios);
        double sum=neginf;
        for(size_t s=0;s<tree.scenarios;++s) { weights[s]=score(tree,s,b,t,p); sum=add(sum,weights[s]); }
        if(!std::isfinite(sum)) throw std::runtime_error("zero likelihood at current parameters; use an interior probability start");
        result.log_likelihood+=tree.root_log_likelihood+tree.prefix[(tree.times.size()-1)*tree.scenarios]-partial(tree,0,b,t)+sum;
        for(auto& w:weights) w=std::exp(w-sum);
        for(size_t j=0;j<b;++j) {
            const double dt=tree.times[j]-std::max(t,tree.times[j+1]);
            const size_t e=tree.epochs[j];
            for(size_t s=0;s<tree.scenarios;++s) {
                const auto& step=tree.steps[j*tree.scenarios+s];
                result.exposure[e]+=weights[s]*pairs(step.lineages[1])*dt;
                result.coalescences[e]+=weights[s]*step.coalescences[1];
                result.base_exposure[e]+=weights[s]*pairs(step.lineages[0])*dt;
                result.base_coalescences[e]+=weights[s]*step.coalescences[0];
            }
        }
        // Below migration all colourings have joined the base population.
        // These intervals are the scenario-independent tail used by evaluate().
        if(tree.times[b]<t) {
            const size_t j=b-1, e=tree.epochs[j];
            const auto& step=tree.steps[j*tree.scenarios];
            result.base_exposure[e]+=pairs(step.lineages[0])*(t-tree.times[b]);
        }
        for(size_t j=b;j+1<tree.times.size();++j) {
            const size_t e=tree.epochs[j];
            const auto& step=tree.steps[j*tree.scenarios];
            result.base_exposure[e]+=pairs(step.lineages[0])*(tree.times[j]-tree.times[j+1]);
            result.base_coalescences[e]+=step.coalescences[0];
        }
    }
    return result;
}
size_t Ensemble::bytes() const { return resident_bytes_; }
size_t Ensemble::scenario_count() const { size_t n=0;for(const auto& t:trees_)n+=t.scenarios;return n; }
}
