#pragma once
#include "epoch_model.h"
#include <string>
namespace epochfit {
struct FitOptions {
    size_t iterations=30, max_evaluations=100;
    double min_size=100, max_size=1e7, tolerance=1e-5;
    bool fit_ghost=true;
    std::array<double,2> time_bounds{}, probability_bounds{1e-6,1-1e-6};
    std::function<void(size_t,double,const Parameters&,double)> iteration;
    std::function<void(size_t,size_t,double)> optimizer;
};
struct FitResult {
    Parameters parameters;
    double log_likelihood=0;
    size_t iterations=0;
    std::string status;
};
FitResult fit(Ensemble&, Parameters, const FitOptions&);
}
