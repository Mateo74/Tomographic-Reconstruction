#ifndef METRICS_H
#define METRICS_H

#include <vector>

struct Metrics
{
    double reconstruction_time;
    double cond_number;
    std::vector<double> psnr;
    unsigned num_eigen_found;
};

#endif