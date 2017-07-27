#ifndef ACKLEY_H
#define ACKLEY_H


#include <cmath>
#include <numeric>

#include "config.h"


template <typename T>
double ackley(T const& x)
{
    static const double a = 20.0;
    static const double b = 0.2;
    static const double c = 2.0 * M_PI;
    auto n = x.size();

    return -a * exp(-b * sqrt(1./n * std::accumulate(
                std::begin(x),
                std::end(x),
                0.0,
                [](double sum, double x) { return sum + std::pow(x, 2); })))
            - exp(1./n * std::accumulate(
                std::begin(x),
                std::end(x),
                0.0,
                [](double sum, double x) { return sum + std::cos(c * x); }))
            + a
            + exp(1.0);
}


#endif // ACKLEY_H
