#ifndef REVOLTBENCHMARK_H
#define REVOLTBENCHMARK_H


#include <vector>
#include <benchmark/benchmark.h>


class REvolTBenchmark : public benchmark::Fixture
{
public:
    double ackley(std::vector<double> const& x) const;
};


#endif // REVOLBENCHMARK_H
