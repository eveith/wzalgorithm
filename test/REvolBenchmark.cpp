#include <cmath>
#include <vector>
#include <algorithm>

#include <benchmark/benchmark.h>

#include "REvol.h"
#include "REvolBenchmark.h"


using std::cos;
using std::pow;
using std::sqrt;
using std::accumulate;

using Winzent::Algorithm::REvol;
using Winzent::Algorithm::detail::Individual;


double REvolBenchmark::ackley(std::vector<double> const& x) const
{
    static double a = 20.0;
    static double b = 0.2;
    static double c = 2.0 * M_PI;
    auto n = x.size();

    return -a * exp(-b * sqrt(1./n * accumulate(
                x.begin(),
                x.end(),
                0.0,
                [](double sum, double x) { return sum + pow(x, 2); })))
            - exp(1./n * accumulate(
                x.begin(),
                x.end(),
                0.0,
                [](double sum, double x) { return sum + cos(c * x); }))
            + a
            + exp(1.0);
}


BENCHMARK_F(REvolBenchmark, ackleyBenchmark)(benchmark::State& state)
{
    while (state.KeepRunning()) {
        REvol revol;
        revol
                .ebmax(5.0)
                .gradientWeight(1.0)
                .populationSize(33)
                .eliteSize(3)
                .successWeight(1.2)
                .measurementEpochs(500)
                .startTTL(100)
                .maxEpochs(5000)
                .maxNoSuccessEpochs(5000);

        Individual i;
        i.parameters = {
            5.0, 5.0, 7.0, 5.0, 5.0, 2.0, 6.42, 2.67, 1.998, 19.3 };
        i.scatter = {
            10.0, 10.0, 10., 10., 10., 10., 10., 10., 10., 10. };

        bool success = false;
        auto result = revol.run(i, [this, &success](Individual& i) {
            double r = this->ackley(i.parameters);
            i.restrictions[0] = r;

            success |= (i.restrictions[0] + 1.0 < 1.0000000001);
            return success;
        });
    }
}
