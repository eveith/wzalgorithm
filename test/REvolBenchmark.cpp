#include <cmath>
#include <vector>
#include <algorithm>

#include <benchmark/benchmark.h>

#include "ackley.h"
#include "crossin.h"

#include "REvol.h"
#include "REvolBenchmark.h"


using std::cos;
using std::pow;
using std::sqrt;
using std::accumulate;

using wzalgorithm::REvol;
using Individual = wzalgorithm::REvol::Individual;


static void REvolAckleyBenchmark(benchmark::State& state)
{
    while (state.KeepRunning()) {
        REvol revol;
        revol
                .ebmax(5.0)
                .gradientWeight(1.0)
                .populationSize(40)
                .eliteSize(4)
                .successWeight(1.0)
                .measurementEpochs(100)
                .startTTL(120)
                .maxEpochs(5000)
                .maxNoSuccessEpochs(5000);

        Individual i;
        i.parameters = {
            5.0, 5.0, 7.0, 5.0, 5.0, 2.0, 6.42, 2.67, 1.998, 19.3 };
        i.scatter = {
            32., 32., 32., 32., 32., 32., 32., 32., 32., 32. };

        bool success = false;
        revol.run(i, [&success](Individual& i) {
            double r = ackley(i.parameters);
            i.restrictions[0] = r;

            success |= (i.restrictions[0] + 1.0 < 1.0000000001);
            return success;
        });
    }
}
BENCHMARK(REvolAckleyBenchmark);


static void REvolCrossInTrayBenchmark(benchmark::State& state)
{
    auto succeeds = [](Individual& i) {
        i.restrictions[0] = ::crossInTray(i.parameters[0], i.parameters[1]);
        return i.restrictions[0] < 2.062;
    };

    while (state.KeepRunning()) {
        REvol revol;
        revol.maxEpochs(5000);

        Individual origin;
        origin.scatter = {0.1, 0.1};
        origin.parameters = {-10.0, 10.0};

        auto result = revol.run(origin, succeeds);
    }
}
BENCHMARK(REvolCrossInTrayBenchmark);
