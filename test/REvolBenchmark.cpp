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
    const size_t dimensions = 10;
    auto succeeds = [](Individual& i) {
        i.restrictions[0] = ackley(i.parameters);
        return (i.restrictions[0] + 1.0 < 1.0000000001);
    };

    while (state.KeepRunning()) {
        REvol revol;
        revol
                .maxEpochs(5000)
                .maxNoSuccessEpochs(5000)
                .run(revol.generateOrigin(dimensions), succeeds);
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
