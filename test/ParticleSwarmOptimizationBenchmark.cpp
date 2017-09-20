#include <cmath>
#include <vector>
#include <algorithm>

#include <benchmark/benchmark.h>

#include "ackley.h"
#include "crossin.h"
#include "ParticleSwarmOptimization.h"
#include "ParticleSwarmOptimizationBenchmark.h"


using std::cos;
using std::pow;
using std::sqrt;
using std::accumulate;


using wzalgorithm::ParticleSwarmOptimization;
using Particle = wzalgorithm::ParticleSwarmOptimization::Particle;



static void PsoAckleyBenchmark(benchmark::State& state)
{
    while (state.KeepRunning()) {
        ParticleSwarmOptimization pso;
        pso.lowerBoundary(-32.0).upperBoundary(32.0).maxIterations(50000);

        bool success = false;
        auto p = pso.run(10, [&success](Particle& particle) {
            double r = ackley(particle.currentPosition);
            success |= (r + 1.0 < 1.0000000001);
            particle.currentFitness = r;
            return success;
        });
    }
}
BENCHMARK(PsoAckleyBenchmark);



static void PsoCrossInTrayBenchmark(benchmark::State& state)
{
    auto succeeds = [](Particle& p) {
        p.currentFitness = ::crossInTray(
                p.currentPosition[0],
                p.currentPosition[1]);
        return p.currentFitness < 2.062;
    };

    while (state.KeepRunning()) {
        ParticleSwarmOptimization pso;
        pso.lowerBoundary(-10.0).upperBoundary(10.0).maxIterations(5000);
        auto result = pso.run(2, succeeds);
    }
}
BENCHMARK(PsoCrossInTrayBenchmark);
