#include <cmath>

#include <gtest/gtest.h>

#include "peaks.h"
#include "ackley.h"
#include "crossin.h"

#include "ParticleSwarmOptimization.h"
#include "ParticleSwarmOptimizationTest.h"


using std::pow;
using std::exp;
using std::cos;
using std::sin;

using wzalgorithm::ParticleSwarmOptimization;
using Particle = wzalgorithm::ParticleSwarmOptimization::Particle;


TEST(ParticleSwarmOptimizationTest, testPeaks)
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-10.0).upperBoundary(10.0).maxIterations(5000);

    bool success = false;
    pso.run(2, [&success](Particle& particle) {
        double r = peaks(
                particle.currentPosition[0],
                particle.currentPosition[1]);
        particle.currentFitness = r;
        success |= (r <= -6.0);
        return success;
    });

    ASSERT_TRUE(success);
}


TEST(ParticleSwarmOptimizationTest, testAckley)
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-32.0).upperBoundary(32.0).maxIterations(50000);

    bool success = false;
    auto p = pso.run(2, [&success](Particle& particle) {
        double r = ackley(particle.currentPosition);
        success |= (r + 1.0 < 1.0000000001);
        particle.currentFitness = r;
        return success;
    });

    ASSERT_NEAR(ackley(wzalgorithm::vector_t{ 0.0, 0.0 }), 0.0, 1e-6);
    ASSERT_TRUE(success);
    ASSERT_TRUE(ackley(p.bestParticle.bestPosition) + 1.0 < 1.0000000001);
}


TEST(ParticleSwarmOptimizationTest, testCrossInTray)
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-10.0).upperBoundary(10.0).maxIterations(5000);

    auto succeeds = [](Particle& p) {
        p.currentFitness = ::crossInTray(
                p.currentPosition[0],
                p.currentPosition[1]);
        return p.currentFitness < -2.0;
    };

    auto result = pso.run(2, succeeds);

    ASSERT_TRUE(succeeds(result.bestParticle));
}
