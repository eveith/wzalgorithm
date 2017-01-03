#include <cmath>

#include <gtest/gtest.h>

#include "ParticleSwarmOptimization.h"
#include "ParticleSwarmOptimizationTest.h"


using std::pow;
using std::exp;
using std::cos;
using std::sin;

using Winzent::Algorithm::detail::Particle;
using Winzent::Algorithm::ParticleSwarmOptimization;


double ParticleSwarmOptimizationTest::peaks(double x, double y)
{
    return 3.0 * pow(1.0 - x, 2) * exp(-pow(x, 2) - pow(y + 1.0, 2))
            - 10.0 * (x / 5.0 - pow(x, 3) - pow(y, 5))
            * exp(- pow(x, 2) - pow(y, 2))
            - 1.0 / 3.0 * exp(- pow(x + 1.0, 2) - pow(y, 2));
}


double ParticleSwarmOptimizationTest::ackley(double x, double y)
{
    static double a = 20.0;
    static double b = 0.2;
    static double c = 2.0 * M_PI;
    return -a * exp(-b * sqrt(0.5 * (pow(x, 2.0) + pow(y, 2.0))))
            - exp(0.5 * (cos(c*x) + cos(c*y)))
            + a
            + exp(1.0);
}


TEST_F(ParticleSwarmOptimizationTest, testPeaks)
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-10.0).upperBoundary(10.0).maxIterations(5000);

    bool success = false;
    pso.run(2, [&success](Particle &particle) {
        double r = peaks(
                particle.currentPosition[0],
                particle.currentPosition[1]);
        particle.currentFitness = r;
        success |= (r <= -6.0);
        return success;
    });

    ASSERT_TRUE(success);
}


TEST_F(ParticleSwarmOptimizationTest, testAckley)
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-32.0).upperBoundary(32.0).maxIterations(50000);

    bool success = false;
    auto p = pso.run(2, [&success](Particle& particle) {
        double r = ackley(
                particle.currentPosition[0],
                particle.currentPosition[1]);
        success |= (r + 1.0 < 1.0000000001);
        particle.currentFitness = r;
        return success;
    });

    ASSERT_NEAR(ackley(0.0, 0.0), 0.0, 1e-6);
    ASSERT_TRUE(success);
    ASSERT_TRUE(ackley(
            p.bestParticle.bestPosition[0],
            p.bestParticle.bestPosition[1]) + 1.0 < 1.0000000001);
}
