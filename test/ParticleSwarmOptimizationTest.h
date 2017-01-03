#ifndef PARTICLESWARMOPTIMIZATIONTEST_H_
#define PARTICLESWARMOPTIMIZATIONTEST_H_


#include <gtest/gtest.h>


class ParticleSwarmOptimizationTest: public ::testing::Test
{
public:
    static double peaks(double x, double y);
    static double ackley(double x, double y);
};


#endif
