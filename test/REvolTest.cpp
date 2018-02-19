#include <cmath>

#include <gtest/gtest.h>

#include "peaks.h"
#include "ackley.h"
#include "crossin.h"

#include "REvol.h"
#include "REvolTest.h"


using std::pow;
using std::exp;
using std::cos;
using std::sin;

using wzalgorithm::REvol;
using Individual = wzalgorithm::REvol::Individual;


TEST(REvolTest, testPeaks)
{
    REvol revol;
    revol
            .ebmax(10.0)
            .gradientWeight(1.0)
            .populationSize(30)
            .eliteSize(3)
            .successWeight(1.0)
            .measurementEpochs(500)
            .startTTL(300)
            .maxEpochs(5000)
            .maxNoSuccessEpochs(4000);

    Individual i;
    i.parameters = { 0.0, 0.0 };
    i.scatter = { 8.0, 8.0 };

    bool success = false;
    auto result = revol.run(i, [&success](Individual &i) {
        double r = peaks(i.parameters[0], i.parameters[1]);

        // Round to see dynamic probability spread in action:

        r *= pow(10, 2);
        r = ceil(r);
        r /= pow(10, 2);

        i.restrictions.resize(1);
        i.restrictions[0] = r;

        success |= (i.restrictions[0] <= -6.0);
        return success;
    });

    ASSERT_TRUE(success);
    ASSERT_TRUE(result.error < -6.0);
}


TEST(REvolTest, testAckley)
{
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
    i.parameters = { 5.0, 5.0 };
    i.scatter = { 10.0, 10.0 };

    bool success = false;
    auto result = revol.run(i, [&success](Individual& i) {
        double r = ::ackley(i.parameters);
        i.restrictions[0] = r;

        success |= (i.restrictions[0] + 1.0 < 1.0000000001);
        return success;
    });

    ASSERT_NEAR(ackley(wzalgorithm::vector_t{ 0.0, 0.0 }), 0.0, 1e-6);
    ASSERT_TRUE(success);
    ASSERT_TRUE(result.error + 1.0 < 1.000000001);
}


TEST(REvolTest, testCrossInTray)
{
    REvol revol;
    revol.maxEpochs(5000);

    Individual origin;
    origin.scatter = {1, 1};
    origin.parameters = {-10.0, 10.0};

    auto succeeds = [](Individual& i) {
        i.restrictions[0] = ::crossInTray(i.parameters[0], i.parameters[1]);
        return i.restrictions[0] < -2.062;
    };

    auto result = revol.run(origin, succeeds);

    ASSERT_TRUE(::crossInTray(result.parameters[0], result.parameters[1])
            < -2.0);
}


TEST(REvolTest, testCompareIndividuals)
{
    Individual i1, i2;

    ASSERT_TRUE(! i1.isBetterThan(i2));

    i2.timeToLive = 1;
    ASSERT_TRUE(! i1.isBetterThan(i2));

    i1.age();
    ASSERT_TRUE(i2.isBetterThan(i1));

    i1.timeToLive = 1;
    i2.timeToLive = 1;

    i1.restrictions.push_back(1.0);
    ASSERT_FALSE(i1.isBetterThan(i2));

    i2.restrictions.push_back(1.0);
    ASSERT_FALSE(i1.isBetterThan(i2));

    i1.restrictions[0] = 5.0;
    ASSERT_FALSE(i1.isBetterThan(i2));

    i1.restrictions[0] = i2.restrictions[0];
    i1.restrictions.push_back(1.0);
    i1.restrictions.push_back(2.0);
    i2.restrictions.push_back(1.0);
    i2.restrictions.push_back(1.0);
    ASSERT_TRUE(i2.isBetterThan(i1));
}


TEST(REvolTest, testSortPopulation)
{
    auto i1 = new Individual(),
            i2 = new Individual(),
            i3 = new Individual();

    i1->timeToLive = 10;
    i1->restrictions.push_back(0.25);

    i2->timeToLive = 10;
    i2->restrictions.push_back(0.5);

    i3->timeToLive = 2;
    i3->restrictions.push_back(0.5);

    ASSERT_TRUE(i1->isBetterThan(*i2));

    wzalgorithm::REvol::Population population;
    population.push_back(i2);
    population.push_back(i1);
    population.push_back(i3);

    ASSERT_EQ(population.front(), *i2);
    population.sort();
    ASSERT_EQ(population.at(0), *i1);
    ASSERT_EQ(population.at(1), *i2);
    ASSERT_EQ(population.at(2), *i3);
}


TEST(REvolTest, testGenerateOrigin)
{
    size_t d = 5;
    auto o = REvol().generateOrigin(d);

    ASSERT_EQ(d, o.scatter.size());
    ASSERT_EQ(o.parameters.size(), o.scatter.size());

    for (size_t i = 0; i != d; ++i) {
        ASSERT_TRUE(o.scatter[i] != 0);
        ASSERT_TRUE(o.parameters[i] != 0);
    }
}
