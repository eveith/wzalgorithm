#include <cmath>

#include <QtTest>
#include <QObject>
#include <QVector>

#include "Testcase.h"
#include "TestcaseManager.h"

#include "REvol.h"
#include "REvolTest.h"


using std::pow;
using std::exp;
using std::cos;
using std::sin;

using Winzent::Algorithm::REvol;
using Winzent::Algorithm::detail::Individual;


qreal REvolTest::peaks(const qreal &x, const qreal &y)
{
    return 3.0 * pow(1.0 - x, 2) * exp(-pow(x, 2) - pow(y + 1.0, 2))
            - 10.0 * (x / 5.0 - pow(x, 3) - pow(y, 5))
            * exp(- pow(x, 2) - pow(y, 2))
            - 1.0 / 3.0 * exp(- pow(x + 1.0, 2) - pow(y, 2));
}


qreal REvolTest::ackley(const qreal &x, const qreal &y)
{
    static qreal a = 20.0;
    static qreal b = 0.2;
    static qreal c = 2.0 * M_PI;
    return -a * exp(-b * sqrt(0.5 * (pow(x, 2.0) + pow(y, 2.0))))
            - exp(0.5 * (cos(c*x) + cos(c*y)))
            + a
            + exp(1.0);
}


REvolTest::REvolTest(QObject *parent): Testcase(parent)
{
}


void REvolTest::testPeaks()
{
    REvol revol;
    revol
            .ebmax(1.0)
            .gradientWeight(1.8)
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
        qreal r = peaks(i.parameters[0], i.parameters[1]);

        // Round to see dynamic probability spread in action:

        r *= pow(10, 2);
        r = ceil(r);
        r /= pow(10, 2);

        i.restrictions.push_front(r);
        i.restrictions.resize(1);

        success |= (i.restrictions[0] <= -6.0);
        return success;
    });

    QVERIFY(success);
    QVERIFY(result.bestIndividual.restrictions[0] < -6.0);
}


void REvolTest::testAckley()
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
    auto result = revol.run(i, [&success](Individual &i) {
        qreal r = ackley(i.parameters[0], i.parameters[1]);
        i.restrictions[0] = r;

        success |= (i.restrictions[0] + 1.0 < 1.0000000001);
        return success;
    });

    QCOMPARE(1.0 + ackley(0.0, 0.0), 1.0);
    QVERIFY(success);
    QVERIFY(result.bestIndividual.restrictions[0] + 1.0 < 1.000000001);
}



void REvolTest::testCompareIndividuals()
{
    Individual i1, i2;

    QCOMPARE(0, i1.compare(i2));

    i2.timeToLive = 1;
    QCOMPARE(0, i1.compare(i2));

    i1.age();
    QCOMPARE(-1, i1.compare(i2));

    i1.timeToLive = 1;
    i2.timeToLive = 1;

    i1.restrictions << 1.0;
    QCOMPARE(0, i1.compare(i2));
    QVERIFY(! i1.isBetterThan(i2));

    i2.restrictions << 1.0;
    QVERIFY(!i1.isBetterThan(i2));
    QCOMPARE(0, i2.compare(i1));

    i1.restrictions[0] = 5.0;
    QVERIFY(! i1.isBetterThan(i2));
    QCOMPARE(1, i2.compare(i1));

    i1.restrictions[0] = i2.restrictions[0];
    i1.restrictions << 1.0 << 2.0;
    i2.restrictions << 1.0 << 1.0;
    QVERIFY(i2.isBetterThan(i1));
}


void REvolTest::testSortPopulation()
{
    Individual *i1 = new Individual();
    Individual *i2 = new Individual();
    Individual *i3 = new Individual();

    i1->timeToLive = 10;
    i1->restrictions.push_back(0.25);

    i2->timeToLive = 10;
    i2->restrictions.push_back(0.5);

    i3->timeToLive = 2;
    i3->restrictions.push_back(0.5);

    QVERIFY(i1->isBetterThan(*i2));

    Winzent::Algorithm::REvol::Population population;
    population.push_back(i2);
    population.push_back(i1);
    population.push_back(i3);

    QCOMPARE(&(population.front()), i2);
    population.sort();
    QCOMPARE(&(population.at(0)), i1);
    QCOMPARE(&(population.at(1)), i2);
    QCOMPARE(&(population.at(2)), i3);
}


REGISTER_TESTCASE(REvolTest);
