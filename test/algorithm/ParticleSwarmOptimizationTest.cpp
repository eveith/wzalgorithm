#include <cmath>

#include <QtTest>
#include <QObject>
#include <QVector>

#include "Testcase.h"
#include "TestcaseManager.h"

#include "ParticleSwarmOptimization.h"
#include "ParticleSwarmOptimizationTest.h"


using std::pow;
using std::exp;
using std::cos;
using std::sin;

using Winzent::Algorithm::ParticleSwarmOptimization;


qreal ParticleSwarmOptimizationTest::peaks(const qreal &x, const qreal &y)
{
    return 3.0 * pow(1.0 - x, 2) * exp(-pow(x, 2) - pow(y + 1.0, 2))
            - 10.0 * (x / 5.0 - pow(x, 3) - pow(y, 5))
            * exp(- pow(x, 2) - pow(y, 2))
            - 1.0 / 3.0 * exp(- pow(x + 1.0, 2) - pow(y, 2));
}


qreal ParticleSwarmOptimizationTest::ackley(const qreal &x, const qreal &y)
{
    static qreal a = 20.0;
    static qreal b = 0.2;
    static qreal c = 2.0 * M_PI;
    return -a * exp(-b * sqrt(0.5 * (pow(x, 2.0) + pow(y, 2.0))))
            - exp(0.5 * (cos(c*x) + cos(c*y)))
            + a
            + exp(1.0);
}


ParticleSwarmOptimizationTest::ParticleSwarmOptimizationTest(QObject *parent):
        Testcase(parent)
{
}


void ParticleSwarmOptimizationTest::testPeaks()
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-10.0).upperBoundary(10.0);

    bool success = false;
    pso.run(2, [&success](const QVector<qreal> &parameters) {
        qreal r = peaks(parameters[0], parameters[1]);
        success |= (r <= -6.0);
        return success;
    });

    QVERIFY(success);
}


void ParticleSwarmOptimizationTest::testAckley()
{
    ParticleSwarmOptimization pso;
    pso.lowerBoundary(-10.0).upperBoundary(10.0);

    bool success = false;
    auto p = pso.run(2, [&success](const QVector<qreal> &parameters) {
        qreal r = ackley(parameters[0], parameters[1]);
        success |= (r + 1.0 < 1.0000000001);
        return success;
    });

    QCOMPARE(1.0 + ackley(0.0, 0.0), 1.0);
    QVERIFY(success);
    QVERIFY(ackley(p[0], p[1]) + 1.0 < 1.000000001);
}

REGISTER_TESTCASE(ParticleSwarmOptimizationTest);
