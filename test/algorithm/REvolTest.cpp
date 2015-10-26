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

using Winzent::Algorithm::REvol;
using Winzent::Algorithm::Individual;


qreal REvolTest::peaks(const qreal &x, const qreal &y)
{
    return 3.0 * pow(1.0 - x, 2) * exp(-pow(x, 2) - pow(y + 1.0, 2))
            - 10.0 * (x / 5.0 - pow(x, 3) - pow(y, 5))
            * exp(- pow(x, 2) - pow(y, 2))
            - 1.0 / 3.0 * exp(- pow(x + 1.0, 2) - pow(y, 2));
}


REvolTest::REvolTest(QObject *parent): Testcase(parent)
{
}


void REvolTest::testPeaks()
{
    REvol revol;
    revol
            .gradientWeight(1.8)
            .populationSize(200)
            .eliteSize(50)
            .successWeight(1.5)
            .measurementEpochs(50)
            .startTTL(20)
            .maxEpochs(500)
            .maxNoSuccessEpochs(500);

    Individual i;
    i.parameters = { 0.0, 0.0 };
    i.scatter = { 8.0, 8.0 };

    bool success = false;
    i = revol.run(i, [&success](Individual &i) {
        i.restrictions.push_front(peaks(i.parameters[0], i.parameters[1]));
        i.restrictions.resize(1);
        success |= (i.restrictions[0] <= -6.0);
        return false; // Keep on training until the end.
    });

    QVERIFY(success);
    QVERIFY(i.restrictions[0] < -6.0);
}


REGISTER_TESTCASE(REvolTest);
