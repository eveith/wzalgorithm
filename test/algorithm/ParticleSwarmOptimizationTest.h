#ifndef PARTICLESWARMOPTIMIZATIONTEST_H_
#define PARTICLESWARMOPTIMIZATIONTEST_H_


#include <QObject>

#include "Testcase.h"


class ParticleSwarmOptimizationTest: public Testcase
{
    Q_OBJECT

public:
    explicit ParticleSwarmOptimizationTest(QObject *parent = 0);
    static qreal peaks(const qreal &x, const qreal &y);
    static qreal ackley(const qreal &x, const qreal &y);

private slots:
    void testPeaks();
    void testAckley();
};


#endif
