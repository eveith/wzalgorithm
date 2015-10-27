#ifndef REVOLTEST_H_
#define REVOLTEST_H_


#include <QObject>

#include "Testcase.h"


class REvolTest: public Testcase
{
    Q_OBJECT

public:
    explicit REvolTest(QObject *parent = 0);
    static qreal peaks(const qreal &x, const qreal &y);
    static qreal ackley(const qreal &x, const qreal &y);

private slots:
    void testPeaks();
    void testAckley();
};


#endif
