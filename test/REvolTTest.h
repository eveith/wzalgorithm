#ifndef REVOLTTEST_H_
#define REVOLTTEST_H_


#include <gtest/gtest.h>


class REvolTTest: public ::testing::Test
{
public:
    static double peaks(double x, double y);
    static double ackley(double x, double y);
};


#endif
