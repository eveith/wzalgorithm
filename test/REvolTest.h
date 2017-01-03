#ifndef REVOLTEST_H_
#define REVOLTEST_H_


#include <gtest/gtest.h>


class REvolTest: public ::testing::Test
{
public:
    static double peaks(double x, double y);
    static double ackley(double x, double y);
};


#endif
