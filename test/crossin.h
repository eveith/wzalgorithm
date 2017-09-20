#ifndef CROSSIN_H
#define CROSSIN_H


#include <cmath>


using std::abs;
using std::pow;
using std::sin;
using std::exp;
using std::sqrt;


static double crossInTray(double x, double y)
{
    return -0.0001 * pow(
            (abs(
                sin(x)
                * sin(y)
                * exp(
                    abs(
                        100 - sqrt(pow(x, 2) * pow(y, 2))/ M_PI)))
            + 1),
            0.1);
}


#endif // CROSSIN_H
