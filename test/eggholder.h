#ifndef EGGHOLDER_H
#define EGGHOLDER_H


#include <cmath>


using std::abs;
using std::sin;
using std::sqrt;


static double eggholder(double x1, double x2)
{
    return -(x2 + 47.) * sin(sqrt(abs(x2 + x1/2. + 47)))
            - x1 * sin(sqrt(abs(x1 - (x2 + 47.))));
}


#endif // EGGHOLDER_H
