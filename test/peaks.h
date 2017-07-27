#ifndef PEAKS_H
#define PEAKS_H


#include <cmath>


using std::pow;
using std::exp;


static double peaks(double x, double y)
{
    return 3.0 * pow(1.0 - x, 2) * exp(-pow(x, 2) - pow(y + 1.0, 2))
            - 10.0 * (x / 5.0 - pow(x, 3) - pow(y, 5))
            * exp(- pow(x, 2) - pow(y, 2))
            - 1.0 / 3.0 * exp(- pow(x + 1.0, 2) - pow(y, 2));
}


#endif // PEAKS_H
