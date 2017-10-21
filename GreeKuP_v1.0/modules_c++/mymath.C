#ifndef MYMATH_CPP
#define MYMATH_CPP
#include <math.h>

double deltagaussian(double e, double deltae)
{
        return (1 / (deltae * sqrt(2 * M_PI))) * exp ( -(e * e) / (2 * deltae * deltae));
}


#endif
