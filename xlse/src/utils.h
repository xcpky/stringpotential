#ifndef UTILS_H
#define UTILS_H

#include <complex.h>

char *formatC(double complex);
// Complex square root function
inline double complex xsqrt(double complex x) { return cimag(x) >= 0 ? csqrt(x + 0I) : -csqrt(x - 0I); }

#endif // UTILS_H
