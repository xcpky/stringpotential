#ifndef UTILS_H
#define UTILS_H

#include <complex.h>
#include <math.h>

char *formatC(double complex);
// Complex square root function
static inline double complex xsqrt(double complex x)
{
    // return csqrt(x);
	return cimag(x) >= 0 ? csqrt(x + 0I):-csqrt(x-0I);
    return cimag(x) < 0 && creal(x) < 0 ? -csqrt(x - 0I) : csqrt(x + 0I);
}
static inline double complex xlog(double complex x)
{
    // return cimag(x) >= 0 ? clog(x + 0I) : clog(x + 0I) + 2 * M_PI * I;
	return clog(x);
    auto res = clog(x);
	// return res;
    // return clog(cabs(x));
    if (cimag(res) < 0 && fabs(cimag(x)) < 1e-7 && creal(x) < 0) {
	res += 2 * M_PI * I;
    }
    return res;
}
#endif // UTILS_H
