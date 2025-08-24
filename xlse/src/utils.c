#include "utils.h"
#include <stdio.h>
char *formatC(double complex x)
{
    static char buf[64];
    sprintf(buf, "%f%+fim", creal(x), cimag(x));
    return buf;
}
