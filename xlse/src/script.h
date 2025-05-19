#ifndef SCRIPT_H
#define SCRIPT_H

#define NTHREADS 32

#include "complex.h"
#include "lse.h"
#include <stddef.h>
#include <threads.h>
#include "interface.h"

typedef struct {
  LSE *lse;
  size_t start;
  size_t len;
  double complex *E;
  onshellElements *res;
  int64_t rs;
} argstruct;

int oT(void *arg);


#endif // !SCRIPT_H
