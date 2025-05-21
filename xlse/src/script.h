#ifndef SCRIPT_H
#define SCRIPT_H

#define NTHREADS 16

#include "complex.h"
#include "lse.h"
#include <stddef.h>
#include <threads.h>

typedef struct {
  void *ose00; // Using void* to avoid complex.h in C++ code
  void *ose01;
  void *ose10;
  void *ose11;
} onshellElements;

// Declare the function with void* to avoid complex.h dependency in C++
double complex *onshellT(double *E, size_t len, size_t pNgauss);
void Free(void* ptr);

void ose_free(onshellElements ose);
typedef struct {
  LSE *lse;
  size_t start;
  size_t len;
  double *E;
  onshellElements *res;
  int64_t rs;
  size_t id;
} argstruct;

int oT(void *arg);

#endif // !SCRIPT_H
