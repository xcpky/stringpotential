#ifndef SCRIPT_H
#define SCRIPT_H

#define NTHREADS 32

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
double complex *onshellT(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex *onshellG(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex *onshellV(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex *Det(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex *Both(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double *Evec(size_t pNgauss, double Lambda, double epsilon);
double complex *Poles(double *Er, size_t rlen, double *Ei, size_t ilen, size_t pNgauss, double Lambda, double epsilon);
void Free(void* ptr);

void ose_free(onshellElements ose);
typedef struct {
  size_t pNgauss;
  double Lambda;
  double epsilon;
  size_t start;
  size_t len;
  double *E;
  void *res;
  RS rs;
  size_t id;
} argstruct;

struct polestruct {
  size_t pNgauss;
  double Lambda;
  double epsilon;
  size_t start;
  size_t len;
  double *Er;
  size_t rlen;
  double *Ei;
  size_t ilen;
  void *res;
  RS rs;
  size_t id;
};

int oT(void *arg);
int oG(void *arg);
int oV(void *arg);
int det(void *arg);
int both(void *arg);
int poles(void *arg);

#endif // !SCRIPT_H
