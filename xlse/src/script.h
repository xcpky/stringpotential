#ifndef SCRIPT_H
#define SCRIPT_H

#include <stddef.h>
#include <threads.h>

#include "complex.h"
#include "lse.h"

typedef struct {
    void* ose00;  // Using void* to avoid complex.h in C++ code
    void* ose01;
    void* ose10;
    void* ose11;
} onshellElements;

// Declare the function with void* to avoid complex.h dependency in C++
double complex* onshell(double* E, size_t len, size_t pNgauss, double Lambda, double epsilon);
double complex* onshellT(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* onshellT_single(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* onshellG(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* onshellV(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* onshellTV(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* traceG(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex* Det(double* E, size_t len, double C[4], uint64_t rs, size_t pNgauss, double Lambda, double epsilon);
double complex* Det_single(double* E, size_t len, double C[4], uint64_t rs, size_t pNgauss, double Lambda, double epsilon);
double complex* Both(double* E, size_t len, double C[4], size_t pNgauss, double Lambda, double epsilon);
double complex V(double E, double complex p, double complex pprime);
double complex* V3d(double E, size_t pNgauss, double Lambda, double epsilon);
double complex* getV(double E, size_t pNgauss, double Lambda, double epsilon);
double complex* getIntegrand();
struct OME* ome_malloc();
void* WF();
double complex getwf(void* wf, double complex p, uint64_t n);
double complex getwfx(void* wf, double r, uint64_t n);
double complex getphi(double r, uint64_t n);
void WFend(void*);
void Evec(double*, double);
void Ntower(size_t*);
double complex* Poles(double* Er, size_t rlen, double* Ei, size_t ilen, double C[4], size_t pNgauss, double Lambda,
                      double epsilon);
double* Fit(double* C, size_t clen, size_t pNgauss, double Lambda, double epsilon);
int thrdfit(void*);
double* Fitsing(double* C, size_t clen, size_t pNgauss, double Lambda, double epsilon);
int thrdfitsing(void*);
double* Cost(double* C, size_t len, double complex resonance, size_t pNgauss, double Lambda, double epsilon);
void Free(void* ptr);

void ose_free(onshellElements ose);
typedef struct {
    size_t pNgauss;
    double Lambda;
    double epsilon;
    size_t start;
    size_t len;
    double* E;
    double C[4];
    void* res;
    RS rs;
    size_t id;
    size_t size;
} arg1d;

struct polestruct {
    size_t pNgauss;
    double Lambda;
    double epsilon;
    size_t start;
    size_t len;
    double* Er;
    size_t rlen;
    double* Ei;
    size_t ilen;
    double C[4];
    void* res;
    RS rs;
    size_t id;
};

int oT(void* arg);
int oTsing(void* arg);
int oG(void* arg);
int oV(void* arg);
int oTV(void* arg);
int trG(void* arg);
int det(void* arg);
int detsing(void* arg);
int both(void* arg);
int poles(void* arg);
int cst(void* arg);

#endif  // !SCRIPT_H
