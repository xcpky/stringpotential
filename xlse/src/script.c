#include "script.h"
#include "lse.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <stdlib.h>

double complex *onshellT(double *E, size_t len, size_t pNgauss) {
  // puts("in C");
  // for (size_t i = 0; i < len; i += 1) {
  //   printf("%f,  ", E[i]);
  // }
  // puts("");
  double complex *res =
      (double complex *)malloc(sizeof(double complex) * len * 4);
  onshellElements onshellBuf = {
      .ose00 = res,
      .ose01 = res + len,
      .ose10 = res + 2 * len,
      .ose11 = res + 3 * len,
  };
  thrd_t tid[NTHREADS];
  argstruct args[NTHREADS];
  size_t ntasks = len / NTHREADS;
  size_t residue = len % NTHREADS;
  for (size_t i = 0; i < NTHREADS; i += 1) {
    args[i].pNgauss = pNgauss;
    args[i].Lambda = 4;
    args[i].epsilon = 1e-7;
    args[i].res = &onshellBuf;
    args[i].rs = 1;
    args[i].E = E;
    args[i].id = i;
    if (i < residue) {
      args[i].len = ntasks + 1;
      args[i].start = (ntasks + 1) * i;
    } else {
      args[i].len = ntasks;
      args[i].start = ntasks * i + residue;
    }
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_create(&tid[i], oT, &args[i]);
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_join(tid[i], NULL);
  }
  return res;
}

double complex *Det(double *E, size_t len, size_t pNgauss) {
  double complex *res = (double complex *)malloc(sizeof(double complex) * len);
  thrd_t tid[NTHREADS];
  argstruct args[NTHREADS];
  size_t ntasks = len / NTHREADS;
  size_t residue = len % NTHREADS;
  for (size_t i = 0; i < NTHREADS; i += 1) {
    args[i].pNgauss = pNgauss;
    args[i].Lambda = 4;
    args[i].epsilon = 1e-7;
    args[i].res = res;
    args[i].rs = 1;
    args[i].E = E;
    args[i].id = i;
    if (i < residue) {
      args[i].len = ntasks + 1;
      args[i].start = (ntasks + 1) * i;
    } else {
      args[i].len = ntasks;
      args[i].start = ntasks * i + residue;
    }
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_create(&tid[i], det, &args[i]);
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_join(tid[i], NULL);
  }
  return res;
}

void ose_free(onshellElements ose) { free(ose.ose00); }

int oT(void *arg) {
  argstruct foo = *(argstruct *)arg;
  LSE *lse = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
  size_t ngauss = foo.pNgauss;
  onshellElements *res = (onshellElements *)foo.res;
  double complex *ose00 = (double complex *)res->ose00;
  double complex *ose01 = (double complex *)res->ose01;
  double complex *ose10 = (double complex *)res->ose10;
  double complex *ose11 = (double complex *)res->ose11;
  // printf("start: %lu, len: %lu\n", foo.start, foo.len);
  for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
    lse_compute(lse, foo.E[i], foo.rs);
    // printf("%f\n", foo.E[i]);
    size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
    ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // ose00[i] = gsl_matrix_complex_get(lse->T, ngauss, ngauss);
    idx += ngauss + 1;
    ose01[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // ose01[i] = matrix_get(lse->T, ngauss, 2 * ngauss + 1);
    idx = (2 * ngauss + 1) * 2 * (ngauss + 1) + ngauss;
    ose10[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // ose10[i] = matrix_get(lse->T, 2 * ngauss + 1, ngauss);
    idx += ngauss + 1;
    ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // ose11[i] = matrix_get(lse->T, 2 * ngauss + 1, 2 * ngauss + 1);
  }
  return 0;
}

int det(void *arg) {
  argstruct foo = *(argstruct *)arg;
  LSE *lse = lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
  double complex *res = (double complex *)foo.res;
  for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
    res[i] = lse_detImVG(lse, foo.E[i]);
  }
  return 0;
}
void Free(void *ptr) { free(ptr); }
