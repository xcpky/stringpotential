#include "script.h"
#include "autofree.h"
#include "lse.h"
#include <gsl/gsl_matrix_complex_double.h>
#include <stdlib.h>
#include <string.h>

double complex *onshellT(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon) {
  // puts("in C");
  // for (size_t i = 0; i < len; i += 1) {
  //   printf("%f,  ", E[i]);
  // }
  // puts("");
  double complex *res = malloc(sizeof(double complex) * len * 4);
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
    args[i].Lambda = Lambda;
    args[i].epsilon = epsilon;
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

double complex *Det(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon) {
  double complex *res = malloc(sizeof(double complex) * len);
  thrd_t tid[NTHREADS];
  argstruct args[NTHREADS];
  size_t ntasks = len / NTHREADS;
  size_t residue = len % NTHREADS;
  for (size_t i = 0; i < NTHREADS; i += 1) {
    args[i].pNgauss = pNgauss;
    args[i].Lambda = Lambda;
    args[i].epsilon = epsilon;
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

double complex *Both(double *E, size_t len, size_t pNgauss, double Lambda, double epsilon) {
  double complex *res = malloc(sizeof(double complex) * len * 5);
  for (size_t i = 0; i < len; i += 1) {
    res[i] = i;
  }
  thrd_t tid[NTHREADS];
  argstruct args[NTHREADS];
  size_t ntasks = len / NTHREADS;
  size_t residue = len % NTHREADS;
  for (size_t i = 0; i < NTHREADS; i += 1) {
    args[i].pNgauss = pNgauss;
    args[i].Lambda = Lambda;
    args[i].epsilon = epsilon;
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
    thrd_create(&tid[i], both, &args[i]);
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_join(tid[i], NULL);
  }
  return res;
}

double *Evec(size_t pNgauss, double Lambda, double epsilon) {
  LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
  double *E = malloc(sizeof(double) * N_MAX);
  memcpy(E, lse->E_vec, sizeof(double[N_MAX]));
  return E;
}

void ose_free(onshellElements ose) { free(ose.ose00); }

int oT(void *arg) {
  argstruct foo = *(argstruct *)arg;
  LSE *lse [[gnu::cleanup(lsefree)]] =
      lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
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
#ifdef TPO
    lse_X(lse);
    lse_XtX(lse);
    ose00[i] = lse->onshellT[0][0];
    ose01[i] = lse->onshellT[0][1];
    ose10[i] = lse->onshellT[1][0];
    ose11[i] = lse->onshellT[1][1];
#else
    auto T = (double complex(*)[2 * ngauss + 2]) lse->TOME->data;
    ose00[i] = T[ngauss][ngauss];
    ose01[i] = T[ngauss][2 * ngauss + 1];
    ose10[i] = T[2 * ngauss + 1][ngauss];
    ose11[i] = T[2 * ngauss + 1][2 * ngauss + 1];
    // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
    // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // idx += ngauss + 1;
    // ose01[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // idx = (2 * ngauss + 1) * 2 * (ngauss + 1) + ngauss;
    // ose10[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
    // idx += ngauss + 1;
    // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx + 1] * I;
#endif
  }
  return 0;
}

int det(void *arg) {
  argstruct foo = *(argstruct *)arg;
  LSE *lse [[gnu::cleanup(lsefree)]] =
      lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
  double complex *res = (double complex *)foo.res;
  for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
    res[i] = lse_detImVG(lse, foo.E[i]);
    // lse_compute(lse, foo.E[i], foo.rs);
    // res[i] = lse->det;
  }
  return 0;
}

int both(void *arg) {
  argstruct foo = *(argstruct *)arg;
  LSE *lse [[gnu::cleanup(lsefree)]] =
      lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
  auto buf = (double complex *)foo.res;
  size_t pNgauss = foo.pNgauss;
  // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
  // for (size_t i = 0; i < 11; i += 1) {
  //   printf("%.2f ", creal(buf[i]));
  // }
  // puts("");
  for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
    lse_compute(lse, foo.E[i], foo.rs);
    auto res = buf + 5 * i;
    auto T = (double complex(*)[2 * pNgauss + 2]) lse->TOME->data;
    // printf("thread %lu, i %lu, det[i] %.1f\n", foo.id, i, creal(det[i]));
    // puts("");
    res[0] = lse->det;
    res[1] = T[pNgauss][pNgauss];
    res[2] = T[pNgauss][2 * pNgauss + 1];
    res[3] = T[2 * pNgauss + 1][pNgauss];
    res[4] = T[2 * pNgauss + 1][2 * pNgauss + 1];
  }
  return 0;
}

void Free(void *ptr) { free(ptr); }
