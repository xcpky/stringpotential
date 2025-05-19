#include "script.h"
#include "interface.h"
#include <stdlib.h>

onshellElements onshellT(void *E, size_t len, size_t pNgauss) {
  double complex *buf =
      (double complex *)malloc(sizeof(double complex) * len * 4);
  onshellElements res = {
      .ose00 = buf,
      .ose01 = buf + len,
      .ose10 = buf + 2 * len,
      .ose11 = buf + 3 * len,
  };
  thrd_t tid[NTHREADS];
  argstruct args[NTHREADS];
  size_t ntasks = len / NTHREADS;
  size_t residue = len % NTHREADS;
  for (size_t i = 0; i < NTHREADS; i += 1) {
    args[i].lse = lse_malloc(pNgauss, 4, 1e-6);
    args[i].res = &res;
    args[i].rs = 1;
    args[i].E = E;
    if (i < residue) {
      args[i].len = ntasks + 1;
      args[i].start = (ntasks + 1) * i;
    } else {
      args[i].len = ntasks;
      args[i].start = ntasks * i + residue;
    }

    thrd_create(&tid[i], oT, &args[i]);
  }
  for (size_t i = 0; i < NTHREADS; i += 1) {
    thrd_join(tid[i], NULL);
  }
  return res;
}

void ose_free(onshellElements ose) {
  free(ose.ose00);
}

int oT(void *arg) {
  argstruct rox = *(argstruct *)arg;
  size_t ngauss = rox.lse->Ngauss;
  double complex *ose00 = (double complex *)rox.res->ose00;
  double complex *ose01 = (double complex *)rox.res->ose01;
  double complex *ose10 = (double complex *)rox.res->ose10;
  double complex *ose11 = (double complex *)rox.res->ose11;
  for (size_t i = rox.start; i < rox.len; i += 1) {
    lse_compute(rox.lse, rox.E[i], rox.rs);
    size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
    ose00[i] =
        rox.lse->T->data[2 * idx] + rox.lse->T->data[2 * idx + 1] * I;
    idx += ngauss + 1;
    ose01[i] =
        rox.lse->T->data[2 * idx] + rox.lse->T->data[2 * idx + 1] * I;
    idx = (2 * ngauss + 1) * (ngauss + 1) + ngauss;
    ose10[i] =
        rox.lse->T->data[2 * idx] + rox.lse->T->data[2 * idx + 1] * I;
    idx += ngauss + 1;
    ose11[i] =
        rox.lse->T->data[2 * idx] + rox.lse->T->data[2 * idx + 1] * I;
  }
  return 0;
}
