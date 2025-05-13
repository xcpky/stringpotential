#include "autofree.h"
#include "lse.h"
#include "wavefunction.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

void testwf() {
  WaveFunction *wf __attribute__((cleanup(auto_wffree))) = WFnew(0, 20, 64);
  for (size_t i = 0; i < N_MAX; i += 1) {
    for (size_t j = 0; j < N_MAX; j += 1) {
      const size_t idx = i * wf->c_solution->tda + j;
      printf("%.4e\t", wf->c_solution->data[idx]);
    }
    puts("");
  }
}

void printmat(matrix *m) {
  __auto_type row = m->size1;
  __auto_type col = m->size2;
  for (uint i = 0; i < row; i += 1) {
    for (uint j = 0; j < col; j += 1) {
      __auto_type val = matrix_get(m, i, j);
      printf("(%.2e) + Im(%.2e) ", creal(val), cimag(val));
    }
    puts("");
  }
}

void testlse() {
  const size_t Ngauss = 2;
  LSE *lse __attribute__((cleanup(auto_lsefree))) = lse_malloc(Ngauss, 4, 1e-9);
  if (!lse) {
    fprintf(stderr, "Failed to create LSE solver\n");
    exit(1);
  }
  lse_compute(lse, 0.2, 1);
  // lse_gmat(lse);
  // lse_vmat(lse);
  // puts("G matrix");
  // gsl_matrix_complex_fprintf(stdout, lse->G, "%.2e");
  puts("V matrix");
  printmat(lse->V);
  puts("T matrix");
  printmat(lse->T);
  // gsl_matrix_complex_fprintf(stdout, lse->T, "%.2e");
  // puts("I - VG matrix");
  // gsl_matrix_complex_fprintf(stdout, lse->iIVG, "%.2e");
  // int result = lse_run(lse);
  // if (result != 0) {
  //   fprintf(stderr, "LSE solver failed with error code %d\n", result);
  //   lse_deinit(lse);
  //   return 1;
  // }
  //
  // printf("LSE solver ran successfully\n");
}

void testonshell() {
  const size_t Ngauss = 2;
  LSE *lse __attribute__((cleanup(auto_lsefree))) = lse_malloc(Ngauss, 4, 1e-9);
  if (!lse) {
    fprintf(stderr, "Failed to create LSE solver\n");
    exit(1);
  }
  lse_refresh(lse, -0.578802970, 1);
  __auto_type q = lse->x0[0];
  printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
  q = lse->x0[1];
  printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
  lse_vmat(lse);
  __auto_type onshellV = matrix_get(lse->V, Ngauss, Ngauss);
  printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
  lse_refresh(lse, -0.578802965, 1);
  lse_vmat(lse);
  q = lse->x0[0];
  printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
  q = lse->x0[1];
  printf("(%.12e) + Im(%.12e)\n", creal(q), cimag(q));
  onshellV = matrix_get(lse->V, Ngauss, Ngauss);
  printf("(%.12e) + Im(%.12e)\n", creal(onshellV), cimag(onshellV));
}

void unitest() {
  double re;
  double im;
  while (scanf("%lf %lf", &re, &im) != 0) {
    __auto_type val = O_00(0.3, re + im * I, 1, m_B);
    printf("E(%.2e): (%.2e) + Im(%.2e)\n", re, creal(val), cimag(val));
  }
}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    puts("what do you want?");
    return EXIT_SUCCESS;
  }
  if (strcmp(argv[1], "lse") == 0) {
    testlse();
  } else if (strcmp(argv[1], "wf") == 0) {
    testwf();
  } else if (strcmp(argv[1], "unit") == 0) {
    unitest();
  } else if (strcmp(argv[1], "onshell") == 0) {
    testonshell();
  }
  return EXIT_SUCCESS;
}
