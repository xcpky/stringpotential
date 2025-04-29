#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>
#define HAVE_INLINE
#include "lse.h"
#include "wavefunction.h"
#include <gsl/gsl_eigen.h>

void testwf() {
  WaveFunction *wf = WFnew(1, 20, RNGAUSS);
  double complex **psi_n_mat =
      (double complex **)malloc(sizeof(double complex *) * N_MAX);
  for (uint64_t i = 0; i < PNGAUSS; i += 1) {
    psi_n_mat[i] = (double complex *)malloc(sizeof(double complex) * PNGAUSS);
  }
  double *pi = (double *)malloc(sizeof(double) * PNGAUSS);
  double *wi = (double *)malloc(sizeof(double) * PNGAUSS);
  gsl_integration_glfixed_table *table =
      gsl_integration_glfixed_table_alloc(PNGAUSS);
  for (uint64_t i = 0; i < PNGAUSS; i += 1) {
    gsl_integration_glfixed_point(0, 4, i, &pi[i], &wi[i], table);
  }
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    psi_n_ft_batch(wf, pi, psi_n_mat[i], PNGAUSS, i + 1);
  }
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    printf("%.3e + Im %.3e\n", creal(psi_n_mat[i][0]), cimag(psi_n_mat[i][0]));
  }
  // printf("Eigenvalues:\n");
  // for (int i = 0; i < N_MAX; i++) {
  //   printf("%d: %f\n", i, gsl_vector_get(wf->E_solution, i));
  // }
  // printf("c_solutions:\n");
  // for (uint64_t i = 0; i < N_MAX; i+=1) {
  //   printf("%.3e\n", gsl_matrix_get(wf->c_solution, i, 0));
  // }
  // double complex psi = psi_n(wf, 1, 1, 0);
  // printf("%f + Im%f\n", GSL_REAL(psi), GSL_IMAG(psi));
  // double r = 0;
  // while (1) {
  //   scanf("%lf", &r);
  //   if (r < 0)
  //     break;
  //   psi = psi_n(wf, r, 1, 0);
  //   printf("%f + Im%f\n", GSL_REAL(psi), GSL_IMAG(psi));
  // }
  free(pi);
  free(wi);
  gsl_integration_glfixed_table_free(table);
  WFfree(wf);
}

void testlse() {
  LSE *lse = lse_malloc(2, 4, 1e-7);
  if (!lse) {
    fprintf(stderr, "Failed to create LSE solver\n");
    exit(1);
  }
  lse_compute(lse, 0.2);
  // lse_gmat(lse);
  // lse_vmat(lse);
  uint64_t dim = 2 * (lse->Ngauss + 1);
  puts("G matrix");
  for (uint64_t i = 0; i < dim; i += 1) {
    for (uint64_t j = 0; j < dim; j += 1) {
      uint64_t pos = i * dim + j;
      printf("%.2e + Im%.2e   ", lse->G->data[2 * pos],
             lse->G->data[2 * pos + 1]);
    }
    printf("\n");
  }
  puts("V matrix");
  for (uint64_t i = 0; i < dim; i += 1) {
    for (uint64_t j = 0; j < dim; j += 1) {
      uint64_t pos = i * dim + j;
      printf("%.2e + Im%.2e   ", lse->V->data[2 * pos],
             lse->V->data[2 * pos + 1]);
    }
    printf("\n");
  }
  puts("T matrix");
  for (uint64_t i = 0; i < dim; i += 1) {
    for (uint64_t j = 0; j < dim; j += 1) {
      uint64_t pos = i * dim + j;
      printf("%.2e + Im%.2e   ", lse->T->data[2 * pos],
             lse->T->data[2 * pos + 1]);
    }
    printf("\n");
  }
  // int result = lse_run(lse);
  // if (result != 0) {
  //   fprintf(stderr, "LSE solver failed with error code %d\n", result);
  //   lse_deinit(lse);
  //   return 1;
  // }
  //
  // printf("LSE solver ran successfully\n");

  // Clean up
  lse_free(lse);
}

int main() {
  testlse();
  return 0;
}
