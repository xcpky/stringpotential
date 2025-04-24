#define HAVE_INLINE
#include "wavefunction.h"
#include <gsl/gsl_eigen.h>

void test() {
  int l = 0;            // Angular momentum quantum number
  double V_0_fit = 0.4; // Potential parameter
  double V_0_l = -0.08858455651140933;
  double sigma_l = 0.0199179973550142;
  double alpha_l = 0.476814032326273;
  double r_0_l = 15;
  double V_0 = V_0_l / a_l;
  double sigma = sigma_l / (a_l * a_l);
  double alpha = alpha_l;
  double r_0 = r_0_l * a_l;

  // Allocate matrices
  gsl_matrix *H_Cornell = gsl_matrix_alloc(N_MAX, N_MAX);
  gsl_matrix *N_n_n_matrix = gsl_matrix_alloc(N_MAX, N_MAX);
  printf("V0FIT: %f\nsigma: %f\nalpha: %f\nl: %d\n", V_0_fit, sigma, alpha, l);

  // Fill matrices
  for (int i = 0; i < N_MAX; i++) {
    int n = i + 1;
    for (int j = 0; j <= i; j++) {
      int n_prime = j + 1;

      double H_ij = -T_n_n_prime(n, n_prime, l) +
                    V_0_fit * N_n_n_prime(n, n_prime, l) +
                    sigma * V_r_n_n_prime(n, n_prime, l) -
                    alpha * V_1_r_n_n_prime(n, n_prime, l);

      double N_nn = N_n_n_prime(n, n_prime, l);

      gsl_matrix_set(H_Cornell, i, j, H_ij);
      gsl_matrix_set(H_Cornell, j, i, H_ij);
      gsl_matrix_set(N_n_n_matrix, i, j, N_nn);
      gsl_matrix_set(N_n_n_matrix, j, i, N_nn);
    }
  }

  // Solve generalized eigenvalue problem
  gsl_vector *eval = gsl_vector_alloc(N_MAX);
  gsl_matrix *evec = gsl_matrix_alloc(N_MAX, N_MAX);

  gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(N_MAX);

  // Make copies of matrices as they will be modified by the eigenvalue solver
  gsl_matrix *H_copy = gsl_matrix_alloc(N_MAX, N_MAX);
  gsl_matrix *N_copy = gsl_matrix_alloc(N_MAX, N_MAX);

  gsl_matrix_memcpy(H_copy, H_Cornell);
  gsl_matrix_memcpy(N_copy, N_n_n_matrix);

  gsl_eigen_gensymmv(H_copy, N_copy, eval, evec, w);

  // Sort eigenvalues
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  // Normalize eigenvectors
  for (int i = 0; i < N_MAX; i++) {
    double norm = 0.0;

    // Calculate c† * N * c
    for (int j = 0; j < N_MAX; j++) {
      for (int k = 0; k < N_MAX; k++) {
        norm += gsl_matrix_get(evec, j, i) *
                gsl_matrix_get(N_n_n_matrix, j, k) * gsl_matrix_get(evec, k, i);
      }
    }

    norm = sqrt(norm);

    // Normalize eigenvector
    for (int j = 0; j < N_MAX; j++) {
      gsl_matrix_set(evec, j, i, gsl_matrix_get(evec, j, i) / norm);
    }
  }

  // Print eigenvalues
  printf("Eigenvalues:\n");
  for (int i = 0; i < N_MAX; i++) {
    printf("%d: %f\n", i, gsl_vector_get(eval, i));
  }

  // Free resources
  gsl_matrix_free(H_Cornell);
  gsl_matrix_free(N_n_n_matrix);
  gsl_matrix_free(H_copy);
  gsl_matrix_free(N_copy);
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
  gsl_eigen_gensymmv_free(w);
}

int main() {
  WaveFunction *wf = WFnew(0, 20, 40);
  // printf("Eigenvalues:\n");
  // for (int i = 0; i < N_MAX; i++) {
  //   printf("%d: %f\n", i, gsl_vector_get(wf->E_solution, i));
  // }
  double complex psi = psi_n(wf, 1, 1, 0);
  printf("%f + Im%f\n", GSL_REAL(psi), GSL_IMAG(psi));
  double r = 0;
  while (1) {
    scanf("%lf", &r);
    if (r < 0)
      break;
    psi = psi_n(wf, r, 1, 0);
    printf("%f + Im%f\n", GSL_REAL(psi), GSL_IMAG(psi));
  }
  WFfree(wf);
  // test();
  return 0;
}
