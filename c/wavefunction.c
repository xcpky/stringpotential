#include "wavefunction.h"
#include <complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include <stddef.h>
#include <sys/types.h>

WaveFunction *WFnew(uint64_t l, double Lambda, uint64_t rNgauss) {
  WaveFunction *self = (WaveFunction *)malloc(sizeof(WaveFunction));
  if (!self) {
    return NULL;
  }
  self->l = l;
  self->Lambda = Lambda;
  self->Ngauss = rNgauss;
  self->table = gsl_integration_glfixed_table_alloc(rNgauss);
  self->xi = (double *)malloc(sizeof(double) * rNgauss);
  self->wi = (double *)malloc(sizeof(double) * rNgauss);
  for (uint64_t i = 0; i < rNgauss; i += 1) {
    gsl_integration_glfixed_point(0, Lambda, i, &self->xi[i], &self->wi[i], self->table);
  }
  self->c_solution = gsl_matrix_alloc(N_MAX, N_MAX);
  self->E_solution = gsl_vector_alloc(N_MAX);
  build(self);
  return self;
};
void WFfree(WaveFunction *self) {
  gsl_integration_glfixed_table_free(self->table);
  gsl_matrix_free(self->c_solution);
  gsl_vector_free(self->E_solution);
  free(self->xi);
  free(self->wi);
};

void build(WaveFunction *self) {
  gsl_matrix *H_Cornell = gsl_matrix_alloc(N_MAX, N_MAX);
  gsl_matrix *N_n_n_matrix = gsl_matrix_alloc(N_MAX, N_MAX);

  uint64_t l = self->l;
  // Fill matrices
  for (int i = 0; i < N_MAX; i++) {
    int n = i + 1;
    for (int j = 0; j <= i; j++) {
      int n_prime = j + 1;

      double H_ij = -T_n_n_prime(n, n_prime, l) +
                    V0FIT * N_n_n_prime(n, n_prime, l) +
                    SIGMA * V_r_n_n_prime(n, n_prime, l) -
                    ALPHA * V_1_r_n_n_prime(n, n_prime, l);

      double N_nn = N_n_n_prime(n, n_prime, l);

      gsl_matrix_set(H_Cornell, i, j, H_ij);
      gsl_matrix_set(H_Cornell, j, i, H_ij);
      gsl_matrix_set(N_n_n_matrix, i, j, N_nn);
      gsl_matrix_set(N_n_n_matrix, j, i, N_nn);
    }
  }
  gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(N_MAX);
  gsl_matrix *N_copy = gsl_matrix_alloc(N_MAX, N_MAX);
  gsl_matrix_memcpy(N_copy, N_n_n_matrix);
  gsl_eigen_gensymmv(H_Cornell, N_copy, self->E_solution, self->c_solution, w);

  // Sort eigenvalues
  gsl_eigen_symmv_sort(self->E_solution, self->c_solution,
                       GSL_EIGEN_SORT_VAL_ASC);
  for (int i = 0; i < N_MAX; i++) {
    double norm = 0.0;

    // Calculate c† * N * c
    for (int j = 0; j < N_MAX; j++) {
      for (int k = 0; k < N_MAX; k++) {
        norm += gsl_matrix_get(self->c_solution, j, i) *
                gsl_matrix_get(N_n_n_matrix, j, k) *
                gsl_matrix_get(self->c_solution, k, i);
      }
    }

    norm = sqrt(norm);

    // Normalize eigenvector
    for (int j = 0; j < N_MAX; j++) {
      gsl_matrix_set(self->c_solution, j, i,
                     gsl_matrix_get(self->c_solution, j, i) / norm);
    }
  }
  gsl_matrix_free(H_Cornell);
  gsl_matrix_free(N_n_n_matrix);
  gsl_matrix_free(N_copy);
  gsl_eigen_gensymmv_free(w);
}

// Get dimensions of c_solution matrix
void WF_get_c_solution_dims(WaveFunction *self, size_t *rows, size_t *cols) {
  *rows = self->c_solution->size1;
  *cols = self->c_solution->size2;
}

// Get a value from c_solution matrix
double *WF_get_c_solution_data(WaveFunction *self) {
  return self->c_solution->data;
}

// Get length of E_solution vector
size_t WF_get_E_solution_length(WaveFunction *self) {
  return self->E_solution->size;
}

// Get a value from E_solution vector
double *WF_get_E_solution_data(WaveFunction *self) {
  return self->E_solution->data;
}

size_t WF_get_c_solution_tda(WaveFunction *self) {
  return self->c_solution->tda;
}

// Wavefunction definition
complex double psi_n(WaveFunction *self, double r, uint64_t n, double theta) {
  complex double psi = 0.0;

  for (int i = 0; i < N_MAX; i++) {
    int n_prime = i + 1;

    complex double Y_lm = gsl_sf_legendre_sphPlm(self->l, 0, cos(theta));

    psi += gsl_matrix_get(self->c_solution, i, n - 1) * N_nl(n_prime, self->l) *
           pow(r, self->l) * exp(-nu_n(n_prime) * r * r) * Y_lm;
  }

  return psi;
}

complex double psi_n_ftcomplex(WaveFunction *self, double complex p, uint64_t n) {
  const uint64_t Ngauss = self->Ngauss;
  const uint64_t l = self->l;
  double *xi = self->xi;
  double *wi = self->wi;
  complex double psi = 0.0;
  for (size_t nidx = 0; nidx < N_MAX; nidx += 1) {
    complex double quad = 0 + 0 * I;
    for (size_t i = 0; i < Ngauss; i += 1) {
      quad += integrand_complex(xi[i], p, nidx + 1, l) * wi[i];
    }
    psi += gsl_matrix_get(self->c_solution, nidx, n - 1) * N_nl(nidx + 1, l) *
           quad;
  }
  return sqrt(2 * l + 1) * 2 * sqrt(PI) * cpow(I, l + 0*I) * psi;
}

complex double psi_n_ft(WaveFunction *self, double p, uint64_t n) {
  const uint64_t Ngauss = self->Ngauss;
  const uint64_t l = self->l;
  double *xi = self->xi;
  double *wi = self->wi;
  complex double psi = 0.0;
  for (size_t nidx = 0; nidx < N_MAX; nidx += 1) {
    complex double quad = 0 + 0 * I;
    for (size_t i = 0; i < Ngauss; i += 1) {
      quad += integrand(xi[i], p, nidx + 1, l) * wi[i];
    }
    psi += gsl_matrix_get(self->c_solution, nidx, n - 1) * N_nl(nidx + 1, l) *
           quad;
  }
  return sqrt(2 * l + 1) * 2 * sqrt(PI) * cpow(I, l + 0*I) * psi;
}

void psi_n_batch(WaveFunction *self, const double *r_values,
                 double complex *results, size_t num_points, uint64_t n,
                 double theta) {
  for (size_t i = 0; i < num_points; i++) {
    results[i] = psi_n(self, r_values[i], n, theta);
  }
}

void psi_n_ft_batch(WaveFunction *self, const double *p_values,
                    double complex *results, size_t num_points, uint64_t n) {
  for (size_t i = 0; i < num_points; i++) {
    results[i] = psi_n_ft(self, p_values[i], n);
  }
}
