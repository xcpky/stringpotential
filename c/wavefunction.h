#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "constants.h"
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector_double.h>
#include <stdint.h>

typedef struct {
  uint64_t l;
  double Lambda;
  uint64_t Ngauss;
  double *xi;
  double *wi;
  gsl_integration_glfixed_table *table;
  gsl_matrix *c_solution;
  gsl_vector *E_solution;
} WaveFunction;

WaveFunction *WFnew(uint64_t l, double Lambda, uint64_t Ngauss);
void WFfree(WaveFunction *self);
void build(WaveFunction *self);
void WF_get_c_solution_dims(WaveFunction *self, size_t *rows, size_t *cols);
double *WF_get_c_solution_data(WaveFunction *self);
size_t WF_get_c_solution_tda(WaveFunction *self);
size_t WF_get_E_solution_length(WaveFunction *self);
double *WF_get_E_solution_data(WaveFunction *self);
complex double psi_n(WaveFunction *self, double r, uint64_t n, double theta);
complex double psi_n_ft(WaveFunction *self, double p, uint64_t n);
void psi_n_batch(WaveFunction *self, const double *r_values,
                 double complex *results, size_t num_points, uint64_t n,
                 double theta);
void psi_n_ft_batch(WaveFunction *self, const double *p_values,
                    double complex *results, size_t num_points, uint64_t n);

// Constants
#define N_MAX 30
#define R_1 (0.02 * 5.068)   // GeV^-1
#define R_N_MAX (30 * 5.068) // GeV^-1
#define C_T ((-1) * 0.19732 * 0.19732 / (0.5 * 4.18) * 5.068)
#define PI 3.14159265358979323846
#define V0FIT 0.4
#define SIGMA_L 0.0199179973550142
#define SIGMA (SIGMA_L / a_l / a_l)
#define ALPHA 0.476814032326273

static inline double doublefactorial(int n) {
  if (n == 0 || n == 1)
    return 1.0;
  return n * doublefactorial(n - 2);
}

// Radial coordinate functions
static inline double r_n(int n) {
  return R_1 * pow(R_N_MAX / R_1, ((double)(n - 1) / (N_MAX - 1)));
}

static inline double nu_n(int n) { return 1.0 / (r_n(n) * r_n(n)); }

// Normalization functions
static inline double N_nl(int n, int l) {
  return sqrt(pow(2, l + 2) * pow(2 * nu_n(n), l + 1.5) /
              (sqrt(PI) * doublefactorial(2 * l + 1)));
}

static inline double N_n_n_prime(int n, int n_prime, int l) {
  return pow((2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime))),
             l + 1.5);
}

// Matrix element functions
static inline double T_n_n_prime(int n, int n_prime, int l) {
  return C_T * (2 * l + 3) * nu_n(n) * nu_n(n_prime) /
         (nu_n(n) + nu_n(n_prime)) * N_n_n_prime(n, n_prime, l);
}

static inline double V_r_n_n_prime(int n, int n_prime, int l) {
  return 1.0 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) *
         doublefactorial(2 * l + 2) / doublefactorial(2 * l + 1) *
         N_n_n_prime(n, n_prime, l);
}

static inline double V_1_r_n_n_prime(int n, int n_prime, int l) {
  return 2.0 / sqrt(PI) * pow(2, l) * gsl_sf_fact(l) /
         doublefactorial(2 * l + 1) * sqrt(nu_n(n) + nu_n(n_prime)) *
         N_n_n_prime(n, n_prime, l);
}

// Modified normalization functions
static inline double N_nl_tilde(int n, int l) {
  return sqrt(pow(2, l + 3) * pow(2 * nu_n(n), l + 2.5) /
              (sqrt(PI) * doublefactorial(2 * l + 3)));
}

static inline double N_n_n_prime_tilde(int n, int n_prime, int l) {
  return pow((2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime))),
             l + 2.5);
}

static inline double V_r_n_n_prime_tilde(int n, int n_prime, int l) {
  return 1.0 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) *
         doublefactorial(2 * l + 4) / doublefactorial(2 * l + 3) *
         N_n_n_prime_tilde(n, n_prime, l);
}

static inline complex double integrand(double r, double p, int n, int l) {
  return gsl_sf_bessel_jl(l, p * r) * pow(r, l+2) * exp(-nu_n(n) * r * r);
}
#endif // !WAVEFUNCTION_H
