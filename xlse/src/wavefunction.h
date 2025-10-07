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
    uint l;
    double rLambda;
    uint64_t rNgauss;
    double *xi;
    double *wi;
    gsl_integration_glfixed_table *table;
    gsl_matrix *c_solution;
    gsl_vector *E_solution;
    // double complex **psi_mat;
    // double complex *psi_raw;
} WaveFunction;

WaveFunction *WFnew(uint64_t l, double rLambda, uint64_t rNgauss);
void WFfree(WaveFunction *self);
void build(WaveFunction *self);
void WF_get_c_solution_dims(WaveFunction *self, size_t *rows, size_t *cols);
double *WF_get_c_solution_data(WaveFunction *self);
size_t WF_get_c_solution_tda(WaveFunction *self);
size_t WF_get_E_solution_length(WaveFunction *self);
double *WF_get_E_solution_data(WaveFunction *self);
void wf_refresh(WaveFunction *, double);
complex double psi_n(WaveFunction *self, double r, uint64_t n, double theta);
complex double psi_n_ft(WaveFunction *self, double p, uint64_t n);
complex double psi_n_ftcomplex(WaveFunction *self, double complex p, uint64_t n);
void psi_n_batch(WaveFunction *self, const double *r_values, double complex *results, size_t num_points, uint64_t n,
                 double theta);
void psi_n_ft_batch(WaveFunction *self, const double *p_values, double complex *results, size_t num_points, uint64_t n);

// Radial coordinate functions
static inline double r_n(uint n) { return R_1 * pow(R_N_MAX / R_1, ((double)(n - 1) / (N_MAX - 1))); }

static inline double nu_n(uint n) { return 1.0 / (r_n(n) * r_n(n)); }

// Normalization functions
static inline double N_nl(uint n, uint l) {
    return sqrt(pow(2, l + 2) * pow(2 * nu_n(n), l + 1.5) / (sqrt(PI) * gsl_sf_doublefact(2 * l + 1)));
}

static inline double N_n_n_prime(uint n, uint n_prime, uint l) {
    return pow((2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime))), l + 1.5);
}

// Matrix element functions
static inline double T_n_n_prime(uint n, uint n_prime, uint l) {
    return C_T * (2 * l + 3) * nu_n(n) * nu_n(n_prime) / (nu_n(n) + nu_n(n_prime)) * N_n_n_prime(n, n_prime, l);
}

static inline double V_r_n_n_prime(uint n, uint n_prime, uint l) {
    return 1.0 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) * gsl_sf_doublefact(2 * l + 2) / gsl_sf_doublefact(2 * l + 1) *
           N_n_n_prime(n, n_prime, l);
}

static inline double V_1_r_n_n_prime(uint n, uint n_prime, uint l) {
    return 2.0 / sqrt(PI) * pow(2, l) * gsl_sf_fact(l) / gsl_sf_doublefact(2 * l + 1) * sqrt(nu_n(n) + nu_n(n_prime)) *
           N_n_n_prime(n, n_prime, l);
}

// Modified normalization functions
static inline double N_nl_tilde(uint n, uint l) {
    return sqrt(pow(2, l + 3) * pow(2 * nu_n(n), l + 2.5) / (sqrt(PI) * gsl_sf_doublefact(2 * l + 3)));
}

static inline double N_n_n_prime_tilde(uint n, uint n_prime, uint l) {
    return pow((2 * sqrt(nu_n(n) * nu_n(n_prime)) / (nu_n(n) + nu_n(n_prime))), l + 2.5);
}

static inline double V_r_n_n_prime_tilde(uint n, uint n_prime, uint l) {
    return 1.0 / sqrt(2 * (nu_n(n) + nu_n(n_prime))) * gsl_sf_doublefact(2 * l + 4) / gsl_sf_doublefact(2 * l + 3) *
           N_n_n_prime_tilde(n, n_prime, l);
}

static inline complex double integrand_complex(double r, double complex p, uint n, uint l) {
    if (l == 0) {
        return csin(r * p) / p * r * exp(-nu_n(n) * r * r);
    } else {
        return (csin(r * p) / p / p - ccos(r * p) * r / p) * r * exp(-nu_n(n) * r * r);
    }
}

static inline complex double integrand(double r, double p, uint n, uint l) {
    if (l == 0) {
        return sin(r * p) / p * r * exp(-nu_n(n) * r * r);
    } else {
        return (sin(r * p) / p / p - cos(r * p) * r / p) * r * exp(-nu_n(n) * r * r);
    }
}

static inline double complex psi_test(double complex p) { return I * p * cexp(-p * p / 4 / M_PI); }

// static inline complex double integrand(double r, double p, int n, int l) {
//   return gsl_sf_bessel_jl(l, p * r) * pow(r, l + 2) * exp(-nu_n(n));
// }
#endif // !WAVEFUNCTION_H
