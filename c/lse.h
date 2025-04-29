#ifndef LSE_H
#define LSE_H
#include "constants.h"
#include "wavefunction.h"
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

typedef gsl_matrix_complex matrix;
#define inverse gsl_complex_inverse
#define add gsl_complex_add
#define mul gsl_complex_mul
#define mul_real gsl_complex_mul_real
#define sub gsl_complex_sub
#define matrix_alloc gsl_matrix_complex_alloc
#define matrix_free gsl_matrix_complex_free
#define matrix_set_zero gsl_matrix_complex_set_zero
#define matrix_set gsl_matrix_complex_set
#define matrix_memcpy gsl_matrix_complex_memcpy
#define matrix_scale gsl_matrix_complex_scale
#define matrix_get gsl_matrix_complex_get

typedef struct {
  size_t Ngauss;
  double Lambda;
  double epsilon;
  double E;
  matrix *T;
  matrix *V;
  matrix *G;
  double *xi;
  double *wi;
  double complex x0[2];
  gsl_integration_glfixed_table *table;
  WaveFunction *wf;
  double complex **psi_n_mat;
  double *E_vec;
} LSE;

LSE *lse_malloc(size_t Ngauss, double Lambda, double epsilon);
int lse_compute(LSE *app, double E);
void lse_free(LSE *app);
double complex *lse_get_g_data(LSE *app);
void lse_get_g_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_v_data(LSE *app);
void lse_get_v_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_t_data(LSE *app);
void lse_get_t_size(LSE *app, unsigned int *rows, unsigned int *cols);

// LSE methods
int lse_gmat(LSE *self);
int lse_vmat(LSE *self);
int lse_tmat(LSE *self);
void lse_refresh(LSE *self, double E);

// Complex square root function
static inline complex double xsqrt(complex double x) {
  if (cimag(x) >= 0)
    return csqrt(x);
  else
    return -csqrt(x - 0 * I);
}

static inline double square(double x) { return x * x; }

static inline double csquare(double complex x) { return x * x; }

// Contact term functions
static inline double Ctct_00(double g_C) { return -3 * 2 * g_C; }

static inline double Ctct_01(double g_C) {
  return pow(1.0 / 2.0, 3.0 / 2.0) * 4 * g_C;
}

static inline double Ctct_11(double g_C) { return g_C; }

// Omega functions
static inline double omega_00(double complex p, double complex pprime) {
  return 2 * m_B + (p * p + pprime * pprime) / (2 * m_B);
}

static inline double omega_01(double complex p, double complex pprime) {
  return m_B + pprime * pprime / (2 * m_B) + m_B_s + p * p / (2 * m_B_s);
}

static inline double omega_10(double complex p, double complex pprime) {
  return m_B_s + pprime * pprime / (2 * m_B_s) + m_B + p * p / (2 * m_B);
}

static inline double omega_11(double complex p, double complex pprime) {
  return 2 * m_B_s + (p * p + pprime * pprime) / (2 * m_B_s);
}

static inline double omegaprime_00(double complex p, double complex pprime) {
  return 2 * m_B_star + (p * p + pprime * pprime) / (2 * m_B_star);
}

static inline double omegaprime_01(double complex p, double complex pprime) {
  return m_B_star + pprime * pprime / (2 * m_B_star) + m_B_star_s +
         p * p / (2 * m_B_star_s);
}

static inline double omegaprime_10(double complex p, double complex pprime) {
  return m_B_star_s + pprime * pprime / (2 * m_B_star_s) + m_B_star +
         p * p / (2 * m_B_star);
}

static inline double omegaprime_11(double complex p, double complex pprime) {
  return 2 * m_B_star_s + (p * p + pprime * pprime) / (2 * m_B_star_s);
}

static inline double complex O_00(double E, double complex p, double complex pprime, double m) {
  if (fabs(p) <= 1e-8) {
    p += 1e-6;
  }
  if (fabs(pprime) <= 1e-8) {
    pprime += 1e-6;
  }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_00(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_00(p, pprime) + 0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_00(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_00(p, pprime) +
               0 * I)) *
         4 * square(g_pi) / square(f_pi);
}

static inline double complex O_01(double E, double complex p, double complex pprime, double m) {
  if (fabs(p) <= 1e-8) {
    p += 1e-6;
  }
  if (fabs(pprime) <= 1e-8) {
    pprime += 1e-6;
  }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_01(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_01(p, pprime) + 0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_01(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_01(p, pprime) +
               0 * I)) *
         4 * square(g_pi) / square(f_pi);
}

static inline double complex O_10(double E, double complex p, double complex pprime, double m) {
  if (fabs(p) <= 1e-8) {
    p += 1e-6;
  }
  if (fabs(pprime) <= 1e-8) {
    pprime += 1e-6;
  }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_10(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_10(p, pprime) + 0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_10(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_10(p, pprime) +
               0 * I)) *
         4 * square(g_pi) / square(f_pi);
}

static inline double complex O_11(double E, double complex p, double complex pprime, double m) {
  if (fabs(p) <= 1e-8) {
    p += 1e-6;
  }
  if (fabs(pprime) <= 1e-8) {
    pprime += 1e-6;
  }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_11(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_11(p, pprime) + 0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_11(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_11(p, pprime) +
               0 * I)) *
         4 * square(g_pi) / square(f_pi);
}

static inline double complex V_OME_00(double E, double complex p, double complex pprime) {
  return -3 * (3 * O_00(E, p, pprime, m_pi) + O_00(E, p, pprime, m_eta) / 3);
}

static inline double complex V_OME_01(double E, double complex p, double complex pprime) {
  return pow(2, 3. / 2) * O_01(E, p, pprime, m_K);
}

static inline double complex V_OME_10(double E, double complex p, double complex pprime) {
  return pow(2, 3. / 2) * O_10(E, p, pprime, m_K);
}

static inline double complex V_OME_11(double E, double complex p, double complex pprime) {
  return 2. / 3 * O_11(E, p, pprime, m_eta);
}

double complex V_QM(LSE *self, double E, uint64_t p, uint64_t pprime);

#endif // !LSE_H
