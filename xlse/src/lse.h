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
  double complex E;
  matrix *T;
  matrix *V;
  matrix *G;
  matrix *reg;
  double *xi;
  double *wi;
  double complex x0[2];
  gsl_integration_glfixed_table *table;
  WaveFunction *wf;
  double complex **psi_n_mat;
  double complex *psi_raw;
  double *E_vec;
} LSE;

LSE *lse_malloc(size_t pNgauss, double Lambda, double epsilon);
int lse_compute(LSE *app, double complex E, int64_t rs);
void lse_free(LSE *app);
double complex *lse_get_g_data(LSE *app);
void lse_get_g_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_v_data(LSE *app);
void lse_get_v_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_t_data(LSE *app);
void lse_get_t_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_ivg_data(LSE *app);
void lse_get_ivg_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_iivg_data(LSE *app);
void lse_get_iivg_size(LSE *app, unsigned int *rows, unsigned int *cols);
double complex *lse_get_psi(LSE *self);
void lse_get_psi_size(LSE *self, unsigned int *rows, unsigned int *cols);
double *lse_get_E(LSE *self);
void lse_get_E_size(unsigned int *levels);
void lse_get_M_size(LSE *app, unsigned int *rows, unsigned int *cols);

// LSE methods
int lse_gmat(LSE *self);
int lse_vmat(LSE *self);
int lse_tmat(LSE *self);
double complex lse_detImVG(LSE *self, double complex E);
void lse_refresh(LSE *self, double complex E, int64_t rs);

// Complex square root function
static inline double complex xsqrt(double complex x) {
  if (cimag(x) >= 0)
    return csqrt(x);
  else
    return -csqrt(x - 0 * I);
}

static inline double square(double x) { return x * x; }

static inline double complex csquare(double complex x) { return x * x; }

// Contact term functions
static inline double Ctct_00(double g_C) { return -3 * 2 * g_C; }

static inline double Ctct_01(double g_C) {
  return pow(1.0 / 2.0, 3.0 / 2.0) * 4 * g_C;
}

static inline double Ctct_10(double g_C) {
  return pow(1.0 / 2.0, 3.0 / 2.0) * 4 * g_C;
}

static inline double Ctct_11(double g_C) { return g_C; }

// Omega functions
static inline double complex omega_00(double complex p, double complex pprime) {
  return 2 * m_B + (p * p + pprime * pprime) / (2 * m_B);
}

static inline double complex omega_01(double complex p, double complex pprime) {
  return m_B + pprime * pprime / (2 * m_B) + m_B_s + p * p / (2 * m_B_s);
}

static inline double complex omega_10(double complex p, double complex pprime) {
  return m_B_s + pprime * pprime / (2 * m_B_s) + m_B + p * p / (2 * m_B);
}

static inline double complex omega_11(double complex p, double complex pprime) {
  return 2 * m_B_s + (p * p + pprime * pprime) / (2 * m_B_s);
}

static inline double complex omegaprime_00(double complex p,
                                           double complex pprime) {
  return 2 * m_B_star + (p * p + pprime * pprime) / (2 * m_B_star);
}

static inline double complex omegaprime_01(double complex p,
                                           double complex pprime) {
  return m_B_star + pprime * pprime / (2 * m_B_star) + m_B_star_s +
         p * p / (2 * m_B_star_s);
}

static inline double complex omegaprime_10(double complex p,
                                           double complex pprime) {
  return m_B_star_s + pprime * pprime / (2 * m_B_star_s) + m_B_star +
         p * p / (2 * m_B_star);
}

static inline double complex omegaprime_11(double complex p,
                                           double complex pprime) {
  return 2 * m_B_star_s + (p * p + pprime * pprime) / (2 * m_B_star_s);
}

#define DEFINE_O_FUNCTION(suffix)                                              \
  static inline double complex O_##suffix(double complex E, double complex p,  \
                                          double complex pprime, double m) {   \
    if (cabs(p) <= 1e-8) {                                                     \
      p += 1e-6;                                                               \
    }                                                                          \
    if (cabs(pprime) <= 1e-8) {                                                \
      pprime += 1e-6;                                                          \
    }                                                                          \
    return square(g_pi) / square(f_pi) * -1. / 4. / p / pprime *           \
           (clog((E - (m + csquare(p - pprime) / 2 / m) -                      \
                  omega_##suffix(p, pprime)) /                                 \
                 (E - (m + csquare(p + pprime) / 2 / m) -                      \
                  omega_##suffix(p, pprime) + 0 * I)) +                        \
            clog((E - (m + csquare(p - pprime) / 2 / m) -                      \
                  omegaprime_##suffix(p, pprime)) /                            \
                 (E - (m + csquare(p + pprime) / 2 / m) -                      \
                  omegaprime_##suffix(p, pprime) + 0 * I)));                   \
  }

DEFINE_O_FUNCTION(00);
DEFINE_O_FUNCTION(01);
DEFINE_O_FUNCTION(10);
DEFINE_O_FUNCTION(11);

static inline double complex V_OME_00(double complex E, double complex p,
                                      double complex pprime) {
  return -3 * (3 * O_00(E, p, pprime, m_pi) + O_00(E, p, pprime, m_eta) / 3);
}

static inline double complex V_OME_01(double complex E, double complex p,
                                      double complex pprime) {
  return pow(2, 3. / 2) * O_01(E, p, pprime, m_K);
}

static inline double complex V_OME_10(double complex E, double complex p,
                                      double complex pprime) {
  return pow(2, 3. / 2) * O_10(E, p, pprime, m_K);
}

static inline double complex V_OME_11(double complex E, double complex p,
                                      double complex pprime) {
  return 2. / 3 * O_11(E, p, pprime, m_eta);
}

double complex V_QM_00(LSE *self, uint64_t p, uint64_t pprime);
double complex V_QM_01(LSE *self, uint64_t p, uint64_t pprime);
double complex V_QM_10(LSE *self, uint64_t p, uint64_t pprime);
double complex V_QM_11(LSE *self, uint64_t p, uint64_t pprime);

#define DEFINE_VQM(alpha, beta)                                                \
  double complex V_QM_##alpha##beta(LSE *self, uint64_t p, uint64_t pprime) {  \
    double complex res = 0 + 0I;                                               \
    __auto_type E = self->E;                                                   \
    for (size_t i = 0; i < N_MAX; i++) {                                       \
      res -= self->psi_n_mat[i][p] * conj(self->psi_n_mat[i][pprime]) /        \
             (E - self->E_vec[i] + self->epsilon * I);                         \
    }                                                                          \
    return res * g##alpha * g##beta;                                           \
  }

static inline double complex curlO(double complex p, double complex pprime,
                                   double m) {
  if (cabs(p) < 1e-8)
    p += 1e-6;
  if (cabs(pprime) < 1e-8)
    pprime += 1e-6;
  return 4 * square(g_pi) / square(f_pi) *
         (1 - square(m) / (4 * p * pprime) *
                  (clog((csquare(p) + csquare(pprime) + 2 * p * pprime +
                         csquare(m))) -
                   clog((csquare(p) + csquare(pprime) - 2 * p * pprime +
                         csquare(m)))));
}

static inline double complex V_curlOME_00(double complex E, double complex p,
                                          double complex pprime) {
  return -3 * (3 * curlO(p, pprime, m_pi) + curlO(p, pprime, m_eta) / 3);
}

static inline double complex V_curlOME_01(double complex E, double complex p,
                                          double complex pprime) {
  return pow(2, 3. / 2) * curlO(p, pprime, m_K);
}

static inline double complex V_curlOME_10(double complex E, double complex p,
                                          double complex pprime) {
  return pow(2, 3. / 2) * curlO(p, pprime, m_K);
}

static inline double complex V_curlOME_11(double complex E, double complex p,
                                          double complex pprime) {
  return 2. / 3 * curlO(p, pprime, m_eta);
}

#define DEFINE_V_FUNCTION(suffix)                                              \
  static inline gsl_complex V##suffix(LSE *self, double complex p,             \
                                      double complex pprime) {                 \
    __auto_type E = self->E;                                                   \
    return V_curlOME_##suffix(E, p, pprime) + Ctct_##suffix(g_c);                                             \
  }

DEFINE_V_FUNCTION(00);
DEFINE_V_FUNCTION(01);
DEFINE_V_FUNCTION(10);
DEFINE_V_FUNCTION(11);

#endif // !LSE_H
