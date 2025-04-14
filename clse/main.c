/**
 * LSE (Lippmann-Schwinger Equation) Solver in C
 */

#include <complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef gsl_complex cpx;
typedef gsl_matrix_complex matrix;
#define inverse gsl_complex_inverse
#define add gsl_complex_add
#define mul gsl_complex_mul
#define mul_real gsl_complex_mul_real
#define sub gsl_complex_sub
#define comp gsl_complex_rect
#define matrix_alloc gsl_matrix_complex_alloc
#define matrix_free gsl_matrix_complex_free
#define matrix_set_zero gsl_matrix_complex_set_zero
#define matrix_set gsl_matrix_complex_set
#define matrix_memcpy gsl_matrix_complex_memcpy
#define matrix_scale gsl_matrix_complex_scale
#define matrix_get gsl_matrix_complex_get

#define V_OME(alpha, beta, E, p, pprime) V_OME_##alpha##_##beta(E, p, pprime)

#define V_OME_0_0(E, p, pprime)                                                \
  (gsl_complex_mul_real(                                                       \
      gsl_complex_add(gsl_complex_mul_real(O(0, 0, E, p, pprime, m_pi), 3),    \
                      gsl_complex_div_real(O(0, 0, E, p, pprime, m_eta), 3)),  \
      -3))
#define V_OME_0_1(E, p, pprime)                                                \
  (gsl_complex_mul_real(O(0, 1, E, p, pprime, m_K), pow(2, 3. / 2.)))
#define V_OME_1_0(E, p, pprime)                                                \
  (gsl_complex_mul_real(O(1, 0, E, p, pprime, m_K), pow(2, 3. / 2.)))
#define V_OME_1_1(E, p, pprime)                                                \
  (gsl_complex_mul_real(O(1, 1, E, p, pprime, m_eta), 2. / 3.))
#define OMEGA(alpha, beta, p, pprime) OMEGA_##alpha##_##beta(p, pprime)

// Define specialized versions for each alpha/beta combination
#define OMEGA_0_0(p, pprime)                                                   \
  (2 * m_B + ((p) * (p) + (pprime) * (pprime)) / 2 / m_B)
#define OMEGA_0_1(p, pprime)                                                   \
  (m_B + (pprime) * (pprime) / 2 / m_B + m_B_s + (p) * (p) / 2 / m_B_s)
#define OMEGA_1_0(p, pprime)                                                   \
  (m_B_s + (pprime) * (pprime) / 2 / m_B_s + m_B + (p) * (p) / 2 / m_B)
#define OMEGA_1_1(p, pprime)                                                   \
  (2 * m_B_s + ((p) * (p) + (pprime) * (pprime)) / 2 / m_B_s)

#define OMEGA_PRIME(alpha, beta, p, pprime)                                    \
  OMEGA_PRIME_##alpha##_##beta(p, pprime)

// Define specialized versions for each alpha/beta combination
#define OMEGA_PRIME_0_0(p, pprime)                                             \
  (2 * m_B_star + ((p) * (p) + (pprime) * (pprime)) / 2 / m_B_star)
#define OMEGA_PRIME_0_1(p, pprime)                                             \
  (m_B_star + (pprime) * (pprime) / 2 / m_B_star + m_B_star_s +                \
   (p) * (p) / 2 / m_B_star_s)
#define OMEGA_PRIME_1_0(p, pprime)                                             \
  (m_B_star_s + (pprime) * (pprime) / 2 / m_B_star_s + m_B_star +              \
   (p) * (p) / 2 / m_B_star)
#define OMEGA_PRIME_1_1(p, pprime)                                             \
  (2 * m_B_star_s + ((p) * (p) + (pprime) * (pprime)) / 2 / m_B_star_s)

#define O(alpha, beta, E, p, pprime, m) O_##alpha##_##beta(E, p, pprime, m)

// Define specialized versions of O for each alpha/beta combination
#define O_0_0(E, p, pprime, m)                                                 \
  O_IMPL(E, p, pprime, m, OMEGA_0_0(p, pprime), OMEGA_PRIME_0_0(p, pprime))
#define O_0_1(E, p, pprime, m)                                                 \
  O_IMPL(E, p, pprime, m, OMEGA_0_1(p, pprime), OMEGA_PRIME_0_1(p, pprime))
#define O_1_0(E, p, pprime, m)                                                 \
  O_IMPL(E, p, pprime, m, OMEGA_1_0(p, pprime), OMEGA_PRIME_1_0(p, pprime))
#define O_1_1(E, p, pprime, m)                                                 \
  O_IMPL(E, p, pprime, m, OMEGA_1_1(p, pprime), OMEGA_PRIME_1_1(p, pprime))

// Implementation macro that takes omega and omega_prime directly
#define O_IMPL(E, p, pprime, m, omega, omega_prime)                            \
  ({                                                                           \
    double tmp1 = (E) - ((m) + square((p) - (pprime)) / 2 / (m));              \
    double tmp2 = (E) - ((m) + square((p) + (pprime)) / 2 / (m));              \
    gsl_complex log1 = gsl_complex_log(                                        \
        gsl_complex_rect((tmp1 - (omega)) / (tmp2 - (omega)), 0));             \
    gsl_complex log2 = gsl_complex_log(                                        \
        gsl_complex_rect((tmp1 - (omega_prime)) / (tmp2 - (omega_prime)), 0)); \
    gsl_complex_mul_real(gsl_complex_add(log1, log2),                          \
                         -1 / (p) / (pprime) / 4);                             \
  })

// Constants
const bool debug = false;
const double f_pi = 0.092;
const double g_pi = 0.5704;
const double m_pi = 0.138039407;
const double m_K = 0.498;
const double m_eta = 0.548;
const double m_eta_s = 0.690625803;
const double m_B = 5.27934;
const double m_B_star = 5.32471;
const double m_B_s = 5.36692;
const double m_B_star_s = 5.4154;
const double m11 = m_B;
const double m12 = m_B_star;
const double m21 = m_B_s;
const double m22 = m_B_star_s;
const double delta[2] = {0, m21 + m22 - m11 - m12};
const double mu[2] = {m11 * m12 / (m11 + m12), m21 *m22 / (m21 + m22)};

// Utility functions
static inline double square(double x) { return x * x; }

static inline gsl_complex xsqrt(gsl_complex x) {
  if (GSL_IMAG(x) >= 0) {
    return gsl_complex_sqrt(x);
  } else {
    return gsl_complex_negative(gsl_complex_sqrt(x));
  }
}

// Forward declarations for functions
// double Omega(int alpha, int beta, double p, double pprime);
// double Omega_prime(int alpha, int beta, double p, double pprime);
// double O(int alpha, int beta, double E, double p, double pprime, double m);
// double V_OME(int alpha, int beta, double E, double p, double pprime);

// LSE structure
typedef struct {
  size_t Ngauss;
  double Lambda;
  double epsilon;
  double E;
  gsl_matrix_complex *T;
  gsl_matrix_complex *V;
  gsl_matrix_complex *G;
  double *xi;
  double *wi;
  double x0[2];
  gsl_integration_glfixed_table *table;
} LSE;

// Function prototypes
LSE *lse_new(size_t Ngauss, double Lambda, double epsilon);
int lse_compute(LSE *app, double E);
void lse_deinit(LSE *app);
double *lse_get_g_data(LSE *app);
void lse_get_g_size(LSE *app, unsigned int *rows, unsigned int *cols);
double *lse_get_v_data(LSE *app);
void lse_get_v_size(LSE *app, unsigned int *rows, unsigned int *cols);
double *lse_get_t_data(LSE *app);
void lse_get_t_size(LSE *app, unsigned int *rows, unsigned int *cols);

// LSE methods
int lse_gmat(LSE *self, double E);
int lse_vmat(LSE *self, double E);
int lse_tmat(LSE *self, double E);
static inline void lse_refresh(LSE *self, double E);

// Implementation of LSE creation
LSE *lse_new(size_t Ngauss, double Lambda, double epsilon) {
  LSE *self = (LSE *)malloc(sizeof(LSE));
  if (!self)
    return NULL;

  const size_t n = 2 * (Ngauss + 1);

  self->Ngauss = Ngauss;
  self->Lambda = Lambda;
  self->epsilon = epsilon;
  // self->E = E;

  self->T = gsl_matrix_complex_alloc(n, n);
  self->V = gsl_matrix_complex_alloc(n, n);
  self->G = gsl_matrix_complex_alloc(n, n);

  if (!self->T || !self->V || !self->G) {
    if (self->T)
      gsl_matrix_complex_free(self->T);
    if (self->V)
      gsl_matrix_complex_free(self->V);
    if (self->G)
      gsl_matrix_complex_free(self->G);
    free(self);
    return NULL;
  }

  self->table = gsl_integration_glfixed_table_alloc(Ngauss);
  self->xi = (double *)malloc(Ngauss * sizeof(double));
  self->wi = (double *)malloc(Ngauss * sizeof(double));

  if (!self->table || !self->xi || !self->wi) {
    if (self->T)
      gsl_matrix_complex_free(self->T);
    if (self->V)
      gsl_matrix_complex_free(self->V);
    if (self->G)
      gsl_matrix_complex_free(self->G);
    if (self->table)
      gsl_integration_glfixed_table_free(self->table);
    if (self->xi)
      free(self->xi);
    if (self->wi)
      free(self->wi);
    free(self);
    return NULL;
  }

  for (size_t i = 0; i < Ngauss; i++) {
    gsl_integration_glfixed_point(0, Lambda, i, &self->xi[i], &self->wi[i],
                                  self->table);
  }

  // Initialize x0
  // for (int i = 0; i < 2; i++) {
  //   const double dE = E - delta[i];
  //   const double mU = mu[i];
  //   const double tmp = sqrt(2 * mU * dE);
  //   self->x0[i] = tmp;
  // }

  return self;
}

// Refresh LSE parameters
static inline void lse_refresh(LSE *self, double E) {
  // self->Lambda = Lambda;
  self->E = E;

  for (int i = 0; i < 2; i++) {
    const double dE = E - delta[i];
    const double mU = mu[i];
    const double tmp = sqrt(2 * mU * dE);
    self->x0[i] = tmp;
  }
}

// Free LSE resources
void lse_deinit(LSE *self) {
  if (!self)
    return;

  gsl_matrix_complex_free(self->T);
  gsl_matrix_complex_free(self->V);
  gsl_matrix_complex_free(self->G);
  gsl_integration_glfixed_table_free(self->table);
  free(self->xi);
  free(self->wi);
  free(self);
}

// Calculate G matrix
int lse_gmat(LSE *self, double E) {
  matrix_set_zero(self->G);
  lse_refresh(self, E);

  for (int i = 0; i < 2; i++) {
    const double dE = self->E - delta[i];
    const double mU = mu[i];
    const double x0 = self->x0[i];
    if (dE >= 0) {

      cpx int_val = comp(0, 0);

      for (size_t j = 0; j < self->Ngauss; j++) {
        double x = self->xi[j];
        double w = self->wi[j];

        cpx new_val = comp(dE - square(x) / 2 / mU, self->epsilon);
        new_val = inverse(new_val);
        new_val = mul_real(new_val, w);
        int_val = add(int_val, new_val);
      }

      cpx tmp = comp(mU * x0 * log((self->Lambda + x0) / (self->Lambda - x0)),
                     -mU * x0 * M_PI);
      int_val = mul_real(int_val, square(x0));
      tmp = sub(tmp, int_val);

      const size_t ii = self->Ngauss + i * (self->Ngauss + 1);
      matrix_set(self->G, ii, ii, mul_real(tmp, 1 / square(M_PI) / 2));
    }

    for (size_t m = 0; m < self->Ngauss; m++) {
      const size_t pos = m + i * (self->Ngauss + 1);
      cpx denominator = comp(dE - square(self->xi[m]) / 2 / mU, self->epsilon);
      cpx inv = inverse(denominator);
      cpx ele =
          mul_real(inv, square(self->xi[m]) * self->wi[m] / 2 / square(M_PI));

      matrix_set(self->G, pos, pos, ele);
    }
  }

  return 0;
}

// Calculate V matrix
int lse_vmat(LSE *self, double E) {
  lse_refresh(self, E);
  gsl_matrix_complex_set_zero(self->V);

  double **x = (double **)malloc(2 * sizeof(double *));
  if (!x)
    return -1;

  for (int i = 0; i < 2; i++) {
    x[i] = (double *)malloc((self->Ngauss + 1) * sizeof(double));
    if (!x[i]) {
      for (int j = 0; j < i; j++) {
        free(x[j]);
      }
      free(x);
      return -1;
    }

    // Copy xi values
    for (size_t j = 0; j < self->Ngauss; j++) {
      x[i][j] = self->xi[j];
    }

    // Set last element to x0
    x[i][self->Ngauss] = self->x0[i];
  }
  // Unrolled loops for alpha=0, beta=0
  for (size_t idx = 0; idx <= self->Ngauss; idx++) {
    double p[4];
    p[0] = x[0][idx];
    p[1] = x[0][idx];
    p[2] = x[1][idx];
    p[3] = x[1][idx];
    for (size_t jdx = 0; jdx <= self->Ngauss; jdx++) {
      double pprime = x[0][jdx];
      size_t i = idx + 0 * (self->Ngauss + 1);
      size_t j = jdx + 0 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_0_0(self->E, p[0], pprime));
      pprime = x[1][jdx];
      i = idx + 0 * (self->Ngauss + 1);
      j = jdx + 1 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_0_1(self->E, p[1], pprime));
      pprime = x[0][jdx];
      i = idx + 1 * (self->Ngauss + 1);
      j = jdx + 0 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_1_0(self->E, p[2], pprime));
      pprime = x[1][jdx];
      i = idx + 1 * (self->Ngauss + 1);
      j = jdx + 1 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_1_1(self->E, p[3], pprime));
    }
  }
  //
  // // Unrolled loops for alpha=0, beta=1
  // for (size_t idx = 0; idx <= self->Ngauss; idx++) {
  //   double p = x[0][idx];
  //   for (size_t jdx = 0; jdx <= self->Ngauss; jdx++) {
  //     double pprime = x[1][jdx];
  //     const size_t i = idx + 0 * (self->Ngauss + 1);
  //     const size_t j = jdx + 1 * (self->Ngauss + 1);
  //     gsl_matrix_complex_set(self->V, i, j, V_OME_0_1(self->E, p, pprime));
  //   }
  // }
  //
  // // Unrolled loops for alpha=1, beta=0
  // for (size_t idx = 0; idx <= self->Ngauss; idx++) {
  //   double p = x[1][idx];
  //   for (size_t jdx = 0; jdx <= self->Ngauss; jdx++) {
  //     double pprime = x[0][jdx];
  //     const size_t i = idx + 1 * (self->Ngauss + 1);
  //     const size_t j = jdx + 0 * (self->Ngauss + 1);
  //     gsl_matrix_complex_set(self->V, i, j, V_OME_1_0(self->E, p, pprime));
  //   }
  // }
  //
  // // Unrolled loops for alpha=1, beta=1
  // for (size_t idx = 0; idx <= self->Ngauss; idx++) {
  //   double p = x[1][idx];
  //   for (size_t jdx = 0; jdx <= self->Ngauss; jdx++) {
  //     double pprime = x[1][jdx];
  //     const size_t i = idx + 1 * (self->Ngauss + 1);
  //     const size_t j = jdx + 1 * (self->Ngauss + 1);
  //     gsl_matrix_complex_set(self->V, i, j, V_OME_1_1(self->E, p, pprime));
  //   }
  // }
  //
  // Clean up
  for (int i = 0; i < 2; i++) {
    free(x[i]);
  }
  free(x);

  return 0;
}

// Calculate T matrix
int lse_tmat(LSE *self, double E) {
  const size_t n = 2 * (self->Ngauss + 1);
  lse_refresh(self, E);

  // Step 1: Compute VG = V * G
  gsl_matrix_complex *VG = gsl_matrix_complex_alloc(n, n);
  if (!VG)
    return -1;

  gsl_matrix_complex_set_zero(VG);

  gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
  gsl_complex beta = gsl_complex_rect(0.0, 0.0);

  // Perform matrix multiplication: VG = V * G
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)n, (int)n, (int)n,
              &GSL_REAL(alpha), self->V->data, (int)self->V->tda, self->G->data,
              (int)self->G->tda, &GSL_REAL(beta), VG->data, (int)VG->tda);

  // Step 2: Compute I - VG
  gsl_matrix_complex *I_minus_VG = gsl_matrix_complex_alloc(n, n);
  if (!I_minus_VG) {
    gsl_matrix_complex_free(VG);
    return -1;
  }

  gsl_matrix_complex_memcpy(I_minus_VG, VG);
  gsl_matrix_complex_scale(I_minus_VG,
                           gsl_complex_rect(-1.0, 0.0)); // I_minus_VG = -VG

  // Add identity matrix: I_minus_VG = I - VG
  for (size_t i = 0; i < n; i++) {
    gsl_complex diag = gsl_matrix_complex_get(I_minus_VG, i, i);
    gsl_complex one = gsl_complex_rect(1.0, 0.0);
    gsl_matrix_complex_set(I_minus_VG, i, i, gsl_complex_add(diag, one));
  }

  // Step 3: Invert (I - VG) using LU decomposition
  gsl_permutation *perm = gsl_permutation_alloc(n);
  if (!perm) {
    gsl_matrix_complex_free(VG);
    gsl_matrix_complex_free(I_minus_VG);
    return -1;
  }

  int signum;
  if (gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) != GSL_SUCCESS) {
    gsl_matrix_complex_free(VG);
    gsl_matrix_complex_free(I_minus_VG);
    gsl_permutation_free(perm);
    return -1;
  }

  gsl_matrix_complex *inv_I_minus_VG = gsl_matrix_complex_alloc(n, n);
  if (!inv_I_minus_VG) {
    gsl_matrix_complex_free(VG);
    gsl_matrix_complex_free(I_minus_VG);
    gsl_permutation_free(perm);
    return -1;
  }

  if (gsl_linalg_complex_LU_invert(I_minus_VG, perm, inv_I_minus_VG) !=
      GSL_SUCCESS) {
    gsl_matrix_complex_free(VG);
    gsl_matrix_complex_free(I_minus_VG);
    gsl_matrix_complex_free(inv_I_minus_VG);
    gsl_permutation_free(perm);
    return -1;
  }

  // Step 4: Compute T = inv(I - VG) * V
  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)n, (int)n, (int)n,
              &GSL_REAL(alpha), inv_I_minus_VG->data, (int)inv_I_minus_VG->tda,
              self->V->data, (int)self->V->tda, &GSL_REAL(beta), self->T->data,
              (int)self->T->tda);

  // Clean up
  gsl_matrix_complex_free(VG);
  gsl_matrix_complex_free(I_minus_VG);
  gsl_matrix_complex_free(inv_I_minus_VG);
  gsl_permutation_free(perm);

  return 0;
}

// Run the LSE solver
int lse_compute(LSE *app, double E) {
  if (lse_gmat(app, E) != 0)
    return -1;
  if (lse_vmat(app, E) != 0)
    return -1;
  if (lse_tmat(app, E) != 0)
    return -1;
  return 0;
}

// Compute with new parameters
// Get pointers to matrix data
double *lse_get_g_data(LSE *app) { return (double *)app->G->data; }

void lse_get_g_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->G->size1;
  *cols = (unsigned int)app->G->size2;
}

double *lse_get_v_data(LSE *app) { return (double *)app->V->data; }

void lse_get_v_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->V->size1;
  *cols = (unsigned int)app->V->size2;
}

double *lse_get_t_data(LSE *app) { return (double *)app->T->data; }

void lse_get_t_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->T->size1;
  *cols = (unsigned int)app->T->size2;
}

// Physics calculation functions
// double Omega(int alpha, int beta, double p, double pprime) {
//   if (alpha == 0) {
//     if (beta == 0) {
//       return 2 * m_B + (p * p + pprime * pprime) / 2 / m_B;
//     } else if (beta == 1) {
//       return m_B + pprime * pprime / 2 / m_B + m_B_s + p * p / 2 / m_B_s;
//     }
//   } else if (alpha == 1) {
//     if (beta == 0) {
//       return m_B_s + pprime * pprime / 2 / m_B_s + m_B + p * p / 2 / m_B;
//     } else if (beta == 1) {
//       return 2 * m_B_s + (p * p + pprime * pprime) / 2 / m_B_s;
//     }
//   }
//
//   // Default case (should not happen)
//   fprintf(stderr, "Invalid alpha/beta in Omega\n");
//   return 0;
// }

// double Omega_prime(int alpha, int beta, double p, double pprime) {
//   if (alpha == 0) {
//     if (beta == 0) {
//       return 2 * m_B_star + (p * p + pprime * pprime) / 2 / m_B_star;
//     } else if (beta == 1) {
//       return m_B_star + pprime * pprime / 2 / m_B_star + m_B_star_s +
//              p * p / 2 / m_B_star_s;
//     }
//   } else if (alpha == 1) {
//     if (beta == 0) {
//       return m_B_star_s + pprime * pprime / 2 / m_B_star_s + m_B_star +
//              p * p / 2 / m_B_star;
//     } else if (beta == 1) {
//       return 2 * m_B_star_s + (p * p + pprime * pprime) / 2 / m_B_star_s;
//     }
//   }
//
//   // Default case (should not happen)
//   fprintf(stderr, "Invalid alpha/beta in Omega_prime\n");
//   return 0;
// }

// gsl_complex O(int alpha, int beta, double E, double p, double pprime, double
// m) {
//     double omega = Omega(alpha, beta, p, pprime);
//     double omega_prime = Omega_prime(alpha, beta, p, pprime);
//     double tmp[2] = {
//         E - (m + square(p - pprime) / 2 / m),
//         E - (m + square(p + pprime) / 2 / m)
//     };
//
//     return gsl_complex_mul_real(add(gsl_complex_log(comp((tmp[0] - omega) /
//     (tmp[1] - omega), 0)), gsl_complex_log(comp((tmp[0] - omega_prime) /
//     (tmp[1] - omega_prime),0))),-1 / p / pprime / 4 );
// }

// double V_OME(int alpha, int beta, double E, double p, double pprime) {
//   if (alpha == 0) {
//     if (beta == 0) {
//       return -3 * (3 * O(0, 0, E, p, pprime, m_pi) +
//                    O(0, 0, E, p, pprime, m_eta) / 3);
//     } else if (beta == 1) {
//       return pow(2, 3.0 / 2.0) * O(0, 1, E, p, pprime, m_K);
//     }
//   } else if (alpha == 1) {
//     if (beta == 0) {
//       return pow(2, 3.0 / 2.0) * O(1, 0, E, p, pprime, m_K);
//     } else if (beta == 1) {
//       return 2.0 / 3.0 * O(1, 1, E, p, pprime, m_eta);
//     }
//   }
//
//   // Default case (should not happen)
//   fprintf(stderr, "Invalid alpha/beta in V_OME\n");
//   return 0;
// }

// Main function for testing
int main() {
  // Test the LSE solver
  LSE *lse = lse_new(2, 4, 1e-7);
  if (!lse) {
    fprintf(stderr, "Failed to create LSE solver\n");
    return 1;
  }
  lse_gmat(lse, 0.2);
  uint64_t dim = 2 * (lse->Ngauss + 1);
  for (uint64_t i = 0; i < dim; i += 1) {
    for (uint64_t j = 0; j < dim; j += 1) {
      uint64_t pos = i * dim + j;
      printf("%.2e + Im%.2e   ", lse->G->data[2 * pos],
             lse->G->data[2 * pos + 1]);
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
  lse_deinit(lse);

  return 0;
}
