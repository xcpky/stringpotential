#include "lse.h"
#include "autofree.h"
#include "wavefunction.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_permutation.h>
#include <stdint.h>
#include <string.h>

// double complex V_QM(LSE *self, uint64_t p, uint64_t pprime) {
//   double complex res = 0 + 0 * I;
//   __auto_type E = self->E;
//   for (uint64_t i = 0; i < N_MAX; i += 1) {
//     // res += self->psi_n_mat[i][p] * conj(self->psi_n_mat[i][pprime]) /
//     //        (E - self->E_vec[i]);
//     res -= 1 / (E - self->E_vec[i] + self->epsilon*I);
//   }
//   return res * g1*g2;
// }
//
DEFINE_VQM(0, 0)
DEFINE_VQM(0, 1)
DEFINE_VQM(1, 0)
DEFINE_VQM(1, 1)

LSE *lse_malloc(size_t pNgauss, double Lambda, double epsilon) {
  LSE *self = (LSE *)malloc(sizeof(LSE));
  if (!self)
    return NULL;

  const size_t n = 2 * (pNgauss + 1);

  self->wf = WFnew(partialwave, RLAMBDA, RNGAUSS);
  self->Ngauss = pNgauss;
  self->Lambda = Lambda;
  self->epsilon = epsilon;
  // self->E = E;

  self->T = matrix_alloc(n, n);
  self->V = matrix_alloc(n, n);
  self->G = matrix_alloc(n, n);
  self->reg = matrix_alloc(n, n);

  if (!self->T || !self->V || !self->G) {
    if (self->T)
      matrix_free(self->T);
    if (self->V)
      matrix_free(self->V);
    if (self->G)
      matrix_free(self->G);
    free(self);
    return NULL;
  }

  self->table = gsl_integration_glfixed_table_alloc(pNgauss);
  self->xi = (double *)malloc(pNgauss * sizeof(double));
  self->wi = (double *)malloc(pNgauss * sizeof(double));

  if (!self->table || !self->xi || !self->wi) {
    if (self->T)
      matrix_free(self->T);
    if (self->V)
      matrix_free(self->V);
    if (self->G)
      matrix_free(self->G);
    if (self->table)
      gsl_integration_glfixed_table_free(self->table);
    if (self->xi)
      free(self->xi);
    if (self->wi)
      free(self->wi);
    free(self);
    return NULL;
  }

  for (size_t i = 0; i < pNgauss; i++) {
    gsl_integration_glfixed_point(0, Lambda, i, &self->xi[i], &self->wi[i],
                                  self->table);
  }

  self->psi_raw =
      (double complex *)malloc(sizeof(double complex) * N_MAX * (pNgauss + 2));
  self->psi_n_mat = (double complex **)malloc(sizeof(double complex *) * N_MAX);
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    self->psi_n_mat[i] = &self->psi_raw[i * (pNgauss + 2)];
    psi_n_ft_batch(self->wf, self->xi, self->psi_n_mat[i], pNgauss, i + 1);
  }
  self->E_vec = (double *)malloc(sizeof(double) * N_MAX);
  memcpy(self->E_vec, self->wf->E_solution->data, N_MAX * sizeof(double));
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
void lse_refresh(LSE *self, double complex E, int64_t rs) {
  // self->Lambda = Lambda;
  self->E = E;

  for (int i = 0; i < 2; i++) {
    const double complex dE = E - delta[i];
    const double mU = mu[i];
    const double complex tmp = xsqrt(2 * mU * dE);
    self->x0[i] = rs * tmp;
    for (uint64_t j = 0; j < N_MAX; j += 1) {
      self->psi_n_mat[j][self->Ngauss + i] =
          psi_n_ftcomplex(self->wf, rs * tmp, j + 1);
    }
  }
}

// Free LSE resources
void lse_free(LSE *self) {
  if (!self)
    return;

  matrix_free(self->T);
  matrix_free(self->V);
  matrix_free(self->G);
  gsl_integration_glfixed_table_free(self->table);
  WFfree(self->wf);
  free(self->xi);
  free(self->wi);
  free(self->E_vec);
  free(self->psi_raw);
  free(self->psi_n_mat);
  free(self);
}

// Calculate G matrix
int lse_gmat(LSE *self) {
  matrix_set_zero(self->G);

  for (int i = 0; i < 2; i++) {
    const double dE = self->E - delta[i];
    const double mU = mu[i];
    const double complex x0 = self->x0[i];
    double complex int_val = 0 + 0 * I;
    for (size_t j = 0; j < self->Ngauss; j++) {
      double x = self->xi[j];
      double w = self->wi[j];

      double complex denom = dE - square(x) / 2 / mU + self->epsilon * I;
      int_val += w / denom;
    }

    double complex tmp =
        mU * x0 * clog((self->Lambda + x0) / (self->Lambda - x0)) -
        mU * x0 * M_PI * I;
    tmp = tmp - int_val * x0 * x0;

    const size_t ii = self->Ngauss + i * (self->Ngauss + 1);
    matrix_set(self->G, ii, ii, tmp * (1 / square(M_PI) / 2));

    for (size_t m = 0; m < self->Ngauss; m++) {
      const size_t pos = m + i * (self->Ngauss + 1);
      double complex denominator =
          dE - square(self->xi[m]) / 2 / mU + self->epsilon * I;
      double complex ele =
          square(self->xi[m]) * self->wi[m] / 2 / square(M_PI) / denominator;

      matrix_set(self->G, pos, pos, ele);
    }
  }

  return 0;
}

// Calculate V matrix
int lse_vmat(LSE *self) {
  matrix_set_zero(self->V);

  for (size_t idx = 0; idx < self->Ngauss; idx++) {
    double complex p = self->xi[idx];
    for (size_t jdx = 0; jdx < self->Ngauss; jdx++) {
      double complex pprime = self->xi[jdx];

      // Calculate matrix indices
      size_t i00 = idx + 0 * (self->Ngauss + 1);
      size_t j00 = jdx + 0 * (self->Ngauss + 1);
      size_t i01 = idx + 0 * (self->Ngauss + 1);
      size_t j01 = jdx + 1 * (self->Ngauss + 1);
      size_t i10 = idx + 1 * (self->Ngauss + 1);
      size_t j10 = jdx + 0 * (self->Ngauss + 1);
      size_t i11 = idx + 1 * (self->Ngauss + 1);
      size_t j11 = jdx + 1 * (self->Ngauss + 1);

      // Set matrix elements for regular cases
      matrix_set(self->V, i00, j00, V00(self, p, pprime));

      matrix_set(self->V, i01, j01, V01(self, p, pprime));

      matrix_set(self->V, i10, j10, V10(self, p, pprime));

      matrix_set(self->V, i11, j11, V11(self, p, pprime));
    }
  }

  // Handle edge cases separately

  // Case 1: idx = Ngauss (special x0 value for first dimension)
  size_t idx = self->Ngauss;

  for (size_t jdx = 0; jdx < self->Ngauss; jdx++) {
    double complex pprime = self->xi[jdx];

    // Calculate matrix indices
    size_t i00 = idx + 0 * (self->Ngauss + 1);
    size_t j00 = jdx + 0 * (self->Ngauss + 1);
    size_t i01 = idx + 0 * (self->Ngauss + 1);
    size_t j01 = jdx + 1 * (self->Ngauss + 1);
    size_t i10 = idx + 1 * (self->Ngauss + 1);
    size_t j10 = jdx + 0 * (self->Ngauss + 1);
    size_t i11 = idx + 1 * (self->Ngauss + 1);
    size_t j11 = jdx + 1 * (self->Ngauss + 1);

    // Set matrix elements for edge cases
    matrix_set(self->V, i00, j00, V00(self, self->x0[0], pprime));

    matrix_set(self->V, i01, j01, V01(self, self->x0[0], pprime));

    matrix_set(self->V, i10, j10, V10(self, self->x0[1], pprime));

    matrix_set(self->V, i11, j11, V11(self, self->x0[1], pprime));
  }

  // Case 2: jdx = Ngauss (special x0 value for second dimension)
  size_t jdx = self->Ngauss;

  for (size_t idx = 0; idx < self->Ngauss; idx++) {
    double complex p = self->xi[idx];

    // Calculate matrix indices
    size_t i00 = idx + 0 * (self->Ngauss + 1);
    size_t j00 = jdx + 0 * (self->Ngauss + 1);
    size_t i01 = idx + 0 * (self->Ngauss + 1);
    size_t j01 = jdx + 1 * (self->Ngauss + 1);
    size_t i10 = idx + 1 * (self->Ngauss + 1);
    size_t j10 = jdx + 0 * (self->Ngauss + 1);
    size_t i11 = idx + 1 * (self->Ngauss + 1);
    size_t j11 = jdx + 1 * (self->Ngauss + 1);

    // Set matrix elements for edge cases
    matrix_set(self->V, i00, j00, V00(self, p, self->x0[0]));
    matrix_set(self->V, i01, j01, V01(self, p, self->x0[1]));

    matrix_set(self->V, i10, j10, V10(self, p, self->x0[0]));

    matrix_set(self->V, i11, j11, V11(self, p, self->x0[1]));
  }

  // Case 3: Both idx = Ngauss and jdx = Ngauss (special x0 values for both
  // dimensions)
  idx = self->Ngauss;
  jdx = self->Ngauss;

  // Calculate matrix indices for the corner case
  size_t i00 = idx + 0 * (self->Ngauss + 1);
  size_t j00 = jdx + 0 * (self->Ngauss + 1);
  size_t i01 = idx + 0 * (self->Ngauss + 1);
  size_t j01 = jdx + 1 * (self->Ngauss + 1);
  size_t i10 = idx + 1 * (self->Ngauss + 1);
  size_t j10 = jdx + 0 * (self->Ngauss + 1);
  size_t i11 = idx + 1 * (self->Ngauss + 1);
  size_t j11 = jdx + 1 * (self->Ngauss + 1);

  // Set matrix elements for the corner case
  matrix_set(self->V, i00, j00, V00(self, self->x0[0], self->x0[0]));

  matrix_set(self->V, i01, j01, V01(self, self->x0[0], self->x0[1]));

  matrix_set(self->V, i10, j10, V10(self, self->x0[1], self->x0[0]));

  matrix_set(self->V, i11, j11, V11(self, self->x0[1], self->x0[1]));
#ifdef DEBUG
  __auto_type val = V00(self, self->x0[0], self->x0[0]);
  printf("(%f) + Im(%f)\n", creal(val), cimag(val));
#endif /* ifdef DEBUG */

  return 0;
}

// Calculate T matrix
int lse_tmat(LSE *self) {
  const size_t n = 2 * (self->Ngauss + 1);

  // Step 1: Compute VG = V * G
  gsl_matrix_complex *VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!VG)
    return -1;
  gsl_matrix_complex_set_zero(VG);

  // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
  gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
  gsl_complex beta = gsl_complex_rect(0.0, 0.0);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, self->V, self->G, beta, VG);

  // Step 2: Compute I - VG
  gsl_matrix_complex *I_minus_VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!I_minus_VG) {
    return -1;
  }

  gsl_matrix_complex_memcpy(I_minus_VG, VG);
  gsl_matrix_complex_scale(I_minus_VG,
                           gsl_complex_rect(-1.0, 0.0)); // I_minus_VG = -VG

  // Add identity matrix: I_minus_VG = I - VG
  gsl_matrix_complex_add_diagonal(I_minus_VG, 1 + 0I);
  gsl_matrix_complex_memcpy(self->reg, I_minus_VG);
  // for (size_t i = 0; i < n; i++) {
  //   gsl_complex diag = gsl_matrix_complex_get(I_minus_VG, i, i);
  //   gsl_complex one = gsl_complex_rect(1.0, 0.0);
  //   gsl_matrix_complex_set(I_minus_VG, i, i, gsl_complex_add(diag, one));
  // }

  // Step 3: Invert (I - VG) using LU decomposition
  gsl_permutation *perm __attribute__((cleanup(auto_permfree))) =
      gsl_permutation_alloc(n);
  if (!perm) {
    return -1;
  }

  int signum;
  if (gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) != GSL_SUCCESS) {
    return -1;
  }

  gsl_matrix_complex *inv_I_minus_VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!inv_I_minus_VG) {
    return -1;
  }

  if (gsl_linalg_complex_LU_invert(I_minus_VG, perm, inv_I_minus_VG) !=
      GSL_SUCCESS) {
    return -1;
  }
  // Step 4: Compute T = inv(I - VG) * V using GSL BLAS wrapper
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, inv_I_minus_VG, self->V,
                 beta, self->T);

  // Clean up
  return 0;
}

double complex *lse_get_iivg_data(LSE *self) {
  lse_gmat(self);
  lse_vmat(self);
  const size_t n = 2 * (self->Ngauss + 1);
  // Step 1: Compute VG = V * G
  gsl_matrix_complex *VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!VG)
    exit(-1);
  gsl_matrix_complex_set_zero(VG);

  // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
  gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
  gsl_complex beta = gsl_complex_rect(0.0, 0.0);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, self->V, self->G, beta, VG);

  // Step 2: Compute I - VG
  gsl_matrix_complex *I_minus_VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!I_minus_VG) {
    exit(-1);
  }

  gsl_matrix_complex_memcpy(I_minus_VG, VG);
  gsl_matrix_complex_scale(I_minus_VG,
                           gsl_complex_rect(-1.0, 0.0)); // I_minus_VG = -VG

  // Add identity matrix: I_minus_VG = I - VG
  gsl_matrix_complex_add_diagonal(I_minus_VG, 1 + 0I);
  // gsl_matrix_complex_memcpy(self->reg, I_minus_VG);
  // for (size_t i = 0; i < n; i++) {
  //   gsl_complex diag = gsl_matrix_complex_get(I_minus_VG, i, i);
  //   gsl_complex one = gsl_complex_rect(1.0, 0.0);
  //   gsl_matrix_complex_set(I_minus_VG, i, i, gsl_complex_add(diag, one));
  // }

  // Step 3: Invert (I - VG) using LU decomposition
  gsl_permutation *perm __attribute__((cleanup(auto_permfree))) =
      gsl_permutation_alloc(n);
  if (!perm) {
    exit(-1);
  }

  int signum;
  if (gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) != GSL_SUCCESS) {
    exit(-1);
  }

  gsl_matrix_complex *inv_I_minus_VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!inv_I_minus_VG) {
    exit(-1);
  }

  if (gsl_linalg_complex_LU_invert(I_minus_VG, perm, inv_I_minus_VG) !=
      GSL_SUCCESS) {
    exit(-1);
  }
  gsl_matrix_complex_memcpy(self->reg, inv_I_minus_VG);
  return (double complex *)self->reg->data;
}

double complex lse_detImVG(LSE *self, double complex E) {
  lse_refresh(self, E, 1);
  lse_gmat(self);
  lse_vmat(self);
  const size_t n = 2 * (self->Ngauss + 1);

  // Step 1: Compute VG = V * G
  gsl_matrix_complex *VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!VG)
    return -1;
  gsl_matrix_complex_set_zero(VG);

  // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
  gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
  gsl_complex beta = gsl_complex_rect(0.0, 0.0);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, self->V, self->G, beta, VG);

  // Step 2: Compute I - VG
  gsl_matrix_complex *I_minus_VG __attribute__((cleanup(auto_matfree))) =
      gsl_matrix_complex_alloc(n, n);
  if (!I_minus_VG) {
    return -1;
  }

  gsl_matrix_complex_memcpy(I_minus_VG, VG);
  gsl_matrix_complex_scale(I_minus_VG,
                           gsl_complex_rect(-1.0, 0.0)); // I_minus_VG = -VG

  // Add identity matrix: I_minus_VG = I - VG
  gsl_matrix_complex_add_diagonal(I_minus_VG, 1 + 0I);
  gsl_matrix_complex_memcpy(self->reg, I_minus_VG);

  // Step 3: Invert (I - VG) using LU decomposition
  gsl_permutation *perm __attribute__((cleanup(auto_permfree))) =
      gsl_permutation_alloc(n);
  if (!perm) {
    return -1;
  }

  int signum;
  if (gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) != GSL_SUCCESS) {
    return -1;
  }
  return gsl_linalg_complex_LU_det(I_minus_VG, signum);
}

// Run the LSE solver
int lse_compute(LSE *app, double complex E, int64_t rs) {
  lse_refresh(app, E, rs);
  if (lse_gmat(app) != 0)
    return -1;
  if (lse_vmat(app) != 0)
    return -1;
  if (lse_tmat(app) != 0)
    return -1;
  return 0;
}

double complex *lse_get_g_data(LSE *app) {
  return (double complex *)app->G->data;
}

void lse_get_g_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->G->size1;
  *cols = (unsigned int)app->G->size2;
}

double complex *lse_get_v_data(LSE *app) {
  return (double complex *)app->V->data;
}

void lse_get_v_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->V->size1;
  *cols = (unsigned int)app->V->size2;
}

double complex *lse_get_t_data(LSE *app) {
  return (double complex *)app->T->data;
}

void lse_get_t_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->T->size1;
  *cols = (unsigned int)app->T->size2;
}

double complex *lse_get_ivg_data(LSE *app) {
  return (double complex *)app->reg->data;
}

void lse_get_ivg_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->reg->size1;
  *cols = (unsigned int)app->reg->size2;
}

void lse_get_iivg_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)app->V->size1;
  *cols = (unsigned int)app->V->size2;
}

void lse_get_m_size(LSE *app, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)2 * (app->Ngauss + 1);
  *cols = (unsigned int)2 * (app->Ngauss + 1);
}

double complex *lse_get_psi(LSE *self) { return self->psi_raw; }

void lse_get_psi_size(LSE *self, unsigned int *rows, unsigned int *cols) {
  *rows = (unsigned int)N_MAX;
  *cols = (unsigned int)(self->Ngauss + 2);
}

double *lse_get_E(LSE *self) { return self->E_vec; }
void lse_get_E_size(unsigned int *levels) { *levels = (unsigned int)N_MAX; }

// complex double V(uint64_t alpha, uint64_t beta, double E, double complex p,
//                  double complex pprime) {
//   if (alpha + beta == 0) {
//     return V00(E, p, pprime);
//   } else if (alpha + beta == 2) {
//     return V11(E, p, pprime);
//   } else if (alpha == 1) {
//     return V10(E, p, pprime);
//   } else if (alpha == 0) {
//     return V01(E, p, pprime);
//   }
//   return 0;
// }
//
double complex noninline_O_00(double E, double complex p, double complex pprime,
                              double m) {
  // if (fabs(p) <= 1e-8) {
  //   p += 1e-6;
  // }
  // if (fabs(pprime) <= 1e-8) {
  //   pprime += 1e-6;
  // }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_00(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_00(p, pprime) +
               0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_00(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_00(p, pprime) +
               0 * I));
}

double complex noninline_O_01(double E, double complex p, double complex pprime,
                              double m) {
  // if (fabs(p) <= 1e-8) {
  //   p += 1e-6;
  // }
  // if (fabs(pprime) <= 1e-8) {
  //   pprime += 1e-6;
  // }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_01(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_01(p, pprime) +
               0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_01(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_01(p, pprime) +
               0 * I));
}

double complex noninline_O_10(double E, double complex p, double complex pprime,
                              double m) {
  // if (fabs(p) <= 1e-8) {
  //   p += 1e-6;
  // }
  // if (fabs(pprime) <= 1e-8) {
  //   pprime += 1e-6;
  // }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_10(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_10(p, pprime) +
               0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_10(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_10(p, pprime) +
               0 * I));
}

double complex noninline_O_11(double E, double complex p, double complex pprime,
                              double m) {
  // if (fabs(p) <= 1e-8) {
  //   p += 1e-6;
  // }
  // if (fabs(pprime) <= 1e-8) {
  //   pprime += 1e-6;
  // }
  return -1. / 4. / p / pprime *
         (clog((E - (m + csquare(p - pprime) / 2 / m) - omega_11(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omega_11(p, pprime) +
               0 * I) +
          clog((E - (m + csquare(p - pprime) / 2 / m) -
                omegaprime_11(p, pprime)) /
                   E -
               (m + csquare(p + pprime) / 2 / m) - omegaprime_11(p, pprime) +
               0 * I));
}

double complex noninline_V_OME_00(double E, double complex p,
                                  double complex pprime) {
  return -3 * (3 * noninline_O_00(E, p, pprime, m_pi) +
               noninline_O_00(E, p, pprime, m_eta) / 3);
}

double complex noninline_V_OME_01(double E, double complex p,
                                  double complex pprime) {
  return pow(2, 3. / 2) * noninline_O_01(E, p, pprime, m_K);
}

double complex noninline_V_OME_10(double E, double complex p,
                                  double complex pprime) {
  return pow(2, 3. / 2) * noninline_O_10(E, p, pprime, m_K);
}

double complex noninline_V_OME_11(double E, double complex p,
                                  double complex pprime) {
  return 2. / 3 * noninline_O_11(E, p, pprime, m_eta);
}

double complex O00(double complex E, double complex p, double complex pprime,
                   double m) {
  return O_00(E, p, pprime, m);
}
double complex O01(double complex E, double complex p, double complex pprime,
                   double m) {
  return O_01(E, p, pprime, m);
}
double complex O10(double complex E, double complex p, double complex pprime,
                   double m) {
  return O_10(E, p, pprime, m);
}
double complex O11(double complex E, double complex p, double complex pprime,
                   double m) {
  return O_11(E, p, pprime, m);
}
