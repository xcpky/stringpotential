#include "lse.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

LSE *lse_malloc(size_t Ngauss, double Lambda, double epsilon) {
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
void lse_free(LSE *self) {
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

      double complex int_val = 0 + 0 * I;

      for (size_t j = 0; j < self->Ngauss; j++) {
        double x = self->xi[j];
        double w = self->wi[j];

        double complex green = dE - square(x) / 2 / mU + self->epsilon * I;
        int_val += w / green;
      }

      double complex tmp =
          mU * x0 * log((self->Lambda + x0) / (self->Lambda - x0)) -
          mU * x0 * M_PI * I;
      tmp = tmp - int_val * square(x0);

      const size_t ii = self->Ngauss + i * (self->Ngauss + 1);
      matrix_set(self->G, ii, ii, tmp * (1 / square(M_PI) / 2));
    }

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
      gsl_matrix_complex_set(self->V, i, j, V_OME_00(self->E, p[0], pprime));
      pprime = x[1][jdx];
      i = idx + 0 * (self->Ngauss + 1);
      j = jdx + 1 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_01(self->E, p[1], pprime));
      pprime = x[0][jdx];
      i = idx + 1 * (self->Ngauss + 1);
      j = jdx + 0 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_10(self->E, p[2], pprime));
      pprime = x[1][jdx];
      i = idx + 1 * (self->Ngauss + 1);
      j = jdx + 1 * (self->Ngauss + 1);
      gsl_matrix_complex_set(self->V, i, j, V_OME_11(self->E, p[3], pprime));
    }
  }
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
