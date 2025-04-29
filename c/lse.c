#include "lse.h"
#include "wavefunction.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <stdint.h>
#include <string.h>

double complex V_QM(LSE *self, double E, uint64_t p, uint64_t pprime) {
  double complex res = 0 + 0 * I;
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    res += self->psi_n_mat[i][p] * self->psi_n_mat[i][pprime] /
           (E - self->E_vec[i]);
  }
  return res * g_qm;
}

LSE *lse_malloc(size_t pNgauss, double Lambda, double epsilon) {
  LSE *self = (LSE *)malloc(sizeof(LSE));
  if (!self)
    return NULL;

  const size_t n = 2 * (pNgauss + 1);

  self->wf = WFnew(1, RLAMBDA, RNGAUSS);
  self->Ngauss = pNgauss;
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

  self->table = gsl_integration_glfixed_table_alloc(pNgauss);
  self->xi = (double *)malloc(pNgauss * sizeof(double));
  self->wi = (double *)malloc(pNgauss * sizeof(double));

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

  for (size_t i = 0; i < pNgauss; i++) {
    gsl_integration_glfixed_point(0, Lambda, i, &self->xi[i], &self->wi[i],
                                  self->table);
  }

  self->psi_n_mat = (double complex **)malloc(sizeof(double complex *) * N_MAX);
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    self->psi_n_mat[i] =
        (double complex *)malloc(sizeof(double complex) * (pNgauss + 3));
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
void lse_refresh(LSE *self, double E) {
  // self->Lambda = Lambda;
  self->E = E;

  for (int i = 0; i < 2; i++) {
    const double dE = E - delta[i];
    const double mU = mu[i];
    const double complex tmp = csqrt(2 * mU * dE);
    self->x0[i] = tmp;
    for (uint64_t j = 0; j < N_MAX; j += 1) {
      self->psi_n_mat[j][self->Ngauss + i] = psi_n_ftcomplex(self->wf, tmp, j + 1);
    }
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
  WFfree(self->wf);
  free(self->xi);
  free(self->wi);
  free(self->E_vec);
  for (uint64_t i = 0; i < N_MAX; i += 1) {
    free(self->psi_n_mat[i]);
  }
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
      tmp = tmp - int_val * x0 * x0;

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
int lse_vmat(LSE *self) {
  gsl_matrix_complex_set_zero(self->V);

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
      gsl_matrix_complex_set(self->V, i00, j00,
                           V_OME_00(self->E, p, pprime) + Ctct_00(g_c));
      
      gsl_matrix_complex_set(self->V, i01, j01,
                           V_OME_01(self->E, p, pprime) + Ctct_01(g_c));
      
      gsl_matrix_complex_set(self->V, i10, j10,
                           V_OME_10(self->E, p, pprime) + Ctct_01(g_c));
      
      gsl_matrix_complex_set(self->V, i11, j11,
                           V_OME_11(self->E, p, pprime) + Ctct_11(g_c));
    }
  }
  
  // Handle edge cases separately
  
  // Case 1: idx = Ngauss (special x0 value for first dimension)
  size_t idx = self->Ngauss;
  double complex p_edge = self->x0[0]; // Use x0 for the edge case
  
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
    gsl_matrix_complex_set(self->V, i00, j00,
                         V_OME_00(self->E, p_edge, pprime) + Ctct_00(g_c));
    
    gsl_matrix_complex_set(self->V, i01, j01,
                         V_OME_01(self->E, p_edge, pprime) + Ctct_01(g_c));
    
    p_edge = self->x0[1];
    gsl_matrix_complex_set(self->V, i10, j10,
                         V_OME_10(self->E, p_edge, pprime) + Ctct_01(g_c));
    
    gsl_matrix_complex_set(self->V, i11, j11,
                         V_OME_11(self->E, p_edge, pprime) + Ctct_11(g_c));
  }
  
  // Case 2: jdx = Ngauss (special x0 value for second dimension)
  size_t jdx = self->Ngauss;
  double complex pprime_edge = self->x0[0]; // Use x0 for the edge case
  
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
    gsl_matrix_complex_set(self->V, i00, j00,
                         V_OME_00(self->E, p, pprime_edge) + Ctct_00(g_c));
    
    pprime_edge = self->x0[1]; // Use x0[1] for the second dimension
    gsl_matrix_complex_set(self->V, i01, j01,
                         V_OME_01(self->E, p, pprime_edge) + Ctct_01(g_c));
    
    pprime_edge = self->x0[0]; // Use x0[0] for the first dimension
    gsl_matrix_complex_set(self->V, i10, j10,
                         V_OME_10(self->E, p, pprime_edge) + Ctct_01(g_c));
    
    pprime_edge = self->x0[1]; // Use x0[1] for the second dimension
    gsl_matrix_complex_set(self->V, i11, j11,
                         V_OME_11(self->E, p, pprime_edge) + Ctct_11(g_c));
  }
  
  // Case 3: Both idx = Ngauss and jdx = Ngauss (special x0 values for both dimensions)
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
  gsl_matrix_complex_set(self->V, i00, j00,
                       V_OME_00(self->E, self->x0[0], self->x0[0]) + Ctct_00(g_c));
  
  gsl_matrix_complex_set(self->V, i01, j01,
                       V_OME_01(self->E, self->x0[0], self->x0[1]) + Ctct_01(g_c));
  
  gsl_matrix_complex_set(self->V, i10, j10,
                       V_OME_10(self->E, self->x0[1], self->x0[0]) + Ctct_01(g_c));
  
  gsl_matrix_complex_set(self->V, i11, j11,
                       V_OME_11(self->E, self->x0[1], self->x0[1]) + Ctct_11(g_c));


  return 0;
}

// Calculate T matrix
int lse_tmat(LSE *self) {
  const size_t n = 2 * (self->Ngauss + 1);

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
  lse_refresh(app, E);
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
