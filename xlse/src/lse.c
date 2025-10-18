#include "lse.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_vector_double.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "autofree.h"
#include "constants.h"
#include "ome.h"
#include "wavefunction.h"

constexpr int max_iter = 4000;
#define IPVG
#define TESTQM
// double complex V_QM(LSE *self, size_t p, size_t pprime) {
//   double complex res = 0 + 0 * I;
//   __auto_type E = self->E;
//   for (size_t i = 0; i < N_MAX; i += 1) {
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
    LSE *self = malloc(sizeof(LSE));
    if (!self)
        return NULL;

    const size_t n = NCHANNELS * (pNgauss + 1);

    self->ome.set = 0;
    ome_build(&self->ome);
    self->wf = WFnew(partialwave, RLAMBDA, RNGAUSS);
    self->pNgauss = pNgauss;
    self->Lambda = Lambda;
    self->epsilon = epsilon;
    // self->E = E;

    self->TOME = matrix_alloc(n, n);
    self->VOME = matrix_alloc(n, n);
    self->G = matrix_alloc(n, n);
    self->reg = matrix_alloc(n, n);

    if (!self->TOME || !self->VOME || !self->G) {
        if (self->TOME)
            matrix_free(self->TOME);
        if (self->VOME)
            matrix_free(self->VOME);
        if (self->G)
            matrix_free(self->G);
        free(self);
        return NULL;
    }

    self->table = gsl_integration_glfixed_table_alloc(pNgauss);
    self->xi = malloc(pNgauss * sizeof(double));
    self->wi = malloc(pNgauss * sizeof(double));

    if (!self->table || !self->xi || !self->wi) {
        if (self->TOME)
            matrix_free(self->TOME);
        if (self->VOME)
            matrix_free(self->VOME);
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
    self->psi_n_mat =
        malloc(sizeof(double complex) * NCHANNELS * (N_MAX + 1) * (pNgauss + 1));
    double complex(*psi)[N_MAX + 1][pNgauss + 1] = self->psi_n_mat;
#ifdef TESTQM
    for (size_t i = 0; i < N_MAX; i += 1) {
        for (size_t pi = 0; pi < pNgauss; pi += 1) {
            psi[0][i][pi] = psi_test(self->xi[pi]);
        }
        memcpy(psi[1][i], psi[0][i], sizeof(double complex) * pNgauss);
    }
#elifdef CONSTQM
    for (size_t i = 0; i < N_MAX; i += 1) {
        for (size_t pi = 0; pi < pNgauss; pi += 1) {
            psi[0][i][pi] = 1;
        }
        memcpy(psi[1][i], psi[0][i], sizeof(double complex) * pNgauss);
    }
#else
    for (size_t i = 0; i < N_MAX; i += 1) {
        psi_n_ft_batch(self->wf, self->xi, psi[0][i], pNgauss, i + 1);
        memcpy(psi[1][i], psi[0][i], sizeof(double complex) * pNgauss);
    }
#endif
    for (size_t ch = 0; ch < NCHANNELS; ch += 1) {
        for (size_t pi = 0; pi < pNgauss + 1; pi += 1) {
            psi[ch][N_MAX][pi] = 1;
        }
    }
    self->E_vec = malloc(sizeof(double) * N_MAX);
    memcpy(self->E_vec, self->wf->E_solution->data, N_MAX * sizeof(double));

    size_t ntp = N_MAX + 1;
    self->Xin = malloc(sizeof(double complex) * 4 * ntp * (pNgauss + 1));
    self->Xout = malloc(sizeof(double complex) * 4 * ntp * (pNgauss + 1));
    self->v = malloc(sizeof(double complex) * 4 * ntp);

    self->sigmat =
        malloc(4 * (N_MAX + 1) * (N_MAX + 1) * sizeof(double complex));

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
void lse_refresh(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    // self->Lambda = Lambda;
    self->E = E;
    self->V0 = C[0];
    self->C00 = C[1];
    self->C01 = C[2];
    self->C10 = C[2];
    self->C11 = C[3];
    self->g0 = g[0];
    self->g1 = g[1];

    double complex(*psi)[N_MAX + 1][self->pNgauss + 1] = self->psi_n_mat;
#ifdef TPO

    double complex(*v)[2][N_MAX + 1] = self->v;
    v[0][0][N_MAX] = C[1];
    v[0][1][N_MAX] = C[2];
    v[1][0][N_MAX] = C[2];
    v[1][1][N_MAX] = C[3];
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t i = 0; i < N_MAX; i += 1) {
                v[alpha][beta][i] =
                    gab[alpha] * gab[beta] / (E - self->E_vec[i] - self->V0);
            }
        }
    }
#endif
    for (int i = 0; i < NCHANNELS; i++) {
        const double complex dE = E - delta[i];
        const double mU = mu[i];
        const double complex tmp = csqrt(2 * mU * dE);
        if (((rs >> i) & 1) == 1) {
            self->x0[i] = tmp;
        } else {
            self->x0[i] = -(tmp);
        }
        for (size_t j = 0; j < N_MAX; j += 1) {
#ifdef TESTQM
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] = psi_test(self->x0[i]);
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
#elifdef CONSTQM
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] = 1;
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
#else
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] =
                    psi_n_ftcomplex(self->wf, self->x0[i], j + 1);
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
            // psi[i][j][self->pNgauss] = psi_n_ftcomplex(self->wf, self->x0[i], j + 1);
#endif
        }
    }
}

void lse_refreshm(LSE *self, double complex p, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS]) {
    // self->Lambda = Lambda;
    self->E = p * p / 2 / mu[1];
    self->V0 = C[0];
    self->C00 = C[1];
    self->C01 = C[2];
    self->C10 = C[2];
    self->C11 = C[3];
    self->g0 = g[0];
    self->g1 = g[1];

    double complex(*psi)[N_MAX + 1][self->pNgauss + 1] = self->psi_n_mat;
#ifdef TPO

    double complex(*v)[2][N_MAX + 1] = self->v;
    v[0][0][N_MAX] = C[1];
    v[0][1][N_MAX] = C[2];
    v[1][0][N_MAX] = C[2];
    v[1][1][N_MAX] = C[3];
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t i = 0; i < N_MAX; i += 1) {
                v[alpha][beta][i] =
                    gab[alpha] * gab[beta] / (E - self->E_vec[i] - self->V0);
            }
        }
    }
#endif /* ifdef TPO                                                              \
                                                                               \ \
   */
    for (int i = 0; i < 1; i++) {
        self->x0[i] = p;
        for (size_t j = 0; j < N_MAX; j += 1) {
#ifdef TESTQM
            auto dE = self->E;
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] = psi_test(self->x0[i]);
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
                psi[i][j][self->pNgauss] = psi_test(self->x0[i]);
#elifdef CONSTQM
            auto dE = self->E;
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] = 1;
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
            // psi[i][j][self->pNgauss] = 1;
#else
            auto dE = self->E;
            if (creal(dE) > 0 && fabs(cimag(dE)) < 1e-7) {
                psi[i][j][self->pNgauss] =
                    psi_n_ftcomplex(self->wf, self->x0[i], j + 1);
            } else {
                psi[i][j][self->pNgauss] = 0;
            }
#endif
        }
    }
}

// Free LSE resources
void lse_free(LSE *self) {
    if (!self)
        return;

    matrix_free(self->TOME);
    matrix_free(self->VOME);
    matrix_free(self->G);
    matrix_free(self->reg);
    gsl_integration_glfixed_table_free(self->table);
    WFfree(self->wf);
    free(self->xi);
    free(self->wi);
    free(self->psi_n_mat);
    free(self->E_vec);
    free(self->Xin);
    free(self->Xout);
    free(self->v);
    free(self->sigmat);
    free(self);
}

// Calculate G matrix
int lse_gmat(LSE *self) {
    matrix_set_zero(self->G);

    for (size_t i = 0; i < NCHANNELS; i++) {
        const double complex dE = self->E - delta[i];
        const double mU = mu[i];
        const double complex x0 = self->x0[i];
        double complex int_val = 0 + 0 * I;
        for (size_t j = 0; j < self->pNgauss; j++) {
            double x = self->xi[j];
            double w = self->wi[j];

            double complex denom = dE - fsquare(x) / 2 / mU + self->epsilon * I;
            int_val += w / denom;
        }

        double complex tmp =
            mU * x0 *
            (clog((self->Lambda + x0) / (self->Lambda - x0)) - M_PI * I);
        // printf("%.4f\n",
        //        creal((clog((self->Lambda + x0) / (self->Lambda - x0))
        //        - M_PI * I)));
        // printf("%.4f\n",
        //        cimag((clog((self->Lambda + x0) / (self->Lambda - x0))
        //        - M_PI * I)));
        tmp = tmp - int_val * x0 * x0;
#ifdef IPVG
        tmp = -tmp;
#endif

        const size_t ii = self->pNgauss + i * (self->pNgauss + 1);
        matrix_set(self->G, ii, ii, tmp * (1 / fsquare(M_PI) / 2));

        for (size_t m = 0; m < self->pNgauss; m++) {
            const size_t pos = m + i * (self->pNgauss + 1);
            double complex denominator =
                dE - fsquare(self->xi[m]) / 2 / mU + self->epsilon * I;
#ifdef IPVG
            double complex ele = -fsquare(self->xi[m]) * self->wi[m] / 2 /
                                 fsquare(M_PI) / denominator;
#elifdef IMVG
            double complex ele = fsquare(self->xi[m]) * self->wi[m] / 2 /
                                 fsquare(M_PI) / denominator;
#endif
#ifdef DEBUG
            if (m == 12 && i == 0) {
                puts("\n---------G---------\n");
                auto p = self->xi[m];
                auto step = self->wi[m];
                printf("p: %s\n", formatC(p));
                printf("step: %s\n", formatC(step));
                printf("E-Delta: %s\n", formatC(dE));
                printf("mu: %s\n", formatC(mU));
                printf("denominator: %s\n", formatC(denominator));
                printf("G: %s\n", formatC(ele));
                puts("-------\n");
            }
#endif /* ifdef DEBUG */

            matrix_set(self->G, pos, pos, ele);
        }
    }

    return 0;
}

// Calculate V matrix
int lse_vmat(LSE *self) {
    size_t pNgauss = self->pNgauss;
    matrix_set_zero(self->VOME);
#define matrindex(idx, block) ((idx) + (block) * (pNgauss + 1))
#define setv(idx, jdx, ib, jb) \
    matrix_set(self->VOME, matrindex(idx, ib), matrindex(jdx, jb), V##ib##jb(self, self->xi[idx], self->xi[jdx], idx, jdx))
#define setvonr(idx, jdx, ib, jb) \
    matrix_set(self->VOME, matrindex(idx, ib), matrindex(jdx, jb), V##ib##jb(self, self->xi[idx], self->x0[jb], idx, matrindex(jdx, jb)))
#define setvonl(idx, jdx, ib, jb) \
    matrix_set(self->VOME, matrindex(idx, ib), matrindex(jdx, jb), V##ib##jb(self, self->x0[ib], self->xi[jdx], matrindex(idx, ib), jdx))
#define setvonf(idx, jdx, ib, jb) \
    matrix_set(self->VOME, matrindex(idx, ib), matrindex(jdx, jb), V##ib##jb(self, self->x0[ib], self->x0[jb], matrindex(idx, ib), matrindex(jdx, jb)))
    // gsl_matrix_complex_set_identity(self->VOME);
    // return 0;

    for (size_t idx = 0; idx < self->pNgauss; idx++) {
        for (size_t jdx = 0; jdx < self->pNgauss; jdx++) {
            setv(idx, jdx, 0, 0);
            setv(idx, jdx, 0, 1);
            setv(idx, jdx, 1, 0);
            setv(idx, jdx, 1, 1);
        }
    }

    // Handle edge cases separately

    // Case 1: idx = Ngauss (special x0 value for first dimension)
    size_t idx = self->pNgauss;

    for (size_t jdx = 0; jdx < self->pNgauss; jdx++) {
        setvonl(idx, jdx, 0, 0);
        setvonl(idx, jdx, 0, 1);
        setvonl(idx, jdx, 1, 0);
        setvonl(idx, jdx, 1, 1);
    }

    // Case 2: jdx = Ngauss (special x0 value for second dimension)
    size_t jdx = self->pNgauss;

    for (size_t idx = 0; idx < self->pNgauss; idx++) {
        setvonr(idx, jdx, 0, 0);
        setvonr(idx, jdx, 0, 1);
        setvonr(idx, jdx, 1, 0);
        setvonr(idx, jdx, 1, 1);
    }

    // Case 3: Both idx = Ngauss and jdx = Ngauss (special x0 values for
    // both dimensions)
    idx = self->pNgauss;
    jdx = self->pNgauss;

    setvonf(idx, jdx, 0, 0);
    setvonf(idx, jdx, 0, 1);
    setvonf(idx, jdx, 1, 0);
    setvonf(idx, jdx, 1, 1);

#ifdef DEBUG
#endif /* ifdef DEBUG */
#undef matrindex
#undef setv
#undef setvonr
#undef setvonl
#undef setvonf

    return 0;
}

// Calculate T matrix
//
int lse_tmat(LSE *self) {
    const size_t n = NCHANNELS * (self->pNgauss + 1);
    gsl_matrix_complex *IVG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!IVG)
        return -1;
    gsl_matrix_complex_set_identity(IVG);

#ifdef IMVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, self->VOME, self->G, 1, IVG);
#elif defined(IPVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, self->VOME, self->G, 1, IVG);
#endif

    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(IVG, perm, &signum) !=
        GSL_SUCCESS) {
        return -1;
    }

    self->det = gsl_linalg_complex_LU_det(IVG, signum);
    gsl_matrix_complex *inv_I_minus_VG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!inv_I_minus_VG) {
        return -1;
    }

    if (gsl_linalg_complex_LU_invert(IVG, perm, inv_I_minus_VG) !=
        GSL_SUCCESS) {
        return -1;
    }
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, inv_I_minus_VG,
                   self->VOME, 0, self->TOME);

    // Clean up
    return 0;
}

int lse_tmat_single(LSE *self) {
    const size_t n = self->pNgauss + 1;
    auto V = gsl_matrix_complex_submatrix(self->VOME, 0, 0, n, n);
    auto G = gsl_matrix_complex_submatrix(self->G, 0, 0, n, n);
    auto T = gsl_matrix_complex_submatrix(self->TOME, 0, 0, n, n);
    gsl_matrix_complex *VG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    // gsl_matrix_complex *I_minus_VG [[gnu::cleanup(matfree)]] =
    // gsl_matrix_complex_alloc(n, n);
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    gsl_matrix_complex *inv_I_minus_VG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!inv_I_minus_VG)
        return -1;
    if (!perm)
        return -1;
    if (!VG)
        return -1;
    //    if (!I_minus_VG)
    // return -1;
    gsl_matrix_complex_set_zero(VG);
    // gsl_matrix_complex_set_zero(I_minus_VG);

    gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta = gsl_complex_rect(0.0, 0.0);
#ifdef IMVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, &V.matrix, &G.matrix, 1, VG);
#elif defined(IPVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, &V.matrix, &G.matrix, 1, VG);

#endif
    int signum;
    if (gsl_linalg_complex_LU_decomp(VG, perm, &signum) != GSL_SUCCESS) {
        return -1;
    }

    self->det = gsl_linalg_complex_LU_det(VG, signum);
    if (gsl_linalg_complex_LU_invert(VG, perm, inv_I_minus_VG) != GSL_SUCCESS) {
        return -1;
    }
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, inv_I_minus_VG, &V.matrix,
                   beta, &T.matrix);
    self->onshellT[0][0] = matrix_get(&T.matrix, n - 1, n - 1);
    return GSL_SUCCESS;
}

double complex lse_invT(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_compute_single(self, E, C, g, rs);
    size_t n = self->pNgauss;
    auto T = gsl_matrix_complex_submatrix(self->TOME, 0, 0, n, n);
    gsl_matrix_complex *tmp [[gnu::cleanup(matfree)]] = matrix_alloc(n, n);
    gsl_matrix_complex *inv [[gnu::cleanup(matfree)]] = matrix_alloc(n, n);
    gsl_matrix_complex_memcpy(tmp, &T.matrix);
    int signum;
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (gsl_linalg_complex_LU_decomp(tmp, perm, &signum) != GSL_SUCCESS)
        exit(1);
    if (gsl_linalg_complex_LU_invert(tmp, perm, inv) != GSL_SUCCESS)
        exit(1);
    auto res = matrix_get(inv, n, n);
    return res;
}

void lse_X(LSE *self) {
    size_t pNgauss = self->pNgauss;
    double complex(*psi)[N_MAX + 1][pNgauss + 1] = self->psi_n_mat;
    double complex(*Xin)[2][N_MAX + 1][pNgauss + 1] = self->Xin;
    double complex(*Xout)[2][N_MAX + 1][pNgauss + 1] = self->Xout;
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            auto xin = Xin[alpha][beta];
            auto xout = Xout[alpha][beta];
            for (size_t i = 0; i < N_MAX + 1; i += 1) {
                for (size_t pi = 0; pi < pNgauss + 1; pi += 1) {
                    xin[i][pi] = 0;
                    xout[i][pi] = 0;
                    for (size_t ki = 0; ki < pNgauss + 1; ki += 1) {
                        double complex G =
                            matrix_get(self->G, ki + alpha * (pNgauss + 1),
                                       ki + alpha * (pNgauss + 1));
                        double complex T =
                            matrix_get(self->TOME, ki + alpha * (pNgauss + 1),
                                       pi + beta * (pNgauss + 1));
                        xin[i][pi] += psi[alpha][i][ki] * G * T;
                        G = matrix_get(self->G, ki + beta * (pNgauss + 1),
                                       ki + beta * (pNgauss + 1));
                        T = matrix_get(self->TOME, pi + alpha * (pNgauss + 1),
                                       ki + beta * (pNgauss + 1));
                        xout[i][pi] += T * G * psi[beta][i][ki];
                    }
                    xin[i][pi] += psi[beta][i][pi] * Delta[alpha][beta];
                    xout[i][pi] += psi[alpha][i][pi] * Delta[alpha][beta];
                }
            }
        }
    }
    double complex(*sigma)[2][N_MAX + 1][N_MAX + 1] = self->sigmat;
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t i = 0; i < N_MAX + 1; i += 1) {
                for (size_t j = 0; j < N_MAX + 1; j += 1) {
                    sigma[alpha][beta][i][j] = 0;
                    for (size_t ki = 0; ki < pNgauss + 1; ki += 1) {
                        double complex G =
                            matrix_get(self->G, ki + alpha * (pNgauss + 1),
                                       ki + alpha * (pNgauss + 1));
                        sigma[alpha][beta][i][j] +=
                            psi[alpha][i][ki] * G * Xout[alpha][beta][j][ki];
                    }
                }
            }
        }
    }
}

void lse_XtX(LSE *self) {
    constexpr size_t n = N_MAX + 1;
    size_t pNgauss = self->pNgauss;
    double complex(*sigma)[2][n][n] = self->sigmat;
    double complex(*Xout)[2][n][pNgauss + 1] = self->Xout;
    [[maybe_unused]] double complex(*Xin)[2][n][pNgauss + 1] = self->Xin;
    matrix *M [[gnu::cleanup(matfree)]] = matrix_alloc(2 * n, 2 * n);
    matrix_set_zero(M);
    double complex(*mat)[2 * n] = (double complex(*)[2 * n]) M->data;
    double complex(*v)[2][n] = self->v;
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t gamma = 0; gamma < 2; gamma += 1) {
                for (size_t i = 0; i < n; i += 1) {
                    for (size_t j = 0; j < n; j += 1) {
                        mat[i + alpha * n][j + beta * n] -=
                            v[alpha][gamma][i] * sigma[gamma][beta][i][j];
                    }
                }
            }
        }
    }
    gsl_matrix_complex_add_diagonal(M, 1);
    matrix *invM [[gnu::cleanup(matfree)]] = matrix_alloc(2 * n, 2 * n);
    gsl_permutation *perm [[gnu::cleanup(permfree)]] =
        gsl_permutation_alloc(2 * n);
    int signum;
    gsl_linalg_complex_LU_decomp(M, perm, &signum);
    gsl_linalg_complex_LU_invert(M, perm, invM);
    double complex(*invmat)[2 * n] = (double complex(*)[2 * n]) invM->data;
    double complex(*invMv)[2][n][n] = malloc(2 * sizeof(*invMv));
    memset(invMv, 0, 2 * sizeof(*invMv));
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t gamma = 0; gamma < 2; gamma += 1) {
                for (size_t i = 0; i < N_MAX + 1; i += 1) {
                    for (size_t j = 0; j < N_MAX + 1; j += 1) {
                        invMv[alpha][beta][i][j] +=
                            invmat[i + alpha * n][j + gamma * n] *
                            v[gamma][beta][j];
                    }
                }
            }
        }
    }
    auto T = self->onshellT;
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            T[alpha][beta] =
                matrix_get(self->TOME, pNgauss + alpha * (pNgauss + 1),
                           pNgauss + beta * (pNgauss + 1));
        }
    }
    for (size_t alpha = 0; alpha < 2; alpha += 1) {
        for (size_t beta = 0; beta < 2; beta += 1) {
            for (size_t gamma = 0; gamma < 2; gamma += 1) {
                for (size_t delta = 0; delta < 2; delta += 1) {
                    for (size_t i = 0; i < N_MAX + 1; i += 1) {
                        for (size_t j = 0; j < N_MAX + 1; j += 1) {
                            T[alpha][beta] += Xout[alpha][gamma][i][pNgauss] *
                                              invMv[gamma][delta][i][j] *
                                              (Xout[delta][beta][j][pNgauss]);
                        }
                    }
                }
            }
        }
    }
    free(invMv);
}

double complex *lse_get_iivg_data(LSE *self) {
    lse_gmat(self);
    lse_vmat(self);
    const size_t n = 2 * (self->pNgauss + 1);
    // Step 1: Compute VG = V * G
    gsl_matrix_complex *VG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!VG)
        exit(-1);
    gsl_matrix_complex_set_zero(VG);

    // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
    gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta = gsl_complex_rect(0.0, 0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, self->VOME, self->G, beta,
                   VG);

    // Step 2: Compute I - VG
    gsl_matrix_complex *I_minus_VG [[gnu::cleanup(matfree)]] =
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
    //   gsl_matrix_complex_set(I_minus_VG, i, i, gsl_complex_add(diag,
    //   one));
    // }

    // Step 3: Invert (I - VG) using LU decomposition
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        exit(-1);
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(I_minus_VG, perm, &signum) !=
        GSL_SUCCESS) {
        exit(-1);
    }

    gsl_matrix_complex *inv_I_minus_VG [[gnu::cleanup(matfree)]] =
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

double complex lse_detImVG(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_refresh(self, E, C, g, rs);
    lse_gmat(self);
    lse_vmat(self);
    const size_t n = NCHANNELS * (self->pNgauss + 1);

    // Step 1: Compute VG = V * G
    gsl_matrix_complex *IVG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!IVG)
        return -1;
    gsl_matrix_complex_set_identity(IVG);

#ifdef IMVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, self->VOME, self->G, 1, IVG);
#elif defined(IPVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, self->VOME, self->G, 1, IVG);
#endif

    // Add identity matrix: I_minus_VG = I - VG
    gsl_matrix_complex_add_diagonal(IVG, 1 + 0I);
    // gsl_matrix_complex_memcpy(self->reg, I_minus_VG);

    // Step 3: Invert (I - VG) using LU decomposition
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(IVG, perm, &signum) !=
        GSL_SUCCESS) {
        return -1;
    }
    // printf("det: %.9f\n", cabs(gsl_linalg_complex_LU_det(I_minus_VG,
    // signum)));
    return gsl_linalg_complex_LU_det(IVG, signum);
}

double complex lse_detImVGsing(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_refresh(self, E, C, g, rs);
    lse_gmat(self);
    lse_vmat(self);
    // const size_t n = 2 * (self->pNgauss + 1);
    size_t n;
    if (creal(self->E) > 0 && fabs(cimag(self->E)) < 1e-5) {
        n = self->pNgauss + 1;
    } else {
        n = self->pNgauss;
    }
    auto V = gsl_matrix_complex_submatrix(self->VOME, 0, 0, n, n);
    auto G = gsl_matrix_complex_submatrix(self->G, 0, 0, n, n);

    // Step 1: Compute VG = V * G
    gsl_matrix_complex *IVG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!IVG)
        return -1;
    gsl_matrix_complex_set_zero(IVG);

    // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
    gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta = gsl_complex_rect(0.0, 0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, &V.matrix, &G.matrix, 1, IVG);
#ifdef IMVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, &V.matrix, &G.matrix, 1, IVG);
#elif defined(IPVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, &V.matrix, &G.matrix, 1, IVG);
#endif

    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(IVG, perm, &signum) !=
        GSL_SUCCESS) {
        return -1;
    }
    // printf("det: %.9f\n", cabs(gsl_linalg_complex_LU_det(I_minus_VG,
    // signum)));
    return gsl_linalg_complex_LU_det(IVG, signum);
}

double complex lse_detImVGsingm(LSE *self, double complex p, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS]) {
    lse_refreshm(self, p, C, g);
    lse_gmat(self);
    lse_vmat(self);
    // const size_t n = 2 * (self->pNgauss + 1);
    const size_t n = self->pNgauss + 1;
    auto V = gsl_matrix_complex_submatrix(self->VOME, 0, 0, n, n);
    auto G = gsl_matrix_complex_submatrix(self->G, 0, 0, n, n);

    gsl_matrix_complex *IVG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!IVG)
        return -1;
    gsl_matrix_complex_set_identity(IVG);

#ifdef IMVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, &V.matrix, &G.matrix, 1, IVG);
#elif defined(IPVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, &V.matrix, &G.matrix, 1, IVG);
#endif

    // Step 3: Invert (I - VG) using LU decomposition
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(IVG, perm, &signum) !=
        GSL_SUCCESS) {
        return -1;
    }
    // printf("det: %.9f\n", cabs(gsl_linalg_complex_LU_det(I_minus_VG,
    // signum)));
    return gsl_linalg_complex_LU_det(IVG, signum);
}

double complex lse_detIJm(LSE *self, _Complex double p, const double *C, const double *g) {
    lse_refreshm(self, p, C, g);
    lse_gmat(self);
    constexpr double E[6] = {m_Xb11P, m_Xb12P, m_Xb13P, m_Xb14P, m_Xb15P, m_Xb16P};
    auto n = self->pNgauss + 1;
    matrix *tau [[gnu::cleanup(matfree)]] = matrix_alloc(6, 6);
	matrix_set_zero(tau);
    auto dE = p * p;
    for (uint64_t i = 0; i < 6; i += 1) {
        matrix_set(tau, i, i, g[0] * g[0] / (dE - E[i]));
    }
    [[maybe_unused]] double complex(*psi)[N_MAX + 1][n] = self->psi_n_mat;
    matrix *J [[gnu::cleanup(matfree)]] = matrix_alloc(6, 6);
    for (uint64_t i = 0; i < 6; i += 1) {
        for (uint64_t j = 0; j < 6; j += 1) {
            double complex ele = 0;
            for (uint64_t pi = 0; pi < n; pi += 1) {
                auto G = matrix_get(self->G, pi, pi);
                ele += G;
            }
            matrix_set(J, i, j, ele);
        }
    }
    matrix *IJ [[gnu::cleanup(matfree)]] = matrix_alloc(6, 6);
    gsl_matrix_complex_set_identity(IJ);
#ifdef IPVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, tau, J, 1, IJ);
#elif defined(IMVG)
        gsl_blas)zgemm(CblasNoTrans, CblasNoTrans, -1, tau, J, 1, IJ);
#endif
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(6);
    if (!perm)
        return -1;

    int signum;
    if (gsl_linalg_complex_LU_decomp(IJ, perm, &signum) != GSL_SUCCESS)
        return -1;
    return gsl_linalg_complex_LU_det(IJ, signum);
}

double complex lse_detVG(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_refresh(self, E, C, g, rs);
    lse_gmat(self);
    lse_vmat(self);
    const size_t n = 2 * (self->pNgauss + 1);

    // Step 1: Compute VG = V * G
    gsl_matrix_complex *VG [[gnu::cleanup(matfree)]] =
        gsl_matrix_complex_alloc(n, n);
    if (!VG)
        return -1;
    gsl_matrix_complex_set_zero(VG);

    // Use GSL BLAS wrapper for matrix multiplication: VG = V * G
    gsl_complex alpha = gsl_complex_rect(1.0, 0.0);
    gsl_complex beta = gsl_complex_rect(0.0, 0.0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, self->VOME, self->G, beta,
                   VG);

    // Step 3: Invert (I - VG) using LU decomposition
    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(n);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(VG, perm, &signum) != GSL_SUCCESS) {
        return -1;
    }
    return gsl_linalg_complex_LU_det(VG, signum);
}

double lse_cost(LSE *self, const double C[NCHANNELS * NCHANNELS], RS rs) {
    auto g = (double[2]){self->g0, self->g1};
    double res = cabs(lse_detImVG(self, m_Xb11P, C, g, rs));
    res += cabs(lse_detImVG(self, m_Xb12P, C, g, rs));
    res += cabs(lse_detImVG(self, m_Xb13P, C, g, rs));
    // res += cabs(lse_detImVG(self, m_Xb14P, C, g,rs));
    res += fabs(creal(lse_detImVG(self, m_Xb14P, C, g, rs)));
    return res;
}

double lse_costsing(LSE *self, const double C[NCHANNELS * NCHANNELS], RS rs) {
    // for (uint64_t i = 0; i < 4; i += 1) {
    // 	printf("C%zu %f\n", i, C[i]);
    // }
    auto g = (double[2]){self->g0, self->g1};
    double res = cabs(lse_detImVGsing(self, m_Xb11P, C, g, rs));
    res += cabs(lse_detImVGsing(self, m_Xb12P, C, g, rs));
    res += cabs(lse_detImVGsing(self, m_Xb13P, C, g, rs));
    res += fabs(creal(lse_detImVGsing(self, m_Xb14P, C, g, rs)));
    return res;
}

// Run the LSE solver
int lse_compute(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_refresh(self, E, C, g, rs);
    if (lse_gmat(self) != 0)
        return -1;
    if (lse_vmat(self) != 0)
        return -1;
    if (lse_tmat(self) != 0)
        return -1;
    return 0;
}

int lse_compute_single(LSE *self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    lse_refresh(self, E, C, g, rs);
    if (lse_gmat(self) != 0)
        return -1;
    if (lse_vmat(self) != 0)
        return -1;
    if (lse_tmat_single(self) != 0)
        return -1;
    return 0;
}

double complex *lse_get_g_data(LSE *self) {
    return (double complex *)self->G->data;
}

void lse_get_g_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)self->G->size1;
    *cols = (unsigned int)self->G->size2;
}

double complex *lse_get_v_data(LSE *self) {
    return (double complex *)self->VOME->data;
}

void lse_get_v_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)self->VOME->size1;
    *cols = (unsigned int)self->VOME->size2;
}

double complex *lse_get_t_data(LSE *self) {
    return (double complex *)self->TOME->data;
}

void lse_get_t_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)self->TOME->size1;
    *cols = (unsigned int)self->TOME->size2;
}

double complex *lse_get_ivg_data(LSE *self) {
    return (double complex *)self->reg->data;
}

void lse_get_ivg_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)self->reg->size1;
    *cols = (unsigned int)self->reg->size2;
}

void lse_get_iivg_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)self->VOME->size1;
    *cols = (unsigned int)self->VOME->size2;
}

void lse_get_M_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)(2 * (self->pNgauss + 1));
    *cols = (unsigned int)(2 * (self->pNgauss + 1));
}

double complex *lse_get_onshellT(LSE *self) {
    return (double complex *)self->onshellT;
};

double complex *lse_get_psi(LSE *self) { return self->psi_n_mat; }

void lse_get_psi_size(LSE *self, unsigned int *rows, unsigned int *cols) {
    *rows = (unsigned int)N_MAX;
    *cols = (unsigned int)(self->pNgauss + 1);
}

double *lse_get_E(LSE *self) { return self->E_vec; }
void lse_get_E_size(unsigned int *levels) { *levels = (unsigned int)N_MAX; }

struct detparams {
    void *lse;
    double C[NCHANNELS * NCHANNELS];
    double g[NCHANNELS];
    RS rs;
};
int detImVG(const gsl_vector *x, void *params, gsl_vector *f) {
    struct detparams *detp = params;
    double complex E = gsl_vector_get(x, 0) + gsl_vector_get(x, 1) * I;
    auto det = lse_detImVG(detp->lse, E, detp->C, detp->g, detp->rs);
    gsl_vector_set(f, 0, creal(det));
    gsl_vector_set(f, 1, cimag(det));
    return GSL_SUCCESS;
}

int detIJ(const gsl_vector *x, void *params, gsl_vector *f) {
    constexpr double E[6] = {m_Xb11P, m_Xb12P, m_Xb13P, m_Xb14P, m_Xb15P, m_Xb16P};
    double complex p = gsl_vector_get(x, 0) + gsl_vector_get(x, 1) * I;
    struct detparams *detp = params;
    LSE *lse = detp->lse;
    const auto n = lse->pNgauss + 1;
    matrix *tau [[gnu::cleanup(matfree)]] = gsl_matrix_complex_alloc(6, 6);
    auto g = detp->g[0];
    gsl_matrix_complex_set_zero(tau);
    auto dE = p * p/2/mu[1];
    for (uint64_t i = 0; i < 6; i += 1) {
        matrix_set(tau, i, i, g * g / (dE - E[i]));
    }
    lse_refreshm(lse, p, detp->C, detp->g);
    lse_gmat(lse);
    double complex(*psi)[N_MAX + 1][n] = lse->psi_n_mat;
    matrix *J [[gnu::cleanup(matfree)]] = matrix_alloc(6, 6);
    for (uint64_t i = 0; i < 6; i += 1) {
        for (uint64_t j = 0; j < 6; j += 1) {
            double complex ele = 0;
            for (uint64_t pi = 0; pi < n; pi += 1) {
                auto G = matrix_get(lse->G, pi, pi);
                ele += conj(psi[0][i][pi])*G*psi[0][j][pi];
            }
            matrix_set(J, i, j, ele);
        }
    }
    matrix *IJ [[gnu::cleanup(matfree)]] = matrix_alloc(6, 6);
    gsl_matrix_complex_set_identity(IJ);
#ifdef IPVG
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, 1, tau, J, 1, IJ);
#elif defined(IMVG)
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, -1, tau, J, 1, IJ);
#endif

    gsl_permutation *perm [[gnu::cleanup(permfree)]] = gsl_permutation_alloc(6);
    if (!perm) {
        return -1;
    }

    int signum;
    if (gsl_linalg_complex_LU_decomp(IJ, perm, &signum) !=
        GSL_SUCCESS) {
        return -1;
    }
    auto det = gsl_linalg_complex_LU_det(IJ, signum);
    gsl_vector_set(f, 0, creal(det));
    gsl_vector_set(f, 1, cimag(det));
    return GSL_SUCCESS;
}

int detImVGm(const gsl_vector *x, void *params, gsl_vector *f) {
    struct detparams *detp = params;
    double complex p = gsl_vector_get(x, 0) + gsl_vector_get(x, 1) * I;
    auto det = lse_detImVGsingm(detp->lse, p, detp->C, detp->g);
    gsl_vector_set(f, 0, creal(det));
    gsl_vector_set(f, 1, cimag(det));
    return GSL_SUCCESS;
}

double complex pole(LSE *lse, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid;
    struct detparams detp = {
        .lse = lse,
        .C = {C[0], C[1], C[2], C[3]},
        .g = {g[0], g[1]},
        .rs = rs,
    };
    gsl_multiroot_function F = {&detImVG, 2, &detp};
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, creal(E));
    gsl_vector_set(x, 1, cimag(E));
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_multiroot_fsolver_set(s, &F, x);
    int status;
    uint iter = 0;
    do {
        iter += 1;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status)
            break;
        status = gsl_multiroot_test_residual(s->f, 1e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    double complex res = NAN * NAN * I;
    if (status == GSL_SUCCESS) {
        double re = gsl_vector_get(s->x, 0);
        double im = gsl_vector_get(s->x, 1);
        if (fabs(re) < 2 && fabs(im) < 1.7) {
            res = re + im * I;
        }
    }
    gsl_vector_free(x);
    gsl_multiroot_fsolver_free(s);
    return res;
}

double complex polem(LSE *lse, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs) {
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid;
    struct detparams detp = {
        .lse = lse,
        .C = {C[0], C[1], C[2], C[3]},
        .g = {g[0], g[1]},
        .rs = rs,
    };
    gsl_multiroot_function F = {&detIJ, 2, &detp};
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector_set(x, 0, creal(E));
    gsl_vector_set(x, 1, cimag(E));
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_multiroot_fsolver_set(s, &F, x);
    int status;
    uint iter = 0;
    do {
        iter += 1;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status)
            break;
        status = gsl_multiroot_test_residual(s->f, 1e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    double complex res = NAN * NAN * I;
    if (status == GSL_SUCCESS) {
        double re = gsl_vector_get(s->x, 0);
        double im = gsl_vector_get(s->x, 1);
        if (fabs(re) < 2 && fabs(im) < 2.) {
            res = re + im * I;
        }
    }
    gsl_vector_free(x);
    gsl_multiroot_fsolver_free(s);
    return res;
}

double cost(const gsl_vector *x, void *params) {
    LSE *lse = params;
    double C[NCHANNELS * NCHANNELS];
    for (size_t i = 0; i < 4; i += 1) {
        C[i] = gsl_vector_get(x, i);
    }
    // C[0] = -1.6 + (-0.9 + 1.6)/(1 + exp(-C[0]));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb11P, C, PP)));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb12P, C, PP)));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb13P, C, PP)));
    return lse_cost(lse, C, PP);
}

double costsing(const gsl_vector *x, void *params) {
    LSE *lse = params;
    double C[NCHANNELS * NCHANNELS];
    for (size_t i = 0; i < 4; i += 1) {
        C[i] = gsl_vector_get(x, i);
    }
    // C[0] = -1.6 + (-0.9 + 1.6)/(1 + exp(-C[0]));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb11P, C, PP)));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb12P, C, PP)));
    // res = max(res, cabs(lse_detImVG(lse, m_Xb13P, C, PP)));
    return lse_costsing(lse, C, PP);
}

void minimize(LSE *lse, const double Cin[NCHANNELS * NCHANNELS], double Cout[NCHANNELS * NCHANNELS + 1]) {
    gsl_vector *x = gsl_vector_alloc(4);
    for (size_t i = 0; i < NCHANNELS * NCHANNELS; i += 1) {
        gsl_vector_set(x, i, Cin[i]);
    }
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_function func;
    func.n = 4;
    func.f = cost;
    func.params = lse;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 4);
    gsl_vector *step = gsl_vector_alloc(4);
    gsl_vector_set_all(step, 0.8);
    gsl_multimin_fminimizer_set(s, &func, x, step);
    int status;
    int iter = 0;
    int iterpercent = max_iter / 20;
    double tolerance = 1e-4;
    double size;
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status) {
            break;
        }
        if (iter % iterpercent == 0) {
            printf("%d-th iteration\n", iter);
            printf("fcn: %12.7e\n", s->fval);
            // printf("C[0]: %f\n\n", carray[0]);
        }
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tolerance);
    } while (status == GSL_CONTINUE && iter < max_iter);
    if (status == GSL_SUCCESS) {
        puts("minimization succeeded");
    } else {
        puts("minimization failed");
    }
    printf("iteration: %d\n", iter);
    printf("fcn: %12.6e\n", s->fval);
    for (size_t i = 0; i < 4; i += 1) {
        Cout[i] = gsl_vector_get(s->x, i);
    }
    Cout[4] = s->fval;
    // res[0] = -1.6 + (-0.9 + 1.6)/(1 + exp(-res[0]));
    puts("res:");
    for (size_t i = 0; i < 4; i += 1) {
        printf("  %12.6e\n", Cout[i]);
    }
    gsl_vector_free(x);
    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(s);
}

void minimizesing(LSE *lse, const double Cin[NCHANNELS * NCHANNELS], double Cout[NCHANNELS * NCHANNELS + 1]) {
    gsl_vector *x = gsl_vector_alloc(4);
    for (size_t i = 0; i < NCHANNELS * NCHANNELS; i += 1) {
        gsl_vector_set(x, i, Cin[i]);
    }
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_function func;
    func.n = 4;
    func.f = costsing;
    func.params = lse;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 4);
    gsl_vector *step = gsl_vector_alloc(4);
    gsl_vector_set_all(step, 0.8);
    gsl_multimin_fminimizer_set(s, &func, x, step);
    int status;
    int iter = 0;
    int max_iter = 6000;
    double tolerance = 1e-4;
    int iterpercent = max_iter / 10;
    double size;
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status) {
            break;
        }
        if (iter % iterpercent == 0) {
            printf("%d-th iteration\n", iter);
            printf("fcn: %12.7e\n", s->fval);
        }
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tolerance);
    } while (status == GSL_CONTINUE && iter < max_iter);
    if (status == GSL_SUCCESS) {
        puts("minimization succeeded");
    } else {
        puts("minimization failed");
    }
    printf("iteration: %d\n", iter);
    printf("fcn: %12.6e\n", s->fval);
    for (size_t i = 0; i < 4; i += 1) {
        Cout[i] = gsl_vector_get(s->x, i);
    }
    Cout[4] = s->fval;
    // res[0] = -1.6 + (-0.9 + 1.6)/(1 + exp(-res[0]));
    puts("res:");
    for (size_t i = 0; i < 4; i += 1) {
        printf("  %12.6e\n", Cout[i]);
    }
    gsl_vector_free(x);
    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(s);
}

// complex double V(size_t alpha, size_t beta, double E, double complex p,
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
