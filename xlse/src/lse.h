#ifndef LSE_H
#define LSE_H
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>

#include "constants.h"
#include "ome.h"
#include "wavefunction.h"

#define NCHANNELS 2

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

// #define TPO

#define vlarray(ptr, ...) ((double complex(*) __VA_ARGS__)(ptr))

typedef struct {
    size_t pNgauss;
    double Lambda;
    double epsilon;
    double complex E;
    double complex det;
    matrix* TOME;
    matrix* VOME;
    matrix* G;
    matrix* reg;
    double complex x0[NCHANNELS];
    gsl_integration_glfixed_table* table;
    WaveFunction* wf;
    double* xi;
    double* wi;
    struct OME ome;
    void* psi_n_mat;
    double* E_vec;
    void* Xin;
    void* Xout;
    double complex onshellT[NCHANNELS][NCHANNELS];
    void* v;
    void* sigmat;
    double V0;
    double C00;
    double C01;
    double C10;
    double C11;
    double g0;
    double g1;
} LSE;

typedef enum : uint64_t {
    NN = 0,
    PN = 1,
    NP = 2,
    PP = 3,
} RS;

LSE* lse_malloc(size_t pNgauss, double Lambda, double epsilon);
int lse_compute(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
int lse_compute_single(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
void lse_free(LSE* self);
double complex* lse_get_g_data(LSE* self);
void lse_get_g_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_v_data(LSE* self);
void lse_get_v_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_t_data(LSE* self);
void lse_get_t_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_ivg_data(LSE* self);
void lse_get_ivg_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_iivg_data(LSE* self);
void lse_get_iivg_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_psi(LSE* self);
void lse_get_psi_size(LSE* self, unsigned int* rows, unsigned int* cols);
double* lse_get_E(LSE* self);
void lse_get_E_size(unsigned int* levels);
void lse_get_M_size(LSE* self, unsigned int* rows, unsigned int* cols);
double complex* lse_get_onshellT(LSE* self);
double complex pole(LSE* lse, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double complex polem(LSE* lse, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double cost(const gsl_vector* x, void* params);
double costsing(const gsl_vector* x, void* params);
void minimize(LSE* lse, const double Cin[NCHANNELS * NCHANNELS], double Cout[NCHANNELS * NCHANNELS + 1]);
void minimizesing(LSE* lse, const double Cin[NCHANNELS * NCHANNELS], double Cout[NCHANNELS * NCHANNELS + 1]);

// LSE methods
int lse_gmat(LSE* self);
int lse_vmat(LSE* self);
int lse_tmat(LSE* self);
int lse_tmat_single(LSE* self);
double complex lse_invT(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double complex lse_detImVG(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double complex lse_detImVGsing(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double complex lse_detImVGsingm(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS]);
double complex lse_detVG(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
double lse_cost(LSE* self, const double C[NCHANNELS * NCHANNELS], RS rs);
double lse_costsing(LSE* self, const double C[NCHANNELS * NCHANNELS], RS rs);
void lse_refresh(LSE* self, double complex E, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS], RS rs);
void lse_refreshm(LSE* self, double complex p, const double C[NCHANNELS * NCHANNELS], const double g[NCHANNELS]);
void lse_X(LSE* self);
void lse_XtX(LSE* self);

static inline double min(double x, double y) { return x > y ? y : x; }

static inline double max(double x, double y) { return x > y ? x : y; }

double complex V_QM_00(LSE* self, size_t p, size_t pprime);
double complex V_QM_01(LSE* self, size_t p, size_t pprime);
double complex V_QM_10(LSE* self, size_t p, size_t pprime);
double complex V_QM_11(LSE* self, size_t p, size_t pprime);

#define DEFINE_VQMTEST(alpha, beta)                                                                 \
    static inline double complex V_QM_TEST_##alpha##beta(LSE* self, [[maybe_unused]] size_t p,      \
                                                         [[maybe_unused]] size_t pprime) {          \
        double E[6] = {m_Xb11P, m_Xb12P, m_Xb13P, m_Xb14P, m_Xb15P, m_Xb16P};                       \
        double complex res = 0;                                                                     \
        size_t pNgauss = self->pNgauss;                                                             \
        auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) self->psi_n_mat;                     \
        size_t chan0 = p / (pNgauss + 1);                                                           \
        size_t chan1 = pprime / (pNgauss + 1);                                                      \
        for (size_t i = 0; i < 6; i += 1) {                                                         \
            res += psi[chan0][i][p % (pNgauss + 1)] * conj(psi[chan1][i][pprime % (pNgauss + 1)]) / \
                   (self->E - E[i]);                                                                \
        }                                                                                           \
        [[maybe_unused]] auto g0 = self->g0;                                                        \
        [[maybe_unused]] auto g1 = self->g1;                                                        \
        return res * g##alpha * g##beta;                                                            \
    }

DEFINE_VQMTEST(0, 0)
DEFINE_VQMTEST(0, 1)
DEFINE_VQMTEST(1, 0)
DEFINE_VQMTEST(1, 1)

#define DEFINE_VQM(alpha, beta)                                                                     \
    double complex V_QM_##alpha##beta(LSE* self, size_t p, size_t pprime) {                         \
        double complex res = 0 + 0I;                                                                \
        auto E = self->E;                                                                           \
        size_t pNgauss = self->pNgauss;                                                             \
        auto psi = (double complex(*)[N_MAX + 1][pNgauss + 1]) self->psi_n_mat;                     \
        size_t chan0 = p / (pNgauss + 1);                                                           \
        size_t chan1 = pprime / (pNgauss + 1);                                                      \
        for (size_t i = 0; i < N_TOWER; i++) {                                                      \
            res += psi[chan0][i][p % (pNgauss + 1)] * conj(psi[chan1][i][pprime % (pNgauss + 1)]) / \
                   (E - self->E_vec[i] - self->V0);                                                 \
        }                                                                                           \
        [[maybe_unused]] auto g0 = self->g0;                                                        \
        [[maybe_unused]] auto g1 = self->g1;                                                        \
        return res * g##alpha * g##beta;                                                            \
    }

#define DEFINE_V_FUNCTION(suffix)                                                                        \
    static inline gsl_complex V##suffix(LSE* self, [[maybe_unused]] double complex p1,                   \
                                        [[maybe_unused]] double complex p2, [[maybe_unused]] size_t p1i, \
                                        [[maybe_unused]] size_t p2i) {                                   \
        [[maybe_unused]] auto E = self->E;                                                               \
        E += m11 + m12;                                                                                  \
        [[maybe_unused]] auto res = -V_QM_TEST_##suffix(self, p1i, p2i);                                 \
        return res;                                                                                      \
    }

DEFINE_V_FUNCTION(00);
DEFINE_V_FUNCTION(01);
DEFINE_V_FUNCTION(10);
DEFINE_V_FUNCTION(11);

#endif  // !LSE_H
