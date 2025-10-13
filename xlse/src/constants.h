#ifndef CONSTANTS_H
#include <complex.h>
#define CONSTANTS_H
#define PNGAUSS 64
#define RNGAUSS 256
#define RLAMBDA 20

#define f_pi (0.092)
#define g_b (0.5704)
#define g_c (-1. / 4.)

#define m_pi (0.138039407)
#define m_K (0.498)
#define m_eta (0.548)
#define m_eta_s (0.69)

#define WIDTH (-1e-6 * I)
// #define DEVI (-0.2)
#define m_B (5.27952)
#define m_B_star (5.32475 + WIDTH)
// #define m_B_star (m_B + m_pi + 0.01)
#define m_B_s (5.36691)
#define m_B_star_s (5.4154 + WIDTH)
#define gamma_B_star (1e-6)
#define gamma_B_star_s (1e-6)
#define m_B0 (5.27952)
#define m_B1 (5.36691)
#define m_B0_star (5.32475 + WIDTH)
#define m_B1_star (5.4154 + WIDTH)

#define m_01 (m_pi)
#define m_02 (m_K)
#define m_10 (m_pi)
#define m_12 (m_K)
#define m_20 (m_K)
#define m_21 (m_K)
#define m_22 (m_eta)
#define m_b (4.183)

#define m11 (5.27952)
#define m12 (5.32475)

#define m21 (5.36691)
#define m22 (5.4154)

#define m_Xb11P (9.909 - (m11 + m12))
#define m_Xb12P (10.224 - (m11 + m12))
#define m_Xb13P (10.529 - (m11 + m12))
#define m_Xb14P (10.852 - (m11 + m12))
#define m_Xb15P (11.114 - (m11 + m12))
#define m_Xb16P (11.395 - (m11 + m12))

#define a_l (0.06426 * 5.068)
#define gcoupling (0.898794378386677 / a_l * 10)
// #define g0 (gcoupling * 0.014926616931653945)
// #define g1 (gcoupling * 0.006467550544943349)
#define G0 0.035
#define G1 0.0183
// #define g1 (36. / 1000 * 1.4142135623730951)
#define delta0 0
#define delta1 (m21 + m22 - m11 - m12)
#define mu0 (m11 * m12 / (m11 + m12))
#define mu1 (m21 * m22 / (m21 + m22))
#define m_c 1.85
#define a_Cornell 1.95
#define g_qm (1)
#define partialwave 1

// wavefunction constants
#define N_MAX 64
#define N_TOWER 4
#define R_1 (0.02 * 5.068)       // GeV^-1
#define R_N_MAX (N_MAX * 5.068)  // GeV^-1
// #define C_T (0.19732 * 0.19732 / (0.5 * 4.18) * 5.068)
#define C_T (1 / m_b * 2)
#define PI 3.14159265358979323846
#define V0FIT 0.4
#define SIGMA_L 0.0199179973550142
// #define SIGMA (SIGMA_L / a_l / a_l)
#define SIGMA 0.198025
#define ALPHA -0.434

constexpr double complex delta[2] = {delta0, delta1};
constexpr double mu[2] = {mu0, mu1};
constexpr double Delta[2][2] = {{1, 0}, {0, 1}};
constexpr double mdiff[2] = {m_Xb12P - m_Xb11P, m_Xb13P - m_Xb11P};
constexpr double G[2] = {G0, G1};

#endif  // !DEBUG
