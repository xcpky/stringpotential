#ifndef OME_H
#define OME_H
#include "constants.h"
#include <math.h>
#include <stddef.h>
#define FSQUARE(F) (F) * (F)
#define EPSILON (1e-9)
#define DIMIM (16)
#define DIMRE (24)
#define ZI (0.3)
#define FACPI (-3 * g_b * g_b / f_pi / f_pi / 24)
#include <complex.h>

static inline double complex csquare(double complex x) { return x * x; }
static inline double fsquare(double x) { return x * x; }
// constexpr double mupisquare00 = FSQUARE(m_pi) - FSQUARE(m_B_star - m_B);
// constexpr double muetasquare00 = FSQUARE(m_eta) - FSQUARE(m_B_star - m_B);
// constexpr double muKsquare01 = FSQUARE(m_K) - FSQUARE(m_B_star_s - m_B);
// constexpr double muetasquare11 = FSQUARE(m_eta) - FSQUARE(m_B_star_s - m_B_s);
constexpr double mupisquare00 = 1;
constexpr double muetasquare00 = 1;
constexpr double muKsquare01 = 1;
constexpr double muetasquare11 = 1;

struct OME {
      double complex xxpiup[2 * DIMIM + DIMRE];
      double complex wwpiup[2 * DIMIM + DIMRE];
      double complex xxpidn[2 * DIMIM + DIMRE];
      double complex wwpidn[2 * DIMIM + DIMRE];
      double xxz[DIMRE];
      double wwz[DIMRE];
      double xx0ii[DIMRE];
      double ww0ii[DIMRE];
};

void ome_build(struct OME *self);
void ome_free(struct OME *self);

static inline double complex Epi(double complex z, double complex p1, double complex p2, double m0)
{
      return csqrt(p1 * p1 + p2 * p2 - 2 * p1 * p2 * z + m0 * m0);
}

static inline double complex Dij(double complex E, double complex z, double complex p1, double complex p2, double complex mi,
				 double complex mj, double m0)
{
      return E - (mi + p1 * p1 / (2 * mi)) - (mj + p2 * p2 / (2 * mj)) - Epi(z, p1, p2, m0) + I * EPSILON;
}

static inline double complex z0(double complex E, double complex m, double complex p1, double complex p2, double m0)
{
      return (csquare(E - 2 * m - (p1 * p1 + p2 * p2) / (2 * m)) - m0 * m0 - p1 * p1 - p2 * p2) / (-2 * p2 * p1);
}

static inline double complex z0E(double complex p1, double complex p2, double m0)
{
      return (p1 * p1 + p2 * p2 + m0 * m0) / (2 * p2 * p1);
}

double complex Vpiu(struct OME ome, double complex E, double complex p1, double complex p2, double m1, double gam1, double m2,
		    double gam2, double m3, double gam3, double m4, double gam4, double m0, double fac);

static inline double complex quadreal(struct OME ome, double complex E, double complex p1, double complex p2, double m1,
				      double gam1, double m2, double gam2, double m3, double gam3, double m4, double gam4,
				      double m0, double fac)
{
      double complex res = 0;
      for (size_t i = 0; i < DIMRE; i += 1) {
	    auto z = ome.xxz[i];
	    auto w = ome.wwz[i];
	    auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2, m0);
	    auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 / 2, m0);
	    auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1, p2, m0)) * (p1 * p1 + p2 * p2 - 2 * p1 * p2 * z);
	    res += Dint * w;
      }
      return res;
}
static inline double complex quadup(struct OME ome, double complex E, double complex p1, double complex p2, double m1,
				    double gam1, double m2, double gam2, double m3, double gam3, double m4, double gam4,
				    double m0, double fac)
{
      double complex res = 0;
      for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
	    auto z = ome.xxpiup[i];
	    auto w = ome.wwpiup[i];
	    auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2, m0);
	    auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 / 2, m0);
	    auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1, p2, m0)) * (p1 * p1 + p2 * p2 - 2 * p1 * p2 * z);
	    res += Dint * w;
      }
      return res;
}
static inline double complex quaddn(struct OME ome, double complex E, double complex p1, double complex p2, double m1,
				    double gam1, double m2, double gam2, double m3, double gam3, double m4, double gam4,
				    double m0, double fac)
{
      double complex res = 0;
      for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
	    auto z = ome.xxpidn[i];
	    auto w = ome.wwpidn[i];
	    auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2, m0);
	    auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 / 2, m0);
	    auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1, p2, m0)) * (p1 * p1 + p2 * p2 - 2 * p1 * p2 * z);
	    res += Dint * w;
      }
      return res;
}

static inline double complex quadii(struct OME ome, double complex E, double complex p1, double complex p2, double m1,
				    double gam1, double m2, double gam2, double m3, double gam3, double m4, double gam4,
				    double m0, double fac, double complex _z, double complex offset)
{
      double complex res = 0;
      for (size_t i = 0; i < DIMRE; i += 1) {
	    auto z = ome.xx0ii[i] * _z + offset;
	    auto w = ome.ww0ii[i] * _z;
	    auto D1 = Dij(E, z, p1, p2, m1 - I * gam1 / 2, m3 - I * gam3 / 2, m0);
	    auto D2 = Dij(E, z, p1, p2, m2 - I * gam2 / 2, m4 - I * gam4 / 2, m0);
	    auto Dint = FACPI * fac * (1 / D1 + 1 / D2) / (2 * Epi(z, p1, p2, m0)) * (p1 * p1 + p2 * p2 - 2 * p1 * p2 * z);
	    res += Dint * w;
      }
      return res;
}

static inline double complex OME_00(struct OME ome, double complex E, double complex p, double complex pprime)
{
      return Vpiu(ome, E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_pi, 3 * mupisquare00) +
	     Vpiu(ome, E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star, gamma_B_star, m_B, 0, m_eta,
		  1. / 3 * muetasquare00);
}

static inline double complex OME_01(struct OME ome, double complex E, double complex p, double complex pprime)
{
      return Vpiu(ome, E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star, gamma_B_star, m_B, 0, m_K,
		  pow(2, 3. / 2) * muKsquare01);
}

static inline double complex OME_10(struct OME ome, double complex E, double complex p, double complex pprime)
{
      return Vpiu(ome, E, p, pprime, m_B_star, gamma_B_star, m_B, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_K,
		  pow(2, 3. / 2) * muKsquare01);
}

static inline double complex OME_11(struct OME ome, double complex E, double complex p, double complex pprime)
{
      return Vpiu(ome, E, p, pprime, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_B_star_s, gamma_B_star_s, m_B_s, 0, m_eta,
		  4. / 3 * muetasquare11);
}
#endif // OME_H
