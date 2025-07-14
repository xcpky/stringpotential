#include "ome.h"
#include <complex.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>

void ome_build(struct OME *self) {
  gsl_integration_glfixed_table *tableim =
      gsl_integration_glfixed_table_alloc(DIMIM);
  gsl_integration_glfixed_table *tablere =
      gsl_integration_glfixed_table_alloc(DIMRE);
  for (size_t i = 0; i < DIMIM; i += 1) {
    double im;
    double wi;
    gsl_integration_glfixed_point(0, ZI, i, &im, &wi, tableim);
    self->xxpiup[i] = -1 + im * I;
    self->wwpiup[i] = I * wi;

    self->xxpiup[i + DIMIM + DIMRE] = 1 + (ZI - im) * I;
    self->wwpiup[i + DIMIM + DIMRE] = -I * wi;

    self->xxpidn[i] = -1 - im * I;
    self->wwpidn[i] = -I * wi;

    self->xxpidn[i + DIMIM + DIMRE] = 1 + (im - ZI) * I;
    self->wwpidn[i + DIMIM + DIMRE] = I * wi;
  }
  for (size_t i = 0; i < DIMRE; i += 1) {
    double xi, wi;
    gsl_integration_glfixed_point(-1, 1, i, &xi, &wi, tablere);
    self->xxpiup[i + DIMIM] = xi + ZI * I;
    self->wwpiup[i + DIMIM] = wi;
    self->xxpidn[i + DIMIM] = xi - ZI * I;
    self->wwpidn[i + DIMIM] = wi;
  }

  for (size_t i = 0; i < DIMRE; i += 1) {
    gsl_integration_glfixed_point(0, 1, i, &self->xx0ii[i], &self->ww0ii[i],
                                  tablere);
  }
  for (size_t i = 0; i < DIMRE; i += 1) {
    gsl_integration_glfixed_point(-1, 1, i, &self->xxz[i], &self->wwz[i],
                                  tablere);
  }
  gsl_integration_glfixed_table_free(tableim);
  gsl_integration_glfixed_table_free(tablere);
}

void ome_free(struct OME *self) { free(self); }

double complex Vpiu(struct OME ome, _Complex double E, _Complex double p1,
                    _Complex double p2, double m1, double gam1, double m2,
                    double gam2, double m3, double gam3, double m4, double gam4,
                    double m0, double fac) {
  double frame = 0.2;
  auto _z0 = z0(E, m1 - I * gam1 / 2, p1, p2, m0);
  auto _z0E = z0E(p1, p2, m0);
  // auto foo = m1 - I*gam1;
  // printf("%f%+f\n", creal(foo), cimag(foo));
  // printf("%f%+f\n", creal(_z0), cimag(_z0));
  if (fabs(creal(_z0)) > 1 || fabs(cimag(_z0)) > frame) {
    if (fabs(creal(_z0E)) > 1 || fabs(cimag(_z0)) > frame) {
      return quadreal(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4,
                      m0, fac);
    } else if (cimag(_z0E) > 0) {
      return quaddn(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac);
    } else {
      return quadup(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac);
    }
  } else if (cimag(_z0) > 0) {
    if (fabs(creal(_z0E)) > 1 || fabs(cimag(_z0E)) > frame || cimag(_z0E) > 0) {
      return quaddn(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac);
    } else {
      double complex _z01, _z02;
      if (creal(_z0) < creal(_z0E)) {
        _z01 = creal(_z0) - 0.1 * I;
        _z02 = creal(_z0E) + 0.1 * I;
      } else {
        _z01 = creal(_z0E) + 0.1 * I;
        _z02 = creal(_z0) - 0.1 * I;
      }
      return quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, _z01 + 1, -1) +
             quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, _z02 - _z01, _z01) +
             quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, 1 - _z02, _z02);
    }
  } else {
    if (fabs(creal(_z0E)) > 1 || fabs(cimag(_z0E)) > frame) {
      return quadup(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac);
    } else if (cimag(_z0E) > 0) {
      double complex _z01, _z02;
      if (creal(_z0) < creal(_z0E)) {
        _z01 = creal(_z0) + 0.1 * I;
        _z02 = creal(_z0E) - 0.1 * I;
      } else {
        _z01 = creal(_z0E) - 0.1 * I;
        _z02 = creal(_z0) + 0.1 * I;
      }
      return quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, _z01 + 1, -1) +
             quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, _z02 - _z01, _z01) +
             quadii(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac, 1 - _z02, _z02);
    } else {
      return quadup(ome, E, p1, p2, m1, gam1, m2, gam2, m3, gam3, m4, gam4, m0,
                    fac);
    }
  }
}
