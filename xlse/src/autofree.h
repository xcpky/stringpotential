#ifndef AUTOFREE_H
#define AUTOFREE_H
#include "lse.h"
#include "wavefunction.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <stdio.h>

static inline void auto_matfree(matrix **mat) { matrix_free(*mat); }

static inline void auto_permfree(gsl_permutation **perm) {
  gsl_permutation_free(*perm);
}

static inline void auto_lsefree(LSE **self) { lse_free(*self); }

static inline void auto_wffree(WaveFunction **wf) { WFfree(*wf); }

static inline void auto_fclose(FILE **fp) { fclose(*fp); }

static inline void auto_ptrfree(double complex **ptr) { free(*ptr); }

#endif // AUTOFREE_H
