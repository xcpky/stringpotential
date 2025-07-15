#ifndef AUTOFREE_H
#define AUTOFREE_H
#include "lse.h"
#include "wavefunction.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <stdio.h>

static inline void matfree(matrix **mat) { matrix_free(*mat); }

static inline void permfree(gsl_permutation **perm) { gsl_permutation_free(*perm); }

static inline void lsefree(LSE **self) { lse_free(*self); }

static inline void wffree(WaveFunction **wf) { WFfree(*wf); }

static inline void auto_fclose(FILE **fp) { fclose(*fp); }

static inline void auto_ptrfree(double complex **ptr) { free(*ptr); }

static inline void autofree(void **ptr) { free(*ptr); }

static inline void close(FILE **ptr) { fclose(*ptr); }

#endif // AUTOFREE_H
