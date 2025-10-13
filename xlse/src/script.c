#include "script.h"

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

#include "autofree.h"
#include "constants.h"
#include "lse.h"
#include "ome.h"
#include "wavefunction.h"

double complex* onshell(double* E, size_t len, size_t pNgauss, double Lambda,
                        double epsilon) {
    double complex(*res)[len] = malloc(sizeof(double complex) * len * 2);
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    for (size_t i = 0; i < len; i += 1) {
        lse_refresh(lse, E[i], (double[4]){0, 0, 0, 0}, G, PP);
        res[0][i] = lse->x0[0];
        res[1][i] = lse->x0[1];
    }
    return (double complex*)res;
}
double complex* onshellT(double* E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
    // puts("in C");
    // for (size_t i = 0; i < len; i += 1) {
    //   printf("%f,  ", E[i]);
    // }
    // puts("");
    double complex* res = malloc(sizeof(double complex) * len * 4);
    onshellElements onshellBuf = {
        .ose00 = res,
        .ose01 = res + len,
        .ose10 = res + 2 * len,
        .ose11 = res + 3 * len,
    };
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = &onshellBuf;
        args[i].rs = PP, args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], oT, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* onshellT_single(double* E, size_t len, double C[4],
                                size_t pNgauss, double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].E = E;
        args[i].id = i;
        args[i].size = len;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], oTsing, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* onshellG(double* E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len * 4);
    onshellElements onshellBuf = {
        .ose00 = res,
        .ose01 = res + len,
        .ose10 = res + 2 * len,
        .ose11 = res + 3 * len,
    };
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = &onshellBuf;
        args[i].rs = PP, args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], oG, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* onshellV(double* E, size_t len, double C[4], size_t pNgauss,
                         double Lambda, double epsilon) {
    // puts("in C");
    // for (size_t i = 0; i < len; i += 1) {
    //   printf("%f,  ", E[i]);
    // }
    // puts("");
    double complex* res = malloc(sizeof(double complex) * len * 4);
    onshellElements onshellBuf = {
        .ose00 = res,
        .ose01 = res + len,
        .ose10 = res + 2 * len,
        .ose11 = res + 3 * len,
    };
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = &onshellBuf;
        args[i].rs = PP;
        args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], oV, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* onshellTV(double* E, size_t len, double C[4], size_t pNgauss,
                          double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len * 8);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].size = len;
        args[i].rs = PP;
        args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], oTV, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

struct OME* ome_malloc() {
    struct OME* ome = malloc(sizeof(*ome));
    ome_build(ome);
    return ome;
}

double complex V(double E, double complex p, double complex pprime) {
    return ANA_00(E, p, pprime, m_pi);
}

double complex Vquad(double E, double complex p, double complex pprime) {
    static struct OME ome = {0};
    ome_build(&ome);
    return Vpiu(ome, E, p, pprime, m_B_star, m_B, m_B_star, m_B, m_pi, 1);
}

double complex* V3d(double E, size_t pNgauss, double Lambda, double epsilon) {
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse_refresh(lse, E, (double[4]){0, 0, 0, 0}, G, PP);
    lse_vmat(lse);
    double complex(*vmat)[2 * pNgauss + 2] =
        malloc(sizeof(double complex) * 2 * (pNgauss + 1) * 2 * (pNgauss + 1));
    for (size_t i = 0; i < pNgauss + 1; i += 1) {
        for (size_t j = 0; j < pNgauss + 1; j += 1) {
            vmat[i][j] = matrix_get(lse->VOME, i, j);
        }
    }
    return (double complex*)vmat;
    // return (double complex*) lse->VOME->data;
}

double complex* traceG(double* E, size_t len, double C[4], size_t pNgauss,
                       double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len * 2);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].E = E;
        args[i].id = i;
        args[i].size = len;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], trG, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* Det(double* E, size_t len, double C[4], uint64_t rs,
                    size_t pNgauss, double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = rs;
        args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], det, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* Det_single(double* E, size_t len, double C[4], uint64_t rs,
                           size_t pNgauss, double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        for (size_t j = 0; j < 4; j += 1) {
            args[i].C[j] = C[j];
        }
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].E = E;
        args[i].size = len;
        args[i].rs = rs;
        args[i].id = i;
        args[i].res = res;
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], detsing, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* Both(double* E, size_t len, double C[4], size_t pNgauss,
                     double Lambda, double epsilon) {
    double complex* res = malloc(sizeof(double complex) * len * 5);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].E = E;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], both, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

void* WF() {
    WaveFunction* wf = WFnew(partialwave, RLAMBDA, RNGAUSS);
    return wf;
}

double complex getwf(void* wf, double complex p, uint64_t n) {
    return psi_n_ftcomplex(wf, p, n);
}

double complex getwfx(void* wf, double r, uint64_t n) {
    return psi_n(wf, r, n, 0);
}

double complex getphi(double r, uint64_t n) {
    double complex phi = 0;
    double complex Ylm = gsl_sf_legendre_sphPlm(partialwave, 0, 1);
    Ylm = 1;
    // printf("%f%+f\n", creal(Ylm), cimag(Ylm));
    phi += N_nl((uint)n, partialwave) * pow(r, partialwave) * exp(-nu_n((uint)n) * r * r) * Ylm;
    return phi;
}

void WFend(void* wf) {
    WFfree(wf);
}

void Evec(double* E, double V0) {
    WaveFunction* wf [[gnu::cleanup(wffree)]] =
        WFnew(partialwave, RLAMBDA, RNGAUSS);
    // LSE *lse [[gnu::cleanup(lsefree)]] = lse_malloc(10, 4, 1e-6);
    for (size_t i = 0; i < N_TOWER; i += 1) {
        E[i] = wf->E_solution->data[i] + V0;
    }
}

void Ntower(size_t* len) { *len = N_TOWER; }

double complex* Poles(double* Er, size_t rlen, double* Ei, size_t ilen, double* g, size_t glen, double C[4], size_t pNgauss, double Lambda, double epsilon) {
    size_t len = rlen * ilen * glen;
    double complex* res = malloc(sizeof(double complex) * len * 4);
    thrd_t tid[NTHREADS];
    struct polestruct args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].Er = Er;
        args[i].rlen = rlen;
        args[i].Ei = Ei;
        args[i].ilen = ilen;
        args[i].g = g;
        args[i].glen = glen;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], poles, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

double complex* Polesm(double* pr, size_t rlen, double* pi, size_t ilen, double* g, size_t glen, double C[4], size_t pNgauss, double Lambda, double epsilon) {
    size_t len = rlen * ilen * glen;
    double complex* res = malloc(sizeof(double complex) * len);
    thrd_t tid[NTHREADS];
    struct polestruct args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;

    auto* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, (uint64_t)time(NULL));

    uint64_t* task = malloc(sizeof(uint64_t) * len);
    for (uint64_t i = 0; i < len; i += 1) {
        task[i] = i;
    }
    gsl_ran_shuffle(rng, task, len, sizeof(uint64_t));

    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
		args[i].task = task;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].Er = pr;
        args[i].rlen = rlen;
        args[i].Ei = pi;
        args[i].ilen = ilen;
        args[i].g = g;
        args[i].glen = glen;
        args[i].id = i;
        for (size_t cc = 0; cc < 4; cc += 1) {
            args[i].C[cc] = C[cc];
        }
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        // printf("id %zu\n", i);
        thrd_create(&tid[i], polesm, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }

    free(task);
    gsl_rng_free(rng);
    return res;
}

double* Fit(double* C, size_t len, const double g[NCHANNELS], size_t pNgauss, double Lambda, double epsilon) {
    void* result = malloc((NCHANNELS * NCHANNELS + 1) * len * sizeof(double));
    double (*Cin)[NCHANNELS * NCHANNELS] = (void*)C;
    double (*Cout)[NCHANNELS * NCHANNELS + 1] = result;
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse->g0 = g[0];
    lse->g1 = g[1];
    for (uint64_t i = 0; i < len; i += 1) {
        minimize(lse, Cin[i], Cout[i]);
    }
    return result;
}

// double* Fit(double* C, size_t len, size_t pNgauss, double Lambda, double epsilon) {
//     double* result = malloc((NCHANNELS * NCHANNELS + 1) * len * sizeof(double));
//     thrd_t tid[NTHREADS];
//     arg1d args[NTHREADS];
//     uint64_t ntasks = len / NTHREADS;
//     uint64_t residu = len % NTHREADS;
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         auto tr = args[thrd];
//         tr.E = C;
//         tr.epsilon = epsilon;
//         tr.Lambda = Lambda;
//         tr.pNgauss = pNgauss;
//         tr.size = len;
//         tr.id = thrd;
//         tr.rs = PP;
//         tr.res = result;
//         if (thrd < residu) {
//             tr.len = ntasks + 1;
//             tr.start = tr.len * thrd;
//         } else {
//             tr.len = ntasks;
//             tr.start = tr.len * thrd + residu;
//         }
//     }
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         thrd_create(&tid[thrd], thrdfit, &args[thrd]);
//     }
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         thrd_join(tid[thrd], NULL);
//     }
//     return result;
// }

int thrdfit(void* cargs) {
    auto arg = *(arg1d*)cargs;
    double (*Cin)[NCHANNELS * NCHANNELS] = (void*)arg.E;
    double (*Cout)[NCHANNELS * NCHANNELS + 1] = arg.res;
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(arg.pNgauss, arg.Lambda, arg.epsilon);
    for (uint64_t i = arg.start; i < arg.start + arg.len; i += 1) {
        minimize(lse, Cin[i], Cout[i]);
    }
    return 0;
}

double* Fitsing(double* C, size_t len, const double g[NCHANNELS], size_t pNgauss, double Lambda, double epsilon) {
    void* result = malloc((NCHANNELS * NCHANNELS + 1) * len * sizeof(double));
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    double (*Cin)[NCHANNELS * NCHANNELS] = (void*)C;
    double (*Cout)[NCHANNELS * NCHANNELS + 1] = result;
    lse->g0 = g[0];
    lse->g1 = g[1];
    for (uint64_t i = 0; i < len; i += 1) {
        minimizesing(lse, Cin[i], Cout[i]);
    }
    return result;
}

// double* Fitsing(double* C, size_t len, size_t pNgauss, double Lambda, double epsilon) {
//     double* result = malloc((NCHANNELS * NCHANNELS + 1) * len * sizeof(double));
//     thrd_t tid[NTHREADS];
//     arg1d args[NTHREADS];
//     uint64_t ntasks = len / NTHREADS;
//     uint64_t residu = len % NTHREADS;
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         auto tr = args[thrd];
//         tr.E = C;
//         tr.epsilon = epsilon;
//         tr.Lambda = Lambda;
//         tr.pNgauss = pNgauss;
//         tr.size = len;
//         tr.id = thrd;
//         tr.rs = PP;
//         tr.res = result;
//         if (thrd < residu) {
//             tr.len = ntasks + 1;
//             tr.start = tr.len * thrd;
//         } else {
//             tr.len = ntasks;
//             tr.start = tr.len * thrd + residu;
//         }
//     }
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         thrd_create(&tid[thrd], thrdfit, &args[thrd]);
//     }
//     for (uint64_t thrd = 0; thrd < NTHREADS; thrd += 1) {
//         thrd_join(tid[thrd], NULL);
//     }
//     return result;
// }

int thrdfitsing(void* cargs) {
    auto arg = *(arg1d*)cargs;
    double (*Cin)[NCHANNELS * NCHANNELS] = (void*)arg.E;
    double (*Cout)[NCHANNELS * NCHANNELS + 1] = arg.res;
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(arg.pNgauss, arg.Lambda, arg.epsilon);
    for (uint64_t i = arg.start; i < arg.start + arg.len; i += 1) {
        minimizesing(lse, Cin[i], Cout[i]);
    }
    return 0;
}

double* Cost(double* c, size_t len, [[maybe_unused]] double complex resonance, size_t pNgauss,
             double Lambda, double epsilon) {
    double* res = malloc(sizeof(double) * len);
    thrd_t tid[NTHREADS];
    arg1d args[NTHREADS];
    size_t ntasks = len / NTHREADS;
    size_t residue = len % NTHREADS;
    for (size_t i = 0; i < NTHREADS; i += 1) {
        args[i].pNgauss = pNgauss;
        args[i].Lambda = Lambda;
        args[i].epsilon = epsilon;
        args[i].res = res;
        args[i].rs = PP;
        args[i].E = c;
        args[i].id = i;
        if (i < residue) {
            args[i].len = ntasks + 1;
            args[i].start = (ntasks + 1) * i;
        } else {
            args[i].len = ntasks;
            args[i].start = ntasks * i + residue;
        }
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_create(&tid[i], cst, &args[i]);
    }
    for (size_t i = 0; i < NTHREADS; i += 1) {
        thrd_join(tid[i], NULL);
    }
    return res;
}

void ose_free(onshellElements ose) { free(ose.ose00); }

int oT(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    onshellElements* res = (onshellElements*)foo.res;
    double complex* ose00 = (double complex*)res->ose00;
    double complex* ose01 = (double complex*)res->ose01;
    double complex* ose10 = (double complex*)res->ose10;
    double complex* ose11 = (double complex*)res->ose11;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, G, foo.rs);
// printf("%f\n", foo.E[i]);
#ifdef TPO
        lse_X(lse);
        lse_XtX(lse);
        ose00[i] = lse->onshellT[0][0];
        ose01[i] = lse->onshellT[0][1];
        ose10[i] = lse->onshellT[1][0];
        ose11[i] = lse->onshellT[1][1];
#else
        // auto T = (double complex(*)[2 * ngauss + 2]) lse->TOME->data;
        ose00[i] = matrix_get(lse->TOME, ngauss - 3, ngauss - 3);
        ose01[i] = matrix_get(lse->TOME, ngauss - 3, 2 * ngauss - 3 + 1);
        ose10[i] = matrix_get(lse->TOME, 2 * ngauss - 3 + 1, ngauss - 3);
        ose11[i] = matrix_get(lse->TOME, 2 * ngauss - 3 + 1, 2 * ngauss - 3 + 1);
        // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
        // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I; idx += ngauss + 1; ose01[i] = lse->TOME->data[2 *
        // idx] + lse->TOME->data[2 * idx + 1] * I; idx = (2 * ngauss +
        // 1) * 2 * (ngauss + 1) + ngauss; ose10[i] = lse->TOME->data[2
        // * idx] + lse->TOME->data[2 * idx + 1] * I; idx += ngauss + 1;
        // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I;
#endif
    }
    return 0;
}

int oTsing(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    double complex(*res) = foo.res;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute_single(lse, foo.E[i], foo.C, G, foo.rs);
// printf("%f\n", foo.E[i]);
#ifdef TPO
        lse_X(lse);
        lse_XtX(lse);
        ose00[i] = lse->onshellT[0][0];
        ose01[i] = lse->onshellT[0][1];
        ose10[i] = lse->onshellT[1][0];
        ose11[i] = lse->onshellT[1][1];
#else
        [[maybe_unused]] auto T = (double complex(*)[2 * ngauss + 2]) lse->TOME->data;
        // res[i] = T[ngauss][ngauss];
        res[i] = matrix_get(lse->TOME, ngauss, ngauss);
        // res[i] = lse->onshellT[0][0];
        // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
        // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I; idx += ngauss + 1; ose01[i] = lse->TOME->data[2 *
        // idx] + lse->TOME->data[2 * idx + 1] * I; idx = (2 * ngauss +
        // 1) * 2 * (ngauss + 1) + ngauss; ose10[i] = lse->TOME->data[2
        // * idx] + lse->TOME->data[2 * idx + 1] * I; idx += ngauss + 1;
        // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I;
#endif
    }
    return 0;
}

int oG(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    size_t ngauss = foo.pNgauss;
    onshellElements* res = (onshellElements*)foo.res;
    double complex* ose00 = (double complex*)res->ose00;
    double complex* ose01 = (double complex*)res->ose01;
    double complex* ose10 = (double complex*)res->ose10;
    double complex* ose11 = (double complex*)res->ose11;
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, G, foo.rs);
        lse_gmat(lse);
        auto G = (double complex(*)[2 * ngauss + 2]) lse->G->data;
        ose00[i] = G[ngauss][ngauss];
        ose01[i] = G[ngauss][2 * ngauss + 1];
        ose10[i] = G[2 * ngauss + 1][ngauss];
        ose11[i] = G[2 * ngauss + 1][2 * ngauss + 1];
        // size_t idx = ngauss * 2 * (ngauss + 1) + ngauss;
        // ose00[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I; idx += ngauss + 1; ose01[i] = lse->TOME->data[2 *
        // idx] + lse->TOME->data[2 * idx + 1] * I; idx = (2 * ngauss +
        // 1) * 2 * (ngauss + 1) + ngauss; ose10[i] = lse->TOME->data[2
        // * idx] + lse->TOME->data[2 * idx + 1] * I; idx += ngauss + 1;
        // ose11[i] = lse->TOME->data[2 * idx] + lse->TOME->data[2 * idx
        // + 1] * I;
    }
    return 0;
}

int oV(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    int64_t ngauss = (int64_t)foo.pNgauss;
    onshellElements* res = (onshellElements*)foo.res;
    double complex* ose00 = (double complex*)res->ose00;
    double complex* ose01 = (double complex*)res->ose01;
    double complex* ose10 = (double complex*)res->ose10;
    double complex* ose11 = (double complex*)res->ose11;
    [[maybe_unused]] int64_t xoffset = 0;
    [[maybe_unused]] int64_t yoffset = 0;
    size_t x0 = (uint64_t)(ngauss + xoffset);
    size_t y0 = (uint64_t)(ngauss + yoffset);
    size_t x1 = (uint64_t)(ngauss + 2 * ngauss + 1 + xoffset);
    size_t y1 = (uint64_t)(ngauss + 2 * ngauss + 1 + yoffset);
    // printf("start: %lu, len: %lu\n", foo.start, foo.len);
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, G, foo.rs);
        lse_vmat(lse);
        ose00[i] = matrix_get(lse->VOME, x0, y0);
        ose01[i] = matrix_get(lse->VOME, x0, y1);
        ose10[i] = matrix_get(lse->VOME, x1, y0);
        ose11[i] = matrix_get(lse->VOME, x1, y1);
    }
    return 0;
}

int oTV(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    // auto buf = (double complex *)foo.res;
    auto buf = (double complex(*)[4][foo.size])foo.res;
    size_t pNgauss = foo.pNgauss;
    // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
    // for (size_t i = 0; i < 11; i += 1) {
    //   printf("%.2f ", creal(buf[i]));
    // }
    // puts("");
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, G, foo.rs);
        auto T = (double complex(*)[2 * pNgauss + 2]) lse->TOME->data;
        auto V = (double complex(*)[2 * pNgauss + 2]) lse->VOME->data;
        // printf("thread %lu, i %lu, det[i] %.1f\n", foo.id, i,
        // creal(det[i])); puts("");
        buf[0][0][i] = T[pNgauss][pNgauss];
        buf[0][1][i] = T[pNgauss][2 * pNgauss + 1];
        buf[0][2][i] = T[2 * pNgauss + 1][pNgauss];
        buf[0][3][i] = T[2 * pNgauss + 1][2 * pNgauss + 1];
        buf[1][0][i] = V[pNgauss][pNgauss];
        buf[1][1][i] = V[pNgauss][2 * pNgauss + 1];
        buf[1][2][i] = V[2 * pNgauss + 1][pNgauss];
        buf[1][3][i] = V[2 * pNgauss + 1][2 * pNgauss + 1];
    }
    return 0;
}

int det(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    double complex* res = (double complex*)foo.res;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        // res[i] = lse_detImVG(lse, foo.E[i], foo.C, PP);
        lse_compute(lse, foo.E[i], foo.C, G, foo.rs);
        res[i] = lse->det;
    }
    return 0;
}

int detsing(void* arg) {
    arg1d foo = *(arg1d*)arg;
    // printf("tr%zu %f\n",foo.id, foo.C[0]);
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    double complex* res = (double complex*)foo.res;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute_single(lse, foo.E[i], foo.C, G, foo.rs);
        res[i] = lse->det;
    }
    return 0;
}

int both(void* arg) {
    arg1d foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    auto buf = (double complex*)foo.res;
    size_t pNgauss = foo.pNgauss;
    // printf("thread %lu, foo.res %p\n", foo.id, foo.res);
    // for (size_t i = 0; i < 11; i += 1) {
    //   printf("%.2f ", creal(buf[i]));
    // }
    // puts("");
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_compute(lse, foo.E[i], foo.C, G, foo.rs);
        auto res = buf + 5 * i;
        auto T = (double complex(*)[2 * pNgauss + 2]) lse->TOME->data;
        // printf("thread %lu, i %lu, det[i] %.1f\n", foo.id, i,
        // creal(det[i])); puts("");
        res[0] = lse->det;
        res[1] = T[pNgauss][pNgauss];
        res[2] = T[pNgauss][2 * pNgauss + 1];
        res[3] = T[2 * pNgauss + 1][pNgauss];
        res[4] = T[2 * pNgauss + 1][2 * pNgauss + 1];
    }
    return 0;
}

int poles(void* arg) {
    auto foo = *(struct polestruct*)arg;
    // size_t len = foo.ilen * foo.rlen;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    // double complex* res = foo.res;
    double complex(*res)[foo.glen][foo.rlen][foo.ilen] = foo.res;
    for (size_t idx = foo.start; idx < foo.start + foo.len; idx += 1) {
        // [g, re, im]
        size_t g = idx / (foo.ilen * foo.rlen);
        size_t re = (idx / foo.ilen) % foo.rlen;
        size_t im = idx % foo.ilen;
        double complex E = foo.Er[re] + foo.Ei[im] * I;
        // res[idx] = pole(lse, E, foo.C, NN);
        // res[idx + len] = pole(lse, E, foo.C, PN);
        // res[idx + len * 2] = pole(lse, E, foo.C, NP);
        // res[idx + len * 3] = pole(lse, E, foo.C, (double[2]){foo.g[g], G[1]}, PP);
        res[PP][g][re][im] = pole(lse, E, foo.C, (double[2]){foo.g[g], G[1]}, PP);
    }
    return 0;
}

int polesm(void* arg) {
    auto foo = *(struct polestruct*)arg;
    // size_t len = foo.ilen * foo.rlen;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    // double complex* res = foo.res;
    double complex(*res)[foo.rlen][foo.ilen] = foo.res;
    // printf("thread[%zu] started, len %zu\n\n", foo.id, foo.len);
    auto mod = foo.len / 20;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
		auto idx = foo.task[i];
        // [g, re, im]
        // size_t g = idx / (foo.ilen * foo.rlen);
        // size_t re = (idx / foo.ilen) % foo.rlen;
        // size_t im = idx % foo.ilen;
        size_t im = idx / (foo.glen * foo.rlen);
        size_t re = (idx / foo.glen) % foo.rlen;
        size_t g = idx % foo.glen;
        double complex E = foo.Er[re] + foo.Ei[im] * I;
        // res[idx] = pole(lse, E, foo.C, NN);
        // res[idx + len] = pole(lse, E, foo.C, PN);
        // res[idx + len * 2] = pole(lse, E, foo.C, NP);
        // res[idx + len * 3] = pole(lse, E, foo.C, (double[2]){foo.g[g], G[1]}, PP);
        res[g][re][im] = polem(lse, E, foo.C, (double[2]){foo.g[g], G[1]}, PP);
        if ((i - foo.start) % mod == 0) {
            printf("id[%zu] at %zu%%\n", foo.id, (i - foo.start) / mod * 5);
        }
    }
    // printf("thread[%zu] finished\n", foo.id);
    return 0;
}

int cst(void* arg) {
    auto foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    double* res = foo.res;
    double C[4] = {0};
    C[0] = -1.010589943548671;
    C[2] = -1.220749787118462;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        C[3] = foo.E[i];
        res[i] = lse_cost(lse, C, PP);
    }
    return 0;
}

int trG(void* arg) {
    auto foo = *(arg1d*)arg;
    LSE* lse [[gnu::cleanup(lsefree)]] =
        lse_malloc(foo.pNgauss, foo.Lambda, foo.epsilon);
    auto res = (double complex(*)[foo.size])foo.res;
    auto pNgauss = foo.pNgauss;
    for (size_t i = foo.start; i < foo.start + foo.len; i += 1) {
        lse_refresh(lse, foo.E[i], foo.C, G, foo.rs);
        lse_gmat(lse);
        auto G = lse->G;
        for (size_t ch = 0; ch < 2; ch += 1) {
            res[ch][i] = 0;
            for (size_t p = 0; p < pNgauss + 1; p += 1) {
                // res[ch][i] += G[p + ch * (pNgauss + 1)][p + ch * (pNgauss + 1)];
                res[ch][i] += matrix_get(G, p + ch * (pNgauss + 1), p + ch * (pNgauss + 1));
            }
        }
    }
    return 0;
}

void Free(void* ptr) { free(ptr); }
double complex* getV(double complex E, size_t pNgauss, double Lambda, double epsilon) {
    // puts("getV get called");
    // printf("E: %f, pNgauss: %zu, Lambda: %f, epsilon: %f\n", E, pNgauss,
    // Lambda, epsilon);
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(pNgauss, Lambda, epsilon);
    lse_refresh(lse, E, (double[4]){0, 0, 0, 0}, G, PP);
    lse_vmat(lse);
    double complex(*V)[2 * pNgauss + 2] = malloc(sizeof(*V) * (2 * pNgauss + 2));
    for (size_t i = 0; i < 2 * pNgauss + 2; i += 1) {
        for (size_t j = 0; j < 2 * pNgauss + 2; j += 1) {
            V[i][j] = matrix_get(lse->VOME, i, j);
        }
    }
    return (double complex*)V;
}

double complex* getIntegrand() {
    double Lambda = 2;
    size_t Ngauss = 64;
    LSE* lse [[gnu::cleanup(lsefree)]] = lse_malloc(Ngauss, Lambda, 1e-6);
    double E = 0.1;
    lse_refresh(lse, E - m11 - m12 + 2 * m_B + m_pi, (double[4]){0, 0, 0, 0}, G, PP);
    double complex p1;
    double complex p2;
    p1 = lse->xi[12];
    p2 = lse->xi[1];
    double complex* res = malloc(sizeof(double complex) * (2 * DIMIM + DIMRE));
    for (size_t i = 0; i < 2 * DIMIM + DIMRE; i += 1) {
        res[i] = TOPTintegrand(E, lse->ome.xxpiup[i], p1, p2, m_B, m_B, m_pi);
    }
    return res;
}
