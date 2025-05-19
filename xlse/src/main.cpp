#include "interface.h"
#include <span>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <matplot/matplot.h>
using namespace std;

double *linspace(double start, double end,
                          size_t len) {
  auto step = (end - start) / (double)len;
  auto res = (double *)malloc(sizeof(double) * len);
  for (size_t i = 0; i < len; i += 1) {
    res[i] = start + step * (double)i;
  }
  return res;
}

int main(int argc, char *argv[]) {
  size_t nsamples = 10;
  vector<double> E(nsamples);
  auto tmp = linspace(-2.5, 0.5, nsamples);
  for (size_t i = 0; i < nsamples; i += 1) {
    E[i] = tmp[i];
  }
  auto ot = onshellT(E.data(), nsamples, 40);
  span<complex<double>> ot00(static_cast<complex<double>*>(ot.ose00), nsamples);
  vector<double> y(nsamples);
  for (size_t i = 0; i < nsamples; i += 1) {
    y[i] = abs(ot00[i]);
  }
  matplot::plot(E, y);
  matplot::show();

  // free(E);
  ose_free(ot);
  return 0;
}
