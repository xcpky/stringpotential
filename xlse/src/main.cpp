#include "interface.h"
#include <span>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <matplot/matplot.h>
using namespace std;

vector<double> linspace(double start, double end,
                          size_t len) {
  auto step = (end - start) / (double)len;
  vector<double> res(len);
  for (size_t i = 0; i < len; i += 1) {
    res[i] = start + step * (double)i;
  }
  return res;
}

int main(int argc, char *argv[]) {
  size_t nsamples = atoi(argv[1]);
  auto E = linspace(-2.5, 0.5, nsamples);
  // puts("in cpp");
  // for (size_t i = 0; i < nsamples; i += 1) {
  //   printf("%f\t", E[i]);
  // }
  // puts("");
  auto ot = onshellT(E.data(), nsamples, 40);
  span<complex<double>> ot00(static_cast<complex<double>*>(ot.ose00), nsamples);
  vector<double> y(nsamples);
  for (size_t i = 0; i < nsamples; i += 1) {
    y[i] = abs(ot00[i]);
  }
  matplot::plot(E, y);
  // matplot::show();
  matplot::save("/home/zhy/code/stringpotential/xlse/onshellT.png");

  // free(E);
  ose_free(ot);
  return 0;
}
