#include <stdio.h>
//#include "src/constants.h"
#include "src/bns_nurates.h"
#include "src/integration/integration.h"
#include "tests/tests.h"

struct FD_params {
  int k;
  double eta;
};
typedef struct FD_params FD_params;

double FD_legendre(double x, void *p) {
  struct FD_params * params = (struct FD_params *)p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(x,k) / (exp(x-eta)+1.);
}

double FD_laguerre(double x, void *p) {
  struct FD_params * params = (struct FD_params *)p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(x,k) / (exp(-eta)+exp(-x)); // FD_legendre divided by weight exp(-x)
}


int main() {

  printf("Hello, World!\n");

  //PrintGaussLegendreQuadrature(10, -1., 1.);
  printf("\n");
  PrintGaussLaguerreQuadrature(5, 0.);

  MyQuadrature quad;
  quad.n = 3;
  quad.dim = 1;
  quad.type = kGaulag;
  quad.x1 = -1.;
  quad.x2 = 1.;
  quad.alpha = 0.;

  SaveQuadrature("/var/home/maitraya/Documents/bns_nurates/src/quadratures/", &quad);

  /*=========================================================================*/

  // Minimum working example of numerical integration with GSL

  int k = 5;
  double eta = 10.;

  FD_params FD_p = {k, eta};

  gsl_function f;
  f.function = &FD_laguerre;
  f.params = &FD_p;

  const int n = 32; // number of quadrature point
  const double a = 0.;
  const double b = 1.;
  const double alpha = 0.;

  printf("\n");

  // Gauss-Laguerre integration (weight = (x-a)^alpha * exp(-b*(x-a)))
  double gLag = GSL_lag_quadr(n, a, b, alpha, &f);
  printf("Gauss-Laguerre integration from 0 to inf:   %.5e\n", gLag);

  f.function = &FD_legendre;

  // Gauss-Legendre integration from a to b (weight = 1)
  // First version
  double gLeg_1 = GSL_leg_quadr(n, a, b, &f);
  printf("Gauss-Legendre integration from a to b (1): %.5e\n", gLeg_1);

  // Second version
  double gLeg_2 = GSL_leg_integ(n, a, b, &f);
  printf("Gauss-Legendre integration from a to b (2): %.5e\n", gLeg_2);

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = 10.;
  double gLeg_inf = GSL_leg_split(n, &f, s);
  printf("Gauss-Legendre integration from 0 to inf:   %.5e\n", gLeg_inf);

  /*=========================================================================*/

  return 0;
}
