#include <stdio.h>
#include "src/constants.h"
#include "src/bns_nurates.h"
#include "src/integration/integration.h"
#include "tests/tests.h"

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
  return 0;
}
