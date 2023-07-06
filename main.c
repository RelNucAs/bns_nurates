#include <stdio.h>
//#include "src/constants.h"
#include "src/bns_nurates.h"
#include "src/integration/integration.h"
#include "tests/tests.h"
#include "src/opacities/kernels/kernels.h"

int main() {

  printf("Hello, World!\n");

  //TestQuadratureInputOutput("/var/home/maitraya/Documents/bns_nurates/src/quadratures/");

  TestGaussLaguerreQuadrature();
  //PrintGaussLegendreQuadrature(10, -1., 1.);
  //printf("\n");
  //PrintGaussLaguerreQuadrature(5, 0.);
  //TestQuadratureWithGSL();
  //TestPairT();
  //printf("\n");
  //TestPairF();
  //printf("\n");
  //TestPairG();
  //printf("\n");
  //TestPairPsi();
  //printf("\n");
  //TestPairPhi();
  //printf("\n");
  //TestPairKernels();
  //printf("\n");
  //TestPairOpacities();
  //printf("\n");
  //TestBremKernelS("/var/home/maitraya/Documents/");

  //BremKernelS(-0.5, 0., 0.31);
  //TestBremKernelG("/var/home/maitraya/Documents/");
  //MyQuadrature quad;
  //quad.n = 3;
  //quad.dim = 1;
  //quad.type = kGaulag;
  //quad.x1 = -1.;
  //quad.x2 = 1.;
  //quad.alpha = 0.;

  //SaveQuadrature("/var/home/maitraya/Documents/bns_nurates/src/quadratures/", &quad);



  return 0;
}
