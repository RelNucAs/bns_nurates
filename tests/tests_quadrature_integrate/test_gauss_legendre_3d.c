//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_gauss_legendre_3d.c
//  \brief Test 3d Gauss-Legendre routines with known values

#include <stdio.h>
#include <math.h>
#include "../../src/bns_nurates.h"
#include "../../src/integration/integration.h"

int main() {

  double max_error = -42.;
  double points[9] = {-(1. / 3.) * sqrt(3.), (1. / 3.) * sqrt(3.), -(1. / 35.) * sqrt(525. + 70. * sqrt(30.)), -(1. / 35.) * sqrt(525. - 70. * sqrt(30.)),
                      (1. / 35.) * sqrt(525. - 70. * sqrt(30.)), (1. / 35.) * sqrt(525. + 70. * sqrt(30.)), -(1. / 5.) * sqrt(15.), 0, (1. / 5.) * sqrt(15.)};
  double weights[9] = {1., 1., (1. / 36.) * (18. - sqrt(30.)), (1. / 36.) * (18. + sqrt(30.)), (1. / 36.) * (18. + sqrt(30.)), (1. / 36.) * (18. - sqrt(30.)), 5. / 9., 8. / 9.,
                       5. / 9.};

  MyQuadrature quad = quadrature_default;
  quad.dim = 3;
  quad.x1 = -1.;
  quad.x2 = 1.;
  quad.nx = 2;
  quad.y1 = -1.;
  quad.y2 = 1.;
  quad.ny = 4;
  quad.y1 = -1.;
  quad.y2 = 1.;
  quad.nz = 3;
  quad.z1 = -1.;
  quad.z2 = 1.;
  // Generate 3d quadratures
  printf("Generating Gauss-Legendre quadratures in 3d ... \n");

  GaussLegendreMultiD(&quad);

  printf("x1 = %f \t x2 = %f \t nx = %d \t y1 = %f \t y2 = %f \t ny = %d \t z1 = %f \t z2 = %f \t nz = %d\n",
         quad.x1, quad.x2, quad.nx, quad.y1, quad.y2, quad.ny, quad.z1, quad.z2, quad.nz);
  printf("\n");

  printf("points \t answer \t weights \t answer\n");
  for (int i = 0; i < quad.nx + quad.ny + quad.nz; i++) {
    if (i == quad.nx || i == quad.nx + quad.ny) {
      printf("\n");
    }

    printf("%f \t %f \t %f \t %f\n", quad.points[i], points[i], quad.w[i], weights[i]);

    if (fabs(quad.points[i] - points[i]) > max_error) {
      max_error = fabs(quad.points[i] - points[i]);
    }
    if (fabs(quad.w[i] - weights[i]) > max_error) {
      max_error = quad.w[i] - weights[i];
    }
  }

  if (max_error < 1e-12) {
    return 0;
  } else {
    return 1;
  }

}