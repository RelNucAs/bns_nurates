//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_quadrature.c
//  \brief perform tests with quadrature generation and IO

#include <stdio.h>
#include <stdlib.h>
#include "../src/bns_nurates.h"
#include "../src/integration/integration.h"

/* Print the 1d Gauss-Legendre quadrature
 *
 * Inputs:
 *      nx:  the number of quadrature points
 *      x1: lower limit of points
 *      x2: upper limit of points
 */
void PrintGaussLegendreQuadrature(const int n, const double x1, const double x2) {
  MyQuadrature quad = {.dim =1, .type = kGauleg, .x1=x1, .x2=x2, .nx = n};

  GaussLegendre(&quad);

  printf("Generating Gauss-Legendre quadratures ... \n");
  printf("x1 = %f \t x2 = %f \t nx = %d\nx", quad.x1, quad.x2, quad.nx);
  printf("\n");
  for (int i = 0; i < quad.nx; i++) {
    printf("%f \t %f\n", quad.points[i], quad.w[i]);
  }

  free(quad.points);
  free(quad.w);
}

/* Generate Gauss-Legendre quadrature table from Wolfram:
 * https://mathworld.wolfram.com/Legendre-GaussQuadrature.html
 */
void TestGaussLegendreQuadrature() {
  PrintGaussLegendreQuadrature(2, -1., 1.);
  printf("\n");
  PrintGaussLegendreQuadrature(3, -1., 1.);
  printf("\n");
  PrintGaussLegendreQuadrature(4, -1., 1.);
  printf("\n");
  PrintGaussLegendreQuadrature(5, -1., 1.);
}

/* Generate Gauss-Legendre quadrature table from Wolfram:
 * https://mathworld.wolfram.com/Legendre-GaussQuadrature.html
 */
void TestGaussLegendreQuadratureMultiD() {
  MyQuadrature quad = quadrature_default;

  // Generate 1d quadratures
  printf("Generating Gauss-Legendre quadratures in 1d ... \n");
  quad.x1 = -1.;
  quad.x2 = 1.;
  quad.nx = 3;

  GaussLegendreMultiD(&quad);

  printf("x1 = %f \t x2 = %f \t nx = %d \t y1 = %f \t y2 = %f \t ny = %d \t z1 = %f \t z2 = %f \t nz = %d\n",
         quad.x1, quad.x2, quad.nx, quad.y1, quad.y2, quad.ny, quad.z1, quad.z2, quad.nz);
  printf("\n");

  for (int i = 0; i < quad.nx + quad.ny + quad.nz; i++) {
    if (i == quad.nx || i == quad.nx + quad.ny) {
      printf("\n");
    }

    printf("%f \t %f\n", quad.points[i], quad.w[i]);

  }

  printf("\n");

  // Generate 2d quadratures
  printf("Generating Gauss-Legendre quadratures in 2d ... \n");
  quad.dim = 2;
  quad.x1 = -1.;
  quad.x2 = 1.;
  quad.nx = 3;
  quad.y1 = -1.;
  quad.y2 = 1.;
  quad.ny = 4;

  GaussLegendreMultiD(&quad);

  printf("x1 = %f \t x2 = %f \t nx = %d \t y1 = %f \t y2 = %f \t ny = %d \t z1 = %f \t z2 = %f \t nz = %d\n",
         quad.x1, quad.x2, quad.nx, quad.y1, quad.y2, quad.ny, quad.z1, quad.z2, quad.nz);
  printf("\n");

  for (int i = 0; i < quad.nx + quad.ny + quad.nz; i++) {
    if (i == quad.nx || i == quad.nx + quad.ny) {
      printf("\n");
    }

    printf("%f \t %f\n", quad.points[i], quad.w[i]);

  }

  printf("\n");

  // Generate 2d quadratures
  printf("Generating Gauss-Legendre quadratures in 3d ... \n");
  quad.dim = 3;
  quad.x1 = -1.;
  quad.x2 = 1.;
  quad.nx = 3;
  quad.y1 = -1.;
  quad.y2 = 1.;
  quad.ny = 4;
  quad.z1 = -1;
  quad.z2 = 1;
  quad.nz = 5;

  GaussLegendreMultiD(&quad);

  printf("x1 = %f \t x2 = %f \t nx = %d \t y1 = %f \t y2 = %f \t ny = %d \t z1 = %f \t z2 = %f \t nz = %d\n",
         quad.x1, quad.x2, quad.nx, quad.y1, quad.y2, quad.ny, quad.z1, quad.z2, quad.nz);
  printf("\n");

  for (int i = 0; i < quad.nx + quad.ny + quad.nz; i++) {
    if (i == quad.nx || i == quad.nx + quad.ny) {
      printf("\n");
    }

    printf("%f \t %f\n", quad.points[i], quad.w[i]);

  }

}
/* Print Gauss-Laguerre quadrature in 1d
 *
 * Inputs:
 *      nx:      the number of quadrature points
 *      alpha:  parameter for the weighting function W(points) = points^alpha e^-points
 */
void PrintGaussLaguerreQuadrature(const int n, const double alpha) {
  MyQuadrature quad = {.dim =1, .type = kGaulag, .nx = n, .alpha = alpha};

  GaussLaguerre(&quad);

  printf("Generating Gauss-Laguerre quadratures ... \n");
  printf("x1 = 0 \t x2 = inf \t nx = %d\nx", quad.nx);
  printf("\n");
  for (int i = 0; i < quad.nx; i++) {
    printf("%f \t %f\n", quad.points[i], quad.w[i]);
  }

  free(quad.points);
  free(quad.w);
}

/* Generate Gauss-Laguerre table from Wolfram:
 * https://mathworld.wolfram.com/Laguerre-GaussQuadrature.html
 */
void TestGaussLaguerreQuadrature() {
  PrintGaussLaguerreQuadrature(2, 0.);
  printf("\n");
  PrintGaussLaguerreQuadrature(3, 0.);
  printf("\n");
  PrintGaussLaguerreQuadrature(4, 0.);
  printf("\n");
  PrintGaussLaguerreQuadrature(5, 0.);
}

/* Test the input/output routines for the quadratures
 * Note: The save only works for 1d quadratures.
 * Input:
 *    filedir:  quadrature save directory path
 */
void TestQuadratureInputOutput(char *filedir) {
  // cheking quadrature creating and saving to disk (Gauss-Legendre)
  MyQuadrature quad = {.dim = 1, .type=kGauleg, .nx = 9, .x1 = -1., .x2 = 1.};
  printf("Saving quadrature information for Gauss-Legendre...\n");
  SaveQuadrature(filedir, &quad);
  printf("Loading saved quadrature ...\n");
  LoadQuadrature(filedir, &quad);
  printf("Print quadrature metadata/data ...\n");
  printf("nx = %d \t x1 = %0.16e \t x2 = %0.16e\nx", quad.nx, quad.x1, quad.x2);
  for (int i = 0; i < quad.nx; i++) {
    printf("%0.16e \t %0.16e\n", quad.points[i], quad.w[i]);
  }

  // checking loading from disk (Gauss-Legendre)
  quad.x1 = 0;
  quad.nx = 6;
  LoadQuadrature(filedir, &quad);
  printf("\n");
  printf("Print quadrature metadata/data for another quadrature ...\n");
  printf("nx = %d \t x1 = %0.16e \t x2 = %0.16e\nx", quad.nx, quad.x1, quad.x2);
  for (int i = 0; i < quad.nx; i++) {
    printf("%0.16e \t %0.16e\n", quad.points[i], quad.w[i]);
  }

  // cheking quadrature creating and saving to disk (Gauss-Laguerre)
  quad.dim = 1;
  quad.type = kGaulag;
  quad.nx = 9;
  quad.alpha = 0.;

  printf("Saving quadrature information for Gauss-Laguerre...\n");
  SaveQuadrature(filedir, &quad);
  printf("Loading saved quadrature ...\n");
  LoadQuadrature(filedir, &quad);
  printf("Print quadrature metadata/data ...\n");
  printf("nx = %d \t x1 = %0.16e \t x2 = %0.16e\nx", quad.nx, quad.x1, quad.x2);
  for (int i = 0; i < quad.nx; i++) {
    printf("%0.16e \t %0.16e\n", quad.points[i], quad.w[i]);
  }

  // checking loading from disk (Gauss-Laguerre)
  quad.nx = 6;
  LoadQuadrature(filedir, &quad);
  printf("\n");
  printf("Print quadrature metadata/data for another quadrature ...\n");
  printf("nx = %d \t x1 = %0.16e \t x2 = %0.16e\nx", quad.nx, quad.x1, quad.x2);
  for (int i = 0; i < quad.nx; i++) {
    printf("%0.16e \t %0.16e\n", quad.points[i], quad.w[i]);
  }

}