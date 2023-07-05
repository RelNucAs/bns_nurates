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

/* Print Gauss-Legendre quadrature
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

/* Print Gauss-Laguerre quadrature
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
 *
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