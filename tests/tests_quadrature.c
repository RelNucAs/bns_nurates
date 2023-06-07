//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_quadrature.c
//  \brief perform tests with the quadrature routines

#include <stdio.h>
#include <stdlib.h>
#include "../src/bns_nurates.h"
#include "../src/integration/integration.h"

void PrintGaussLegendreQuadrature(const int n, const double x1, const double x2) {
  MyQuadrature quad = {.dim =1, .type = kGauleg, .x1=x1, .x2=x2, .n = n};

  GaussLegendre(&quad);

  printf("Generating Gauss-Legendre quadratures ... \n");
  printf("x1 = %f \t x2 = %f \t n = %d\n", quad.x1, quad.x2, quad.n);
  printf("\n");
  for (int i = 0; i < quad.n; i++) {
    printf("%f \t %f\n", quad.x[i], quad.w[i]);
  }

  free(quad.x);
  free(quad.w);
}

void PrintGaussLaguerreQuadrature(const int n, const double alpha) {
  MyQuadrature quad = {.dim =1, .type = kGaulag, .n = n};

  GaussLaguerre(&quad, alpha);

  printf("Generating Gauss-Laguerre quadratures ... \n");
  printf("x1 = 0 \t x2 = inf \t n = %d\n", quad.n);
  printf("\n");
  for (int i = 0; i < quad.n; i++) {
    printf("%f \t %f\n", quad.x[i], quad.w[i]);
  }

  free(quad.x);
  free(quad.w);
}
