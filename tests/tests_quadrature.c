//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_quadrature.c
//  \brief perform tests with the quadrature routines

#include <stdio.h>
#include "../src/bns_nurates.h"
#include "../src/integration/integration.h"

void PrintGaussLegendreQuadrature(const int n) {
  MyQuadrature quad = {.dim =1, .type = kGauleg, .x1=-1., .x2=-1., .n = n};

  GaussLegendre(&quad);

  for (int i = 0; i < quad.n; i++) {
    printf("%f \t %f\n", quad.x[i], quad.w[i]);
  }
}

