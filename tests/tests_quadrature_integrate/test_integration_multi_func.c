//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_multi_func.c
//  \brief test integration routines for multiple functions

#include <stdio.h>
#include <math.h>
#include "../../src/bns_nurates.h"
#include "../../src/integration/integration.h"

MyQuadratureIntegrand TestIntegrand(double *x, void *p) {
  MyQuadratureIntegrand result;

  result.integrand[0] = exp(-x[0]);
  result.integrand[1] = exp(-2.*x[0]);
  result.integrand[2] = exp(-x[0]*x[0]);
  result.integrand[3] = exp(-5.*x[0]);

  return result;
}

int main() {

  printf("==================================================\n");
  printf("Testing Multi-Function 1d integration routines ...\n");
  printf("==================================================\n");

  double max_error = -42.;

  // function and parameters
  MyQuadratureIntegrand integrand_info = {.n = 4};

  MyFunctionMultiD test_function;
  test_function.function = &TestIntegrand;
  test_function.params = NULL;
  test_function.my_quadrature_integrand = integrand_info;

  MyQuadrature quad = quadrature_default;
  quad.dim = 1;
  quad.type = kGauleg;
  quad.x1 = 0.;
  quad.x2 = 1.;
  quad.nx = 45;

  GaussLegendreMultiD(&quad);

  double s = 1.;
  MyQuadratureIntegrand gauleg_inf = GaussLegendreIntegrate1D(&quad, &test_function, s);

  for (int i = 0; i < test_function.my_quadrature_integrand.n; i++) {
    printf("%d: \t %f\n", i, gauleg_inf.integrand[i]);
  }

  return 0;
}