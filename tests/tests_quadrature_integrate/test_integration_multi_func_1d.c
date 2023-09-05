//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_multi_func_1d.c
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

  double correct_answer[4] = {1., 0.5, 0.886226925452758013649083741670, 0.2};
  double error = 42.;
  printf("Num \t Result \t Error\n");
  for (int i = 0; i < test_function.my_quadrature_integrand.n; i++) {
    printf("%d: \t %f \t %f\n", i, gauleg_inf.integrand[i],fabs(gauleg_inf.integrand[i]-correct_answer[i]));
    if(fabs(gauleg_inf.integrand[i]-correct_answer[i]) < fabs(error)) {
      error = fabs(gauleg_inf.integrand[i]-correct_answer[i]);
    }
  }

  if(error < 1e-7) {
    return 0;
  } else {
    return 1;
  }

}