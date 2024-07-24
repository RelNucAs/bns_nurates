//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_multi_func_2d.c
//  \brief test integration routines for multiple functions

#include <stdio.h>
#include <math.h>
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

MyQuadratureIntegrand TestIntegrand(double *x, void *p) {
  (void)p;
  MyQuadratureIntegrand result;

  result.integrand[0] = exp(-x[0]) * exp(-3.*x[1]);
  result.integrand[1] = exp(-2.*x[0]) * exp(-x[1]*x[1]);
  result.integrand[2] = exp(-x[0]*x[0]) * exp(-x[1]*x[1]*x[1]);
  result.integrand[3] = exp(-5.*x[0]) * exp(-x[1]);
  result.integrand[4] = exp(-x[0]) * exp(-x[1]);

  return result;
}

int main() {

  printf("==================================================\n");
  printf("Testing Multi-Function 2d integration routines ...\n");
  printf("==================================================\n");

  // function and parameters
  MyQuadratureIntegrand integrand_info = {.n = 5};

  MyFunctionMultiD test_function;
  test_function.function = &TestIntegrand;
  test_function.params = NULL;
  test_function.my_quadrature_integrand = integrand_info;

  MyQuadrature quad = quadrature_default;
  quad.dim = 2;
  quad.type = kGauleg;
  quad.x1 = 0.;
  quad.x2 = 1.;
  quad.nx = 45;
  quad.y1 = 0.;
  quad.y2 = 1.;
  quad.ny = 45;

  GaussLegendreMultiD(&quad);

  double s = 20.;
  MyQuadratureIntegrand gauleg_inf = GaussLegendreIntegrateFixedSplit2D(&quad, &test_function, s);

  double correct_answer[5] = {1./3., 0.4431134627263790068245418, 0.79138248703032128290432602198572815, 0.2, 1.};
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
