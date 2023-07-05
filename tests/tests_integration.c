//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_integration.c
//  \brief test integration routines

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "tests.h"
#include "../src/bns_nurates.h"
#include "../src/integration/integration.h"

// parameters for test functions
struct FermiDiracParams {
  int k;
  double eta;
};
typedef struct FermiDiracParams FermiDiracParams;

/* Fermi-Dirac function
 *
 * Inputs:
 *      points:  the integration variable
 *      p:  function parameters
 */
double FermiDiracLegendre(double *x, void *p) {
  struct FermiDiracParams *params = (struct FermiDiracParams *) p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(*x, k) / (exp(*x - eta) + 1.);
}

double FermiDiracLegendreGSL(double x, void *p) {
  return FermiDiracLegendre(&x, p);
}
/* Fermi-Dirac function divided by exp(-points) designed for use
 * with Gauss-Laguerre routines
 *
 * Inputs:
 *      points:  the integration variable
 *      p:  function parameters
 */
double FermiDiracLaguerre(double *x, void *p) {
  struct FermiDiracParams *params = (struct FermiDiracParams *) p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(*x, k) / (exp(-eta) + exp(-*x));
}

double FermiDiracLaguerreGSL(double x, void *p) {
  return FermiDiracLaguerre(&x, p);
}

double SimpleLaguerreFunc(double x, void *p) {
  return 1./(x*x*exp(x));
}

void TestQuadratureWithGSL() {

  // function and parameters
  FermiDiracParams fd_p = {.k = 5, .eta = 10};
  gsl_function f;
  f.function = &FermiDiracLaguerreGSL;
  f.params = &fd_p;

  MyFunction fermi;
  fermi.function = &FermiDiracLaguerre;
  fermi.params = &fd_p;

  // test Gauss-Laguerre implementation
  const int n = 35;
  const double alpha = 0.;
  MyQuadrature quad = {.nx = n, .dim = 1, .type = kGaulag, .x1 = 0., .x2 = 1., .alpha = alpha};
  GaussLaguerre(&quad);

  // Gauss-Laguerre integration (weight = (points-a)^alpha * exp(-b*(points-a)))
  const double a = 0.;
  const double b = 1.;
  double gaulag_result_gsl = GSLLagQuadrature(n, a, b, alpha, &f);
  double gaulag_result = GaussLaguerreIntegrateZeroInf(&quad, &fermi);
  printf("Gauss-Laguerre integration from 0 to inf (GSL):   %.5e\n", gaulag_result_gsl);
  printf("Gauss-Laguerre integration from 0 to inf (our):   %.5e\n", gaulag_result);

  // Gauss-Legendre integration tests
  f.function = &FermiDiracLegendreGSL;
  fermi.function = &FermiDiracLegendre;


  // Gauss-Legendre integration from 0 to 1 (weight = 1)
  // first version
  double gauleg_result_gsl_1 = GSLLegQuadrature(n, a, b, &f);
  printf("\n");
  printf("Gauss-Legendre integration from 0 to 1 (GSL type 1): %.5e\n", gauleg_result_gsl_1);
  // second version
  double gauleg_result_gsl_2 = GSLLegIntegrate(n, a, b, &f);
  printf("Gauss-Legendre integration from 0 to 1 (GSL type 2): %.5e\n", gauleg_result_gsl_2);
  // our version
  quad.type = kGauleg;
  GaussLegendre(&quad);
  double gauleg_result = 0.;
  for (int i = 0; i < quad.nx; i++) {
    gauleg_result += FermiDiracLegendre(&quad.points[i], &fd_p) * quad.w[i];
  }
  printf("Gauss-Legendre integration from 0 to 1 (   ours   ): %.5e\n", gauleg_result);

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = 10.;
  double gauleg_inf_gsl = GslLegInfSplit(n, &f, s);
  double gauleg_inf = GaussLegendreIntegrateZeroInf(&quad, &fermi, s);
  printf("\n");
  printf("Gauss-Legendre integration of Fermi function from 0 to inf (GSL):   %.5e\n", gauleg_inf_gsl);
  printf("Gauss-Legendre integration of Fermi function from 0 to inf (our):   %.5e\n", gauleg_inf);

}