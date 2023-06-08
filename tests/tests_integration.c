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
 *      x:  the integration variable
 *      p:  function parameters
 */
double FermiDiracLegendre(double x, void *p) {
  struct FermiDiracParams *params = (struct FermiDiracParams *) p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(x, k) / (exp(x - eta) + 1.);
}

/* Fermi-Dirac function divided by exp(-x) designed for use
 * with Gauss-Laguerre routines
 *
 * Inputs:
 *      x:  the integration variable
 *      p:  function parameters
 */
double FermiDiracLaguerre(double x, void *p) {
  struct FermiDiracParams *params = (struct FermiDiracParams *) p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(x, k) / (exp(-eta) + exp(-x));
}

void TestQuadratureWithGSL() {

  // function and parameters
  FermiDiracParams fd_p = {.k = 5, .eta = 10};
  gsl_function f;
  f.function = &FermiDiracLaguerre;
  f.params = &fd_p;

  MyFunction fermi;
  fermi.function = &FermiDiracLaguerre;
  fermi.params = &fd_p;

  // test Gauss-Laguerre implementation
  const int n = 32;
  const double alpha = 0.;
  MyQuadrature quad = {.n = n, .dim = 1, .type = kGaulag, .alpha = alpha};
  GaussLaguerre(&quad);

  // Gauss-Laguerre integration (weight = (x-a)^alpha * exp(-b*(x-a)))
  const double a = 0.;
  const double b = 1.;
  double gaulag_result_gsl = GSLLagQuadrature(n, a, b, alpha, &f);
  double gaulag_result = GaussLaguerreIntegrateZeroInf(&quad, &fermi);
  printf("Gauss-Laguerre integration from 0 to inf (GSL):   %.5e\n", gaulag_result_gsl);
  printf("Gauss-Laguerre integration from 0 to inf (our):   %.5e\n", gaulag_result);

  f.function = &FermiDiracLegendre;
  fermi.function = &FermiDiracLegendre;
  // Gauss-Legendre integration from a to b (weight = 1)
  // First version
  double gLeg_1 = GSL_leg_quadr(n, a, b, &f);
  printf("Gauss-Legendre integration from a to b (1): %.5e\n", gLeg_1);

  // Second version
  double gLeg_2 = GSL_leg_integ(n, a, b, &f);
  printf("Gauss-Legendre integration from a to b (2): %.5e\n", gLeg_2);

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = 10.;
  double gLeg_inf = GSL_leg_split(n, &f, s);
  printf("Gauss-Legendre integration from 0 to inf:   %.5e\n", gLeg_inf);
}