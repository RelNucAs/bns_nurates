//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_1d.c
//  \brief test 1d integration routine

#ifdef GSL_INCLUDES_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "gsl_routines/integration_gsl.hpp"
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

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
double ModifiedFermiDiracLegendre(double *x, void *p) {
  struct FermiDiracParams *params = (struct FermiDiracParams *) p;
  int k = (params->k);
  double eta = (params->eta);
  return pow(x[0], k) / (exp(x[0] - eta) + 1.);
}

double FermiDiracLegendreGSL(double x, void *p) {
  return ModifiedFermiDiracLegendre(&x, p);
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
  (void)p;
  return 1. / (x * x * exp(x));
}

int main() {

  printf("===================================\n");
  printf("Testing 1d integration routines ...\n");
  printf("===================================\n");

  double max_error = -42.;

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
  MyQuadrature quad = quadrature_default;
  quad.nx = n;
  quad.dim = 1;
  quad.type = kGaulag;
  quad.x1 = 0.;
  quad.x2 = 1.;
  quad.alpha = alpha;
  GaussLaguerre(&quad);

  // Gauss-Laguerre integration (weight = (points-a)^alpha * exp(-b*(points-a)))
  const double a = 0.;
  const double b = 1.;
  double gaulag_result_gsl = GSLLagQuadrature(n, a, b, alpha, &f);
  double gaulag_result = GaussLaguerreIntegrateZeroInf(&quad, &fermi);

  printf("Integrate Fermi-Dirac function ...\n");
  printf("\n");
  printf("Gauss-Laguerre integration from 0 to inf (GSL):   %.5e\n", gaulag_result_gsl);
  printf("Gauss-Laguerre integration from 0 to inf (our):   %.5e\n", gaulag_result);
  printf("Max error for Gauss-Laguerre result: %0.16e\n", fabs(gaulag_result - gaulag_result_gsl));

  if (fabs(gaulag_result - gaulag_result_gsl) > max_error) {
    max_error = fabs(gaulag_result - gaulag_result_gsl);
  }

  // Gauss-Legendre integration tests
  f.function = &FermiDiracLegendreGSL;
  fermi.function = &ModifiedFermiDiracLegendre;

  // Gauss-Legendre integration from 0 to 1 (weight = 1)
  // first version
  double gauleg_result_gsl_1 = GSLLegQuadrature(n, a, b, &f);
  printf("\n");
  printf("Gauss-Legendre integration from 0 to 1 (GSL type 1): %.5e\n", gauleg_result_gsl_1);
  // second version
  double gauleg_result_gsl_2 = GSLLegIntegrate(n, a, b, &f);
  printf("Gauss-Legendre integration from 0 to 1 (GSL type 2): %.5e\n", gauleg_result_gsl_2);
  // our version
  quad = quadrature_default;
  quad.type = kGauleg;
  GaussLegendreMultiD(&quad);
  double gauleg_result = 0.;
  for (int i = 0; i < quad.nx; i++) {
    gauleg_result += ModifiedFermiDiracLegendre(&quad.points[i], &fd_p) * quad.w[i];
  }
  printf("Gauss-Legendre integration from 0 to 1 (   ours   ): %.5e\n", gauleg_result);
  printf("Maximum error between our result & GSL type 1: %0.16e\n", fabs(gauleg_result - gauleg_result_gsl_1));
  printf("Maximum error between our result & GSL type 2: %0.16e\n", fabs(gauleg_result - gauleg_result_gsl_2));

  if (fabs(gauleg_result - gauleg_result_gsl_1) > max_error) {
    max_error = fabs(gauleg_result - gauleg_result_gsl_1);
  }
  if (fabs(gauleg_result - gauleg_result_gsl_2) > max_error) {
    max_error = fabs(gauleg_result - gauleg_result_gsl_2);
  }

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = fmax(fd_p.k,fd_p.eta); //10.;
  double gauleg_inf_gsl = GslLegInfSplit(n, &f, s);
  double gauleg_inf = GaussLegendreIntegrateZeroInf(&quad, &fermi, s);
  printf("\n");
  printf("Gauss-Legendre integration of Fermi function from 0 to inf (GSL):   %.5e\n", gauleg_inf_gsl);
  printf("Gauss-Legendre integration of Fermi function from 0 to inf (our):   %.5e\n", gauleg_inf);
  printf("Maximum error: %0.16e\n", fabs(gauleg_inf - gauleg_inf_gsl));

  if (fabs(gauleg_inf - gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf - gauleg_inf_gsl);
  }

  if (max_error < 1e-5) {
    return 0;
  } else {
    return 1;
  }

}

#endif //GSL_INCLUDES_H_
