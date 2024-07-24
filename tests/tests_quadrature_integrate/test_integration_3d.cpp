//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_3d.c
//  \brief test 3d integration routine

#ifdef GSL_INCLUDES_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "gsl_routines/integration_gsl.hpp"
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

#define M_PI_VAL 3.14159265358979323846264338328 /* pi */

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
  return exp(-x[2]*x[2]) * exp(-x[1]) * pow(x[0], k) / (exp(x[0] - eta) + 1.);
}

double FermiDiracLegendreGSL(double x, void *p) {
  return ModifiedFermiDiracLegendre(&x, p);
}

int main() {

  printf("===================================\n");
  printf("Testing 2d integration routines ...\n");
  printf("===================================\n");

  double max_error = -42.;

  // function and parameters
  FermiDiracParams fd_p = {.k = 5, .eta = 10};

  gsl_function f;
  f.function = &FermiDiracLegendreGSL;
  f.params = &fd_p;

  MyFunction fermi;
  fermi.function = &ModifiedFermiDiracLegendre;
  fermi.params = &fd_p;

  MyQuadrature quad = quadrature_default;
  quad.dim = 3;
  quad.type = kGauleg;
  quad.x1 = 0.;
  quad.x2 = 1.;
  quad.nx = 45;
  quad.y1 = -1.;
  quad.y2 = 1.;
  quad.ny = 45;
  quad.z1 = 0.;
  quad.z2 = 5. * M_PI_VAL;
  quad.nz = 45;
  GaussLegendreMultiD(&quad);

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = 10.;
  double gauleg_inf_gsl = GslLegInfSplit(35, &f, s);
  double gauleg_inf = GaussLegendreIntegrateZeroInf(&quad, &fermi, s);

  printf("Integrate ModifiedFermi(x) dx dy from x = 0 to inf, y = -1 to 1, z = 0 to 5 pi!\n");
  printf("Result (GSL):   %.5e\n", 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  printf("Result (our):   %.5e\n", gauleg_inf);

  if (fabs(gauleg_inf - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  }
  printf("Maximum error: %0.16e\n", max_error);

  if (max_error < 1e-4) {
    return 0;
  } else {
    return 1;
  }

}

#endif //GSL_INCLUDES_H_
