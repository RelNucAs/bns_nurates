//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_integration_special.c
//  \brief test integration routine with special integration routine ('multiple kernel integrator')

#ifdef GSL_INCLUDES_H_

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "gsl_routines/integration_gsl.h"
#include "../../include/bns_nurates.h"
#include "../../include/integration.h"

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
  return exp(-x[2] * x[2]) * exp(-x[1]) * pow(x[0], k) / (exp(x[0] - eta) + 1.);
}

MyKernelQuantity ModifiedFermiDiracLegendreSpecial(double *x, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {
  (void)my_kernel_params;
  struct MyEOSParams *params = (struct MyEOSParams *) my_eos_params;
  int k = (params->nb);
  double eta = (params->mu_e);
  double result = exp(-x[2] * x[2]) * exp(-x[1]) * pow(x[0], k) / (exp(x[0] - eta) + 1.);
  MyKernelQuantity result_special = {.em_x = result, .em_e = result, .abs_x = result, .abs_e = result};

  return result_special;
}

double FermiDiracLegendreGSL(double x, void *p) {
  return ModifiedFermiDiracLegendre(&x, p);
}

int main() {

  printf("===========================================\n");
  printf("Testing 2d integration special routines ...\n");
  printf("===========================================\n");

  double max_error = -42.;

  // function and parameters
  FermiDiracParams fd_p = {.k = 5, .eta = 10.};
  MyEOSParams eos_params = {.nb = 5., .mu_e = 10.};
  MyKernelParams my_kernel_params;

  gsl_function f;
  f.function = &FermiDiracLegendreGSL;
  f.params = &fd_p;

  MyFunctionSpecial fermi;
  fermi.function = &ModifiedFermiDiracLegendreSpecial;
  fermi.eos_params = &eos_params;
  fermi.kernel_params = &my_kernel_params;

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
  quad.z2 = 5. * M_PI;
  quad.nz = 45;
  GaussLegendreMultiD(&quad);

  // Gauss-Legendre integration from 0 to inf (integral split at s)
  double s = 10.;
  double gauleg_inf_gsl = GslLegInfSplit(35, &f, s);
  MyKernelQuantity gauleg_inf = GaussLegendreIntegrateZeroInfSpecial(&quad, &fermi, s);

  printf("Integrate ModifiedFermiSpecial(x) dx dy from x = 0 to inf, y = -1 to 1, z = 0 to 5 pi!\n");
  printf("Result (GSL):   %.5e\n", 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  printf("Result (our):   %.5e, %.5e, %.5e, %.5e\n", gauleg_inf.em_e, gauleg_inf.em_x, gauleg_inf.abs_e, gauleg_inf.abs_x);

  if (fabs(gauleg_inf.em_e - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf.em_e - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  }
  if (fabs(gauleg_inf.em_x - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf.em_x - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  }
  if (fabs(gauleg_inf.abs_e - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf.abs_e - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  }
  if (fabs(gauleg_inf.abs_x - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl) > max_error) {
    max_error = fabs(gauleg_inf.abs_x - 2.08298988126271493703118623899952578251399250678958515806232097 * gauleg_inf_gsl);
  }

  printf("Maximum error: %0.16e\n", max_error);

  if (max_error < 1e-4) {
    return 0;
  } else {
    return 1;
  }

}

#endif //GSL_INCLUDES_H_
