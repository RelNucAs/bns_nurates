// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file m1_opacities.h
//  \brief compute opacities for all processes in the M1 code

#include "../bns_nurates.h"
#include "../constants.h"
#include "kernels/kernels.h"
#include "../functions/functions.h"
#include "../integration/integration.h"
#include "../distribution/distribution.h"

/* Computes the integrand for nubar integration in Eqn (50) of Leonardo's notes
 *
 * I1 = nubar^2 (R^P_{pair} + R^p_{brem})(1 - gbar)
 *
 * where gbar = w_t g_t + w_f g_f from Federico's notes
 *
 * Returns two doubles, one for electron neutrino and another for mu/tau neutrino
 *
 */
void M1DoubleIntegrand(double nu, double nubar, GreyOpacityParams *my_grey_opacity_params) {

  MyEOSParams my_eos_params = my_grey_opacity_params->eos_pars;

  // compute the anti-neutrino distribution function
  double g_t = NuFThick(nubar, &my_grey_opacity_params->distr_pars);
  double g_f = NuFThin(nubar, &my_grey_opacity_params->distr_pars);
  double gbar = my_grey_opacity_params->distr_pars.w_t * g_t + my_grey_opacity_params->distr_pars.w_f * g_f;

  // compute the pair kernels
  MyKernelParams my_kernel_params;
  my_kernel_params.pair_kernel_params.omega = kH * nu;
  my_kernel_params.pair_kernel_params.omega_prime = kH * nubar;
  my_kernel_params.pair_kernel_params.cos_theta = 1.;
  my_kernel_params.pair_kernel_params.filter = 0.;
  my_kernel_params.pair_kernel_params.lmax = 0;
  my_kernel_params.pair_kernel_params.mu = 1.;
  my_kernel_params.pair_kernel_params.mu_prime = 1.;

  MyOpacityQuantity pair_kernels_m1 = PairKernelsM1(&my_eos_params, &my_kernel_params.pair_kernel_params);

  // compute the bremsstrahlung kernels

  double integrand_1_e = nubar * nubar * (1. - gbar) * (pair_kernels_m1.em_e + 0.);
  double integrand_1_x = nubar * nubar * (1. - gbar) * (pair_kernels_m1.em_x + 0.);

  double integrand_2_e = nubar * nubar * gbar * (pair_kernels_m1.abs_e + 0.);
  double integrand_2_x = nubar * nubar * gbar * (pair_kernels_m1.abs_x + 0.);

  double integrand_eta_e = nu * nu * nu * integrand_1_e;
  double integrand_eta_x = nu * nu * nu * integrand_1_x;

  double integrand_kappa_a_e = integrand_eta_e * gbar + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_e + 0.);
  double integrand_kappa_a_x = integrand_eta_x * gbar + + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_x + 0.);

}

void M1SingleIntegrand(double nu, GreyOpacityParams *my_grey_opacity_params) {

}


/* Integrand for the M1 opacities
 *
 */
M1Opacities M1OpacitiesIntegrand(var3d x, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params, GreyOpacityParams *my_grey_opacity_params) {

}

/* Computes the M1 opacities for the M1 code
 *
 */
M1Opacities ComputeM1Opacities(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params, GreyOpacityParams *my_grey_opacity_params) {

}
