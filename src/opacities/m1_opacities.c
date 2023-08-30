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
MyOpacityQuantity M1IntegrandNubar(double nubar, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params, GreyOpacityParams *my_grey_opacity_params) {

  // compute the anti-neutrino distribution function
  double g_t = NuFThick(kH * nubar, &my_grey_opacity_params->distr_pars);
  double g_f = NuFThin(kH * nubar, &my_grey_opacity_params->distr_pars);
  double gbar = my_grey_opacity_params->distr_pars.w_t * g_t + my_grey_opacity_params->distr_pars.w_f * g_f;

  // compute the pair kernels
  MyOpacityQuantity pair_kernels_m1 = PairKernelsM1(my_eos_params, &my_kernel_params->pair_kernel_params);

  double integrand_1_e = nubar * nubar * (1. - gbar) * (pair_kernels_m1.em_e + 0.);
  double integrand_1_x = nubar * nubar * (1. - gbar) * (pair_kernels_m1.em_x + 0.);

  double integrand_2_e = nubar * nubar * gbar * (pair_kernels_m1.abs_e + 0.);
  double integrand_2_x = nubar * nubar * gbar * (pair_kernels_m1.abs_x + 0.);

  MyOpacityQuantity result = {.em_e = integrand_1_e, .em_x = integrand_1_x,
                              .abs_e = integrand_2_e, .abs_x = integrand_2_x};

  return result;
}

void M1Integrand2() {

}

void M1Integrand3() {

}

void M1Integrand4() {

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
