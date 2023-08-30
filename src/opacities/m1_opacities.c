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

/* Computes the integrands from all double integration terms in Leonardo's notes (See Eqn. (50) and (51))
 * Each computation is perfomed for the 'e' & 'x' type neutrinos, so a total of 6 integrands are computed
 *
 * Integrand 1: (4 pi)^2/(hc)^6 nu^3 nubar^2 [R^prod(pair) + R^prod(brem)][1 - g_nubar]
 * Integrand 2: (1/(c J))(4 pi)^2/(hc)^6 nu^3 nubar^2 g_nu [R^prod(pair) + R^prod(brem)][1 - g_nubar]
 * Integrand 3: (1/(c J))(4 pi)^2/(hc)^6 nu^3 nubar^2 g_nu [R^ann(pair) + R^ann(brem)] g_nubar
 */
void M1DoubleIntegrand(double nu, double nubar, GreyOpacityParams *my_grey_opacity_params, double *result) {

  const double four_pi_hc3 = (4. * kPi) / (kH * kH * kH * kClight * kClight * kClight);
  const double four_pi_hc3_sqr_cJ_e = four_pi_hc3 / (kClight * 1.); // @TODO: Add J_e from M1
  const double four_pi_hc3_sqr_cJ_x = four_pi_hc3 / (kClight * 1.); // @TODO: Add J_x from M1

  MyEOSParams my_eos_params = my_grey_opacity_params->eos_pars;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu = TotalNuF(kH * nu, &my_grey_opacity_params->distr_pars);
  double g_nubar = TotalNuF(kH * nubar, &my_grey_opacity_params->distr_pars);

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

  // compute the bremsstrahlung kernels (TODO: Leonardo, check this!)
  my_kernel_params.brem_kernel_params.omega = kH * nu;
  my_kernel_params.brem_kernel_params.omega_prime = kH * nubar;
  MyOpacityQuantity brem_kernels_m1 = BremKernels(&my_kernel_params.brem_kernel_params, &my_eos_params);

  double integrand_1_e = four_pi_hc3 * four_pi_hc3 * nu * nu * nu * nubar * nubar * (1. - g_nubar) * (pair_kernels_m1.em_e + brem_kernels_m1.em_e);
  double integrand_1_x = four_pi_hc3 * four_pi_hc3 * nu * nu * nu * nubar * nubar * (1. - g_nubar) * (pair_kernels_m1.em_x + brem_kernels_m1.em_x);

  double integrand_2_e = integrand_1_e * g_nu * four_pi_hc3_sqr_cJ_e;
  double integrand_2_x = integrand_1_x * g_nu * four_pi_hc3_sqr_cJ_x;

  double integrand_3_e = four_pi_hc3_sqr_cJ_e * four_pi_hc3 * nu * nu * nu * nubar * nubar * g_nu * g_nubar * (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e);
  double integrand_3_x = four_pi_hc3_sqr_cJ_x * four_pi_hc3 * nu * nu * nu * nubar * nubar * g_nu * g_nubar * (pair_kernels_m1.abs_x + brem_kernels_m1.abs_x);

  //double result[6];

  result[0] = integrand_1_e;
  result[1] = integrand_1_x;
  result[2] = integrand_2_e;
  result[3] = integrand_2_x;
  result[4] = integrand_3_e;
  result[5] = integrand_3_x;

}

/* Computes the integrand for all single integrals from Leonardo's notes
 *
 * A total of six integrands are computed
 *
 */
void M1SingleIntegrand(double nu, GreyOpacityParams *my_grey_opacity_params) {

  // compute the neutrino distribution function
  double g_t_nu = NuFThick(nu, &my_grey_opacity_params->distr_pars);
  double g_f_nu = NuFThin(nu, &my_grey_opacity_params->distr_pars);
  double g_nu = my_grey_opacity_params->distr_pars.w_t * g_t_nu + my_grey_opacity_params->distr_pars.w_f * g_f_nu;

  const double four_pi_hc3 = (4. * kPi) / (kH * kH * kH * kClight * kClight * kClight);
  const double four_pi_hc3_cJ_e = four_pi_hc3 / (kClight * 0.); // @TODO: Add J_e from M1
  const double four_pi_hc3_cJ_x = four_pi_hc3 / (kClight * 0.); // @TODO: Add J_x from M1

  double integrand_1_e = 0.;
  double integrand_1_x = 0.;
  double integrand_2_e = 0.;
  double integrand_2_x = 0.;
  double integrand_3_e = four_pi_hc3_cJ_e * 4. * kPi * nu * nu * nu * nu * nu * g_nu * (0.);
  double integrand_3_x = 0.;
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
