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
#include "opacities.h"

/* Compute the 2d integrands for all reactions from Leonardo's notes [Eqns. (50) & (51)]
 * There are a total of two expressions for 'e' and 'x' neutrinos, so 4 integrands in total
 *
 * 1. Contribution to emissivity: (4 pi)^2/(hc)^6 nu^3 nubar^2 [R^prod(pair) + R^prod(brem)][1 - g_nubar]
 * 2. Contribution to absorption coefficient:
 *
 * Note that there are no double integrals for the computation of the scattering coefficient.
 */
MyQuadratureIntegrand M1DoubleIntegrand(double *var, void *p) {

  // energies and parameters
  double nu = var[0];
  double nubar = var[1];
  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *) p;
  MyEOSParams my_eos_params = my_grey_opacity_params->eos_pars;

  // compute some constants
  const double four_pi_hc3 = (4. * kPi) / (kH * kH * kH * kClight * kClight * kClight);
  const double four_pi_hc3_sqr = four_pi_hc3 * four_pi_hc3;

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

  double integrand_1_e = four_pi_hc3_sqr * nu * nu * nu * nubar * nubar * (1. - g_nubar) * (pair_kernels_m1.em_e + brem_kernels_m1.em_e);
  double integrand_1_x = four_pi_hc3_sqr * nu * nu * nu * nubar * nubar * (1. - g_nubar) * (pair_kernels_m1.em_x + brem_kernels_m1.em_x);

  double J_e = my_grey_opacity_params->m1_pars.J;
  double J_x = my_grey_opacity_params->m1_pars.J; // @TODO: check what is to be done for 'e' and 'x' type neutrinos!

  double integrand_2_e = (1./(kClight * J_e)) * (integrand_1_e * g_nu + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_nu * g_nubar);
  double integrand_2_x = (1./(kClight * J_x)) * (integrand_1_x * g_nu + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_x + brem_kernels_m1.abs_x) * g_nu * g_nubar);

  MyQuadratureIntegrand result = {.n = 4};

  result.integrand[0] = integrand_1_e;
  result.integrand[1] = integrand_1_x;
  result.integrand[2] = integrand_2_e;
  result.integrand[3] = integrand_2_x;

  return result;
}

/* Computes the integrand for all single integrals from Leonardo's notes
 *
 * A total of six integrands are computed
 *
 */
void M1SingleIntegrand(double nu, GreyOpacityParams *my_grey_opacity_params) {

  // compute the neutrino & anti-neutrino distribution function
  double g_nu = TotalNuF(kH * nu, &my_grey_opacity_params->distr_pars);

  const double four_pi_hc3 = (4. * kPi) / (kH * kH * kH * kClight * kClight * kClight);
  const double four_pi_hc3_cJ_e = four_pi_hc3 / (kClight * my_grey_opacity_params->m1_pars.J); // @TODO: Add J_e from M1
  const double four_pi_hc3_cJ_x = four_pi_hc3 / (kClight * 0.); // @TODO: Add J_x from M1

  const double iso_scatt = IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

  double integrand_1_e = 0.;
  double integrand_1_x = 0.;
  double integrand_2_e = 0.;
  double integrand_2_x = 0.;
  double integrand_3_e = four_pi_hc3_cJ_e * 4 * kPi * nu * nu * nu * nu * nu * iso_scatt * g_nu;
  double integrand_3_x = integrand_3_e;
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
