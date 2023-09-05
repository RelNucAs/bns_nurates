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

/* Compute the 2d integrands for all reactions from Leonardo's notes [Eqns. (51) & (52)]
 * There are a total of two expressions for 'e' and 'x' neutrinos, so 4 integrands in total
 *
 * 1. Contribution to emissivity: (4 pi)^2/(hc)^6 nu^3 nubar^2 [R^prod(pair) + R^prod(brem)][1 - g_nubar]
 * 2. Contribution to absorption coefficient: (1/(c J)) *(4 pi)^2/(hc)^6 * (nu^3 nubar^2 [R_pro(Pair) + R_pro(Brem)][1 - g_nubar] g_nu
 *                                                                        + nu^3 nubar^2 [R_abs(Pair) + R_abs(Brem)]g_nubar g_nu)
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

  double integrand_2_e = (1. / (kClight * J_e)) * (integrand_1_e * g_nu + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_nu * g_nubar);
  double integrand_2_x = (1. / (kClight * J_x)) * (integrand_1_x * g_nu + nu * nu * nu * nubar * nubar * (pair_kernels_m1.abs_x + brem_kernels_m1.abs_x) * g_nu * g_nubar);

  MyQuadratureIntegrand result = {.n = 4};

  result.integrand[0] = integrand_1_e;
  result.integrand[1] = integrand_1_x;
  result.integrand[2] = integrand_2_e;
  result.integrand[3] = integrand_2_x;

  return result;
}

/* Computes the integrand for all single integrals from Leonardo's notes
 *
 * There are a total of 3 expressions for electron-type neutrinos, electron-type antineutrinos
 * and 'x' neutrinos, so a total of 9 integrands should be computed
 * 
 * However, two of them (those in Eq.(51) and Eq.(52) for 'x' neutrinos) are trivially equal  
 * to zero)
 *
 * 1. Contribution to emissivity: (4 pi /(h c)^3) nu^3 j_x
 * 2. Contribution to absorption coefficient: (1/(c J)) (4 pi /(h c)^3) nu^3 g_nu (j_x + 1/lambda_x)
 * 3. Contribution to scattering coefficient: (1/(c J)) (4 pi)^2 nu^5 g_nu (R_iso(1)/3 - R_iso(0))
 */
MyQuadratureIntegrand M1SingleIntegrand(double *var, void *p) {

  // energy and parameters
  double nu = var[0];
  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *) p;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu = TotalNuF(kH * nu, &my_grey_opacity_params->distr_pars);
  // @TODO: g_nu is not the same for all the neutrino species

  // compute some constants
  const double four_pi_hc3 = (4. * kPi) / (kH * kH * kH * kClight * kClight * kClight);
  const double four_pi_hc3_sqr = four_pi_hc3 * four_pi_hc3;

  const double iso_scatt = IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

  double J_nue = my_grey_opacity_params->m1_pars.J;
  double J_anue = my_grey_opacity_params->m1_pars.J;
  double J_nux = my_grey_opacity_params->m1_pars.J; // @TODO: check what is to be done for 'e' and 'x' type neutrinos!
  // @TODO: M1Quantities structure must be generalized in order to account for different neutrino
  //        species (e.g. double J -> double J[Nsp], where Nsp is the number of different neutrino species)

  MyOpacity abs_em_beta = StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

  double integrand_1_nue = four_pi_hc3 * nu * nu * nu * abs_em_beta.em_nue;
  double integrand_1_anue = four_pi_hc3 * nu * nu * nu * abs_em_beta.em_anue;
  double integrand_2_nue = (1. / (kClight * J_nue)) * four_pi_hc3 * nu * nu * nu * g_nu * abs_em_beta.ab_nue; // @TODO: use the correct g_nu
  double integrand_2_anue = (1. / (kClight * J_anue)) * four_pi_hc3_sqr * nu * nu * nu * g_nu * abs_em_beta.ab_anue; // @TODO: use the correct g_nu
  double integrand_3_nue = (1. / (kClight * J_nue)) * 16. * kPi * kPi * nu * nu * nu * g_nu * iso_scatt; // factor nu^2 already in iso_scatt, @TODO: use the correct g_nu
  double integrand_3_anue = (1. / (kClight * J_anue)) * 16. * kPi * kPi * nu * nu * nu * g_nu * iso_scatt; // @TODO: use the correct g_nu
  double integrand_3_nux  = (1. / (kClight * J_nux)) * 16. * kPi * kPi * nu * nu * nu * g_nu * iso_scatt; // @TODO: use the correct g_nu

  MyQuadratureIntegrand result = {.n = 6};

  result.integrand[0] = integrand_1_nue;
  result.integrand[1] = integrand_1_anue;
  result.integrand[2] = integrand_2_nue;
  result.integrand[3] = integrand_2_anue;
  result.integrand[4] = integrand_3_nue;
  result.integrand[5] = integrand_3_anue;
  //result.integrand... = integrand_3_nux;

  return result;
}

/* Computes the opacities for the M1 code
 *
 */
M1Opacities ComputeM1Opacities(MyQuadrature *quad_1d, MyQuadrature *quad_2d, GreyOpacityParams *my_grey_opacity_params) {

  // set up 1d integration
  MyFunctionMultiD integrand_m1_1d;
  MyQuadratureIntegrand integrand_m1_1d_info = {.n = 6};
  integrand_m1_1d.function = &M1SingleIntegrand;
  integrand_m1_1d.params = my_grey_opacity_params;
  integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

  // set up 2d integration
  MyFunctionMultiD integrand_m1_2d;
  MyQuadratureIntegrand integrand_m1_2d_info = {.n = 4};
  integrand_m1_2d.function = &M1DoubleIntegrand;
  integrand_m1_2d.params = my_grey_opacity_params;
  integrand_m1_2d.my_quadrature_integrand = integrand_m1_2d_info;

  double s = 1. * my_grey_opacity_params->eos_pars.temp; // @TODO: choose this appropriately

  MyQuadratureIntegrand integrals_1d = GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s);
  MyQuadratureIntegrand integrals_2d = GaussLegendreIntegrate2D(quad_2d, &integrand_m1_2d, s);

  M1Opacities m1_opacities;

  m1_opacities.eta_e = integrals_2d.integrand[0] + integrals_1d.integrand[0];
  m1_opacities.eta_x = integrals_2d.integrand[1];
  m1_opacities.kappa_a_e = integrals_2d.integrand[2];
  m1_opacities.kappa_a_x = integrals_2d.integrand[3];
  m1_opacities.kappa_s_e = 0;
  m1_opacities.kappa_s_x = 0;

  return m1_opacities;
}
