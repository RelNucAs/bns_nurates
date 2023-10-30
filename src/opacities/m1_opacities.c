// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file m1_opacities.h
//  \brief compute opacities for all processes in the M1 code

#include "bns_nurates.h"
#include "constants.h"
#include "kernels.h"
#include "functions.h"
#include "integration.h"
#include "distribution.h"
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
  double nu = var[0];     // [MeV]
  double nubar = var[1];  // [MeV]
  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *) p;
  MyEOSParams my_eos_params = my_grey_opacity_params->eos_pars;
  OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

  // compute some constants
  const double four_pi_hc3 = (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
  const double four_pi_hc3_sqr = four_pi_hc3 * four_pi_hc3; // [MeV^-6 cm^-6]

  // compute the neutrino & anti-neutrino distribution function
  double g_nue = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 0);
  double g_anue = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 1);
  double g_nux = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 2);

  double g_nue_bar = TotalNuF(nubar, &my_grey_opacity_params->distr_pars, 0);
  double g_anue_bar = TotalNuF(nubar, &my_grey_opacity_params->distr_pars, 1);
  double g_nux_bar = TotalNuF(nubar, &my_grey_opacity_params->distr_pars, 2);

  // compute the pair kernels
  MyKernelQuantity pair_kernels_m1 = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  MyKernelParams my_kernel_params;
  if (opacity_flags.use_pair) {
    my_kernel_params.pair_kernel_params.omega = nu;
    my_kernel_params.pair_kernel_params.omega_prime = nubar;
    my_kernel_params.pair_kernel_params.cos_theta = 1.;
    my_kernel_params.pair_kernel_params.filter = 0.;
    my_kernel_params.pair_kernel_params.lmax = 0;
    my_kernel_params.pair_kernel_params.mu = 1.;
    my_kernel_params.pair_kernel_params.mu_prime = 1.;
    pair_kernels_m1 = PairKernelsM1(&my_eos_params, &my_kernel_params.pair_kernel_params);
  }

  // compute the bremsstrahlung kernels (TODO: Leonardo, check this!)
  MyKernelQuantity brem_kernels_m1 = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  if (opacity_flags.use_brem) {
    my_kernel_params.brem_kernel_params.omega = nu;
    my_kernel_params.brem_kernel_params.omega_prime = nubar;
    brem_kernels_m1 = BremKernels(&my_kernel_params.brem_kernel_params, &my_eos_params);
  }

  const double pro_term_nue = (pair_kernels_m1.em_e + brem_kernels_m1.em_e) * (1. - g_anue_bar);
  const double pro_term_anue = (pair_kernels_m1.em_e + brem_kernels_m1.em_e) * (1. - g_nue_bar);
  const double pro_term_nux = (pair_kernels_m1.em_x + brem_kernels_m1.em_x) * (1. - g_nux_bar);

  const double ann_term_nue = (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_anue_bar;
  const double ann_term_anue = (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_nue_bar;
  const double ann_term_nux = (pair_kernels_m1.abs_x + brem_kernels_m1.abs_x) * g_nux_bar;

  double n_integrand_1_nue = four_pi_hc3_sqr * nu * nu * nubar * nubar * pro_term_nue;
  double n_integrand_1_anue = four_pi_hc3_sqr * nu * nu * nubar * nubar * pro_term_anue;
  double n_integrand_1_nux = four_pi_hc3_sqr * nu * nu * nubar * nubar * pro_term_nux;

  double e_integrand_1_nue = nu * n_integrand_1_nue;
  double e_integrand_1_anue = nu * n_integrand_1_anue;
  double e_integrand_1_nux = nu * n_integrand_1_nux;
  
  double n_nue = my_grey_opacity_params->m1_pars.n[0];
  double n_anue = my_grey_opacity_params->m1_pars.n[1];  
  double n_nux = my_grey_opacity_params->m1_pars.n[2];

  double J_nue = my_grey_opacity_params->m1_pars.J[0];
  double J_anue = my_grey_opacity_params->m1_pars.J[1];  
  double J_nux = my_grey_opacity_params->m1_pars.J[2];

  double n_integrand_2_nue = four_pi_hc3_sqr * g_nue * nu * nu * nubar * nubar * (pro_term_nue + ann_term_nue);
  double n_integrand_2_anue = four_pi_hc3_sqr * g_anue * nu * nu * nubar * nubar * (pro_term_anue + ann_term_anue);
  double n_integrand_2_nux = four_pi_hc3_sqr * g_nux * nu * nu * nubar * nubar * (pro_term_nux + ann_term_nux);

  double e_integrand_2_nue = nu * n_integrand_2_nue;
  double e_integrand_2_anue = nu * n_integrand_2_anue;
  double e_integrand_2_nux = nu * n_integrand_2_nux;

  MyQuadratureIntegrand result = {.n = 12};

  result.integrand[0] = n_integrand_1_nue;
  result.integrand[1] = n_integrand_1_anue;
  result.integrand[2] = n_integrand_1_nux;
  result.integrand[3] = (1. / (kClight * n_nue)) * n_integrand_2_nue;
  result.integrand[4] = (1. / (kClight * n_anue)) * n_integrand_2_anue;
  result.integrand[5] = (1. / (kClight * n_nux)) * n_integrand_2_nux;
  result.integrand[6] = e_integrand_1_nue;
  result.integrand[7] = e_integrand_1_anue;
  result.integrand[8] = e_integrand_1_nux;
  result.integrand[9] = (1. / (kClight * J_nue)) * e_integrand_2_nue;
  result.integrand[10] = (1. / (kClight * J_anue)) * e_integrand_2_anue;
  result.integrand[11] = (1. / (kClight * J_nux)) * e_integrand_2_nux;

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
  double nu = var[0]; // [MeV]
  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *) p;
  OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

  // compute the neutrino & anti-neutrino distribution function
  double g_nue = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 0);
  double g_anue = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 1);
  double g_nux = TotalNuF(nu, &my_grey_opacity_params->distr_pars, 2);

  // compute some constants
  const double four_pi_hc3 = (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]

  const double iso_scatt = IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

  const double n_nue = my_grey_opacity_params->m1_pars.n[0];
  const double n_anue = my_grey_opacity_params->m1_pars.n[1];
  const double n_nux = my_grey_opacity_params->m1_pars.n[2];

  double J_nue = my_grey_opacity_params->m1_pars.J[0];
  double J_anue = my_grey_opacity_params->m1_pars.J[1];
  double J_nux = my_grey_opacity_params->m1_pars.J[2];

  MyOpacity abs_em_beta = StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars); // [s^-1]

  double n_integrand_1_nue = 0.;
  double n_integrand_1_anue = 0.;
  double n_integrand_2_nue = 0.;
  double n_integrand_2_anue = 0.;

  double e_integrand_1_nue = 0.;
  double e_integrand_1_anue = 0.;
  double e_integrand_2_nue = 0.;
  double e_integrand_2_anue = 0.;
  double e_integrand_3_nue = 0.;
  double e_integrand_3_anue = 0.;
  double e_integrand_3_nux = 0.;

  if (opacity_flags.use_abs_em == 1) {
    n_integrand_1_nue = four_pi_hc3 * nu * nu * abs_em_beta.em_nue; // [MeV^-1 cm^-3 s^-1]
    n_integrand_1_anue = four_pi_hc3 * nu * nu * abs_em_beta.em_anue; // [MeV^-1 cm^-3 s^-1]
    n_integrand_2_nue = four_pi_hc3 * nu * nu * g_nue * abs_em_beta.ab_nue; // ab = em + ab (stimulated absorption)
    n_integrand_2_anue = four_pi_hc3 * nu * nu * g_anue * abs_em_beta.ab_anue;

    e_integrand_1_nue = nu * n_integrand_1_nue; // [cm^-3 s^-1]
    e_integrand_1_anue = nu * n_integrand_1_anue; // [cm^-3 s^-1]
    e_integrand_2_nue = nu * n_integrand_2_nue;
    e_integrand_2_anue = nu * n_integrand_2_anue;
  }

  if (opacity_flags.use_iso == 1) {
    e_integrand_3_nue = (4. * kPi / (kClight * J_nue)) * 4. * kPi * nu * nu * nu * nu * nu * g_nue * iso_scatt;
    e_integrand_3_anue = (4. * kPi / (kClight * J_anue)) * 4. * kPi * nu * nu * nu * nu * nu * g_anue * iso_scatt;
    e_integrand_3_nux = (4. * kPi / (kClight * J_nux)) * 4. * kPi * nu * nu * nu * nu * nu * g_nux * iso_scatt;
  }

  MyQuadratureIntegrand result = {.n = 11};

  result.integrand[0] = n_integrand_1_nue;
  result.integrand[1] = n_integrand_1_anue;
  result.integrand[2] = (1. / (kClight * n_nue)) * n_integrand_2_nue;
  result.integrand[3] = (1. / (kClight * n_anue)) * n_integrand_2_anue;
  result.integrand[4] = e_integrand_1_nue;
  result.integrand[5] = e_integrand_1_anue;
  result.integrand[6] = (1. / (kClight * J_nue)) * e_integrand_2_nue;
  result.integrand[7] = (1. / (kClight * J_anue)) * e_integrand_2_anue;
  result.integrand[8] = e_integrand_3_nue;
  result.integrand[9] = e_integrand_3_anue;
  result.integrand[10] = e_integrand_3_nux;

  return result;
}

/* Computes the opacities for the M1 code
 *
 */
M1Opacities ComputeM1Opacities(MyQuadrature *quad_1d, MyQuadrature *quad_2d, GreyOpacityParams *my_grey_opacity_params) {

  // set up 1d integration
  MyFunctionMultiD integrand_m1_1d;
  MyQuadratureIntegrand integrand_m1_1d_info = {.n = 7};
  integrand_m1_1d.function = &M1SingleIntegrand;
  integrand_m1_1d.dim = 1;
  integrand_m1_1d.params = my_grey_opacity_params;
  integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

  // set up 2d integration
  MyFunctionMultiD integrand_m1_2d;
  MyQuadratureIntegrand integrand_m1_2d_info = {.n = 6};
  integrand_m1_2d.function = &M1DoubleIntegrand;
  integrand_m1_2d.dim = 2;
  integrand_m1_2d.params = my_grey_opacity_params;
  integrand_m1_2d.my_quadrature_integrand = integrand_m1_2d_info;

  double s = 1.5 * my_grey_opacity_params->eos_pars.temp; // @TODO: choose this appropriately
  //double s = 1;

  MyQuadratureIntegrand integrals_1d = GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s);
  MyQuadratureIntegrand integrals_2d = GaussLegendreIntegrate2D(quad_2d, &integrand_m1_2d, s);

  M1Opacities m1_opacities;

  m1_opacities.eta_0_nue = integrals_2d.integrand[0] + integrals_1d.integrand[0];
  m1_opacities.eta_0_anue = integrals_2d.integrand[1] + integrals_1d.integrand[1];
  m1_opacities.eta_0_nux = integrals_2d.integrand[2];
  m1_opacities.kappa_0_a_nue = integrals_2d.integrand[3] + integrals_1d.integrand[2];
  m1_opacities.kappa_0_a_anue = integrals_2d.integrand[4] + integrals_1d.integrand[3];
  m1_opacities.kappa_0_a_nux = integrals_2d.integrand[5];

  m1_opacities.eta_nue = integrals_2d.integrand[6] + integrals_1d.integrand[4];
  m1_opacities.eta_anue = integrals_2d.integrand[7] + integrals_1d.integrand[5];
  m1_opacities.eta_nux = integrals_2d.integrand[8];
  m1_opacities.kappa_a_nue = integrals_2d.integrand[9] + integrals_1d.integrand[6];
  m1_opacities.kappa_a_anue = integrals_2d.integrand[10] + integrals_1d.integrand[7];
  m1_opacities.kappa_a_nux = integrals_2d.integrand[11];
  m1_opacities.kappa_s_nue = integrals_1d.integrand[8];
  m1_opacities.kappa_s_anue = integrals_1d.integrand[9];
  m1_opacities.kappa_s_nux = integrals_1d.integrand[10];

  return m1_opacities;
}
