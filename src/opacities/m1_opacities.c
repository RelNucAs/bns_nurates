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

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species];
  double g_nu_bar[total_num_species];

  for (int idx=0; idx<total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    g_nu_bar[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
  }

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

  // compute the bremsstrahlung kernels
  MyKernelQuantity brem_kernels_m1 = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  if (opacity_flags.use_brem) {
    my_kernel_params.brem_kernel_params.omega = nu;
    my_kernel_params.brem_kernel_params.omega_prime = nubar;
    my_kernel_params.brem_kernel_params.l = 0;
    brem_kernels_m1 = BremKernelsLegCoeff(&my_kernel_params.brem_kernel_params, &my_eos_params);
  }

  const double pro_term_nue = (pair_kernels_m1.em_e + brem_kernels_m1.em_e) * (1. - g_nu_bar[id_anue]);
  const double pro_term_anue = (pair_kernels_m1.em_e + brem_kernels_m1.em_e) * (1. - g_nu_bar[id_nue]);
  const double pro_term_nux = (pair_kernels_m1.em_x + brem_kernels_m1.em_x) * (1. - g_nu_bar[id_nux]);

  const double ann_term_nue = (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_nu_bar[id_anue];
  const double ann_term_anue = (pair_kernels_m1.abs_e + brem_kernels_m1.abs_e) * g_nu_bar[id_nue];
  const double ann_term_nux = (pair_kernels_m1.abs_x + brem_kernels_m1.abs_x) * g_nu_bar[id_nux];

  double n_integrand_1_nue = nu * nu * nubar * nubar * pro_term_nue;
  double n_integrand_1_anue = nu * nu * nubar * nubar * pro_term_anue;
  double n_integrand_1_nux = nu * nu * nubar * nubar * pro_term_nux;

  double e_integrand_1_nue = nu * n_integrand_1_nue;
  double e_integrand_1_anue = nu * n_integrand_1_anue;
  double e_integrand_1_nux = nu * n_integrand_1_nux;

  double n_integrand_2_nue = g_nu[id_nue] * nu * nu * nubar * nubar * (pro_term_nue + ann_term_nue);
  double n_integrand_2_anue = g_nu[id_anue] * nu * nu * nubar * nubar * (pro_term_anue + ann_term_anue);
  double n_integrand_2_nux = g_nu[id_nux] * nu * nu * nubar * nubar * (pro_term_nux + ann_term_nux);

  double e_integrand_2_nue = nu * n_integrand_2_nue;
  double e_integrand_2_anue = nu * n_integrand_2_anue;
  double e_integrand_2_nux = nu * n_integrand_2_nux;

  MyQuadratureIntegrand result = {.n = 12};

  result.integrand[0] = n_integrand_1_nue;
  result.integrand[1] = n_integrand_1_anue;
  result.integrand[2] = n_integrand_1_nux;
  result.integrand[3] = n_integrand_2_nue;
  result.integrand[4] = n_integrand_2_anue;
  result.integrand[5] = n_integrand_2_nux;
  result.integrand[6] = e_integrand_1_nue;
  result.integrand[7] = e_integrand_1_anue;
  result.integrand[8] = e_integrand_1_nux;
  result.integrand[9] = e_integrand_2_nue;
  result.integrand[10] = e_integrand_2_anue;
  result.integrand[11] = e_integrand_2_nux;

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
  double g_nu[total_num_species];

  for (int idx=0; idx<total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
  }
  
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
    MyOpacity abs_em_beta = StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars); // [s^-1]

    n_integrand_1_nue = nu * nu * abs_em_beta.em[id_nue]; // [MeV^-1 cm^-3 s^-1]
    n_integrand_1_anue = nu * nu * abs_em_beta.em[id_anue]; // [MeV^-1 cm^-3 s^-1]
    n_integrand_2_nue = nu * nu * g_nu[id_nue] * abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)
    n_integrand_2_anue = nu * nu * g_nu[id_anue] * abs_em_beta.abs[id_anue];

    e_integrand_1_nue = nu * n_integrand_1_nue; // [cm^-3 s^-1]
    e_integrand_1_anue = nu * n_integrand_1_anue; // [cm^-3 s^-1]
    e_integrand_2_nue = nu * n_integrand_2_nue;
    e_integrand_2_anue = nu * n_integrand_2_anue;
  }

  if (opacity_flags.use_iso == 1) {
    const double iso_scatt = IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

    e_integrand_3_nue = 4. * kPi * nu * nu * nu * nu * nu * g_nu[id_nue] * iso_scatt;
    e_integrand_3_anue = 4. * kPi * nu * nu * nu * nu * nu * g_nu[id_anue] * iso_scatt;
    e_integrand_3_nux = 4. * kPi * nu * nu * nu * nu * nu * g_nu[id_nux] * iso_scatt;
  }

  MyQuadratureIntegrand result = {.n = 11};

  result.integrand[0] = n_integrand_1_nue;
  result.integrand[1] = n_integrand_1_anue;
  result.integrand[2] = n_integrand_2_nue;
  result.integrand[3] = n_integrand_2_anue;
  result.integrand[4] = e_integrand_1_nue;
  result.integrand[5] = e_integrand_1_anue;
  result.integrand[6] = e_integrand_2_nue;
  result.integrand[7] = e_integrand_2_anue;
  result.integrand[8] = e_integrand_3_nue;
  result.integrand[9] = e_integrand_3_anue;
  result.integrand[10] = e_integrand_3_nux;

  return result;
}

/* Computes the opacities for the M1 code
 *
 */
M1Opacities ComputeM1Opacities(MyQuadrature *quad_1d, MyQuadrature *quad_2d, GreyOpacityParams *my_grey_opacity_params) {
  // compute some constants
  static const double four_pi_hc3 = (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
  static const double four_pi_hc3_sqr = four_pi_hc3 * four_pi_hc3; // [MeV^-6 cm^-6]

  // set up 1d integration
  MyFunctionMultiD integrand_m1_1d;
  MyQuadratureIntegrand integrand_m1_1d_info = {.n = 11};
  integrand_m1_1d.function = &M1SingleIntegrand;
  integrand_m1_1d.dim = 1;
  integrand_m1_1d.params = my_grey_opacity_params;
  integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

  // set up 2d integration
  MyFunctionMultiD integrand_m1_2d;
  MyQuadratureIntegrand integrand_m1_2d_info = {.n = 12};
  integrand_m1_2d.function = &M1DoubleIntegrand;
  integrand_m1_2d.dim = 2;
  integrand_m1_2d.params = my_grey_opacity_params;
  integrand_m1_2d.my_quadrature_integrand = integrand_m1_2d_info;

  double n[total_num_species];
  double J[total_num_species];

  for (int idx=0; idx<total_num_species; idx++) {
    n[idx] = my_grey_opacity_params->m1_pars.n[idx];
    J[idx] = my_grey_opacity_params->m1_pars.J[idx];
  }

  double s_2d_x[12];
  double s_2d_y[12];

  for (int i=0; i<11; i++) {
    s_2d_x[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
    s_2d_y[i] = s_2d_x[i];
  }
  // @TODO: choose this appropriately
  
  double s_1d[11];
  for (int i=0; i<11; i++) {
    s_1d[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
  }

  MyQuadratureIntegrand integrals_1d = GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_1d);
  MyQuadratureIntegrand integrals_2d = GaussLegendreIntegrate2D(quad_2d, &integrand_m1_2d, s_2d_x, s_2d_y);

  M1Opacities m1_opacities;

  m1_opacities.eta_0[id_nue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[0] + integrals_1d.integrand[0]);
  m1_opacities.eta_0[id_anue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[1] + integrals_1d.integrand[1]);
  m1_opacities.eta_0[id_nux] = four_pi_hc3_sqr * integrals_2d.integrand[2];

  m1_opacities.kappa_0_a[id_nue] = four_pi_hc3 / (kClight * n[id_nue]) * (four_pi_hc3 * integrals_2d.integrand[3] + integrals_1d.integrand[2]);
  m1_opacities.kappa_0_a[id_anue] = four_pi_hc3 / (kClight * n[id_anue]) * (four_pi_hc3 * integrals_2d.integrand[4] + integrals_1d.integrand[3]);
  m1_opacities.kappa_0_a[id_nux] = four_pi_hc3_sqr / (kClight * n[id_nux]) * integrals_2d.integrand[5];

  m1_opacities.eta[id_nue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[6] + integrals_1d.integrand[4]);
  m1_opacities.eta[id_anue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[7] + integrals_1d.integrand[5]);
  m1_opacities.eta[id_nux] = four_pi_hc3_sqr * integrals_2d.integrand[8];

  m1_opacities.kappa_a[id_nue] = four_pi_hc3 / (kClight * J[id_nue]) * (four_pi_hc3 * integrals_2d.integrand[9] + integrals_1d.integrand[6]);
  m1_opacities.kappa_a[id_anue] = four_pi_hc3 / (kClight * J[id_anue]) * (four_pi_hc3 * integrals_2d.integrand[10] + integrals_1d.integrand[7]);
  m1_opacities.kappa_a[id_nux] = four_pi_hc3_sqr / (kClight * J[id_nux]) * integrals_2d.integrand[11];

  m1_opacities.kappa_s[id_nue] = four_pi_hc3 / (kClight * J[id_nue]) * integrals_1d.integrand[8];
  m1_opacities.kappa_s[id_anue] = four_pi_hc3 / (kClight * J[id_anue]) * integrals_1d.integrand[9];
  m1_opacities.kappa_s[id_nux] = four_pi_hc3 / (kClight * J[id_nux]) * integrals_1d.integrand[10];

  return m1_opacities;
}
