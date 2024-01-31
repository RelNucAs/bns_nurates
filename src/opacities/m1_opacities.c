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
MyQuadratureIntegrand M1CoeffsDoubleIntegrand(double *var, void *p) {

  // energies and parameters
  double nu = var[0];     // [MeV]
  double nu_bar = var[1]; // [MeV]

  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *)p;
  OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species];
  double g_nu_bar[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    g_nu_bar[idx] = TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
  }

  // compute the pair kernels
  MyKernelQuantity pair_kernels_m1 = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  if (opacity_flags.use_pair) {
    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega = nu;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega_prime = nu_bar;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = 1.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.filter = 0.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax = 0;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.mu = 1.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime = 1.;
    // pair_kernels_m1 = PairKernelsM1(&my_grey_opacity_params->eos_pars, &my_grey_opacity_params->kernel_pars.pair_kernel_params);
    pair_kernels_m1 = PairKernelsM1Optimized(&my_grey_opacity_params->eos_pars, &my_grey_opacity_params->kernel_pars.pair_kernel_params);
  }

  // compute the bremsstrahlung kernels
  MyKernelOutput brem_kernels_m1 = {0};
  if (opacity_flags.use_brem) {
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega = nu;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega_prime = nu_bar;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
    brem_kernels_m1 = BremKernelsLegCoeff(&my_grey_opacity_params->kernel_pars.brem_kernel_params, &my_grey_opacity_params->eos_pars);
  }

  // compute the inelastic NES/NPS kernels
  MyKernelOutput inelastic_kernels_m1 = {0};
  if (opacity_flags.use_inelastic_scatt) {
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
    inelastic_kernels_m1 = InelasticScattKernels(&my_grey_opacity_params->kernel_pars.inelastic_kernel_params, &my_grey_opacity_params->eos_pars);
  }

  double pro_term[total_num_species];

  pro_term[id_nue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_nue]) * (1. - g_nu_bar[id_anue]) + inelastic_kernels_m1.em[id_nue] * g_nu_bar[id_nue];
  pro_term[id_anue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_anue]) * (1. - g_nu_bar[id_nue]) + inelastic_kernels_m1.em[id_anue] * g_nu_bar[id_anue];
  pro_term[id_nux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_nux]) * (1. - g_nu_bar[id_anux]) + inelastic_kernels_m1.em[id_nux] * g_nu_bar[id_nux];
  pro_term[id_anux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_anux]) * (1. - g_nu_bar[id_nux]) + inelastic_kernels_m1.em[id_anux] * g_nu_bar[id_anux];

  // pro_term[id_nue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_nue]) + inelastic_kernels_m1.em[id_nue] * g_nu_bar[id_nue];
  // pro_term[id_anue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_anue]) + inelastic_kernels_m1.em[id_anue] * g_nu_bar[id_anue];
  // pro_term[id_nux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_nux]) + inelastic_kernels_m1.em[id_nux] * g_nu_bar[id_nux];
  // pro_term[id_anux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_anux]) + inelastic_kernels_m1.em[id_anux] * g_nu_bar[id_anux];

  double ann_term[total_num_species];

  ann_term[id_nue] = (pair_kernels_m1.abs_e + brem_kernels_m1.abs[id_nue]) * g_nu_bar[id_anue] + inelastic_kernels_m1.abs[id_nue] * (1. - g_nu_bar[id_nue]);
  ann_term[id_anue] = (pair_kernels_m1.abs_e + brem_kernels_m1.abs[id_anue]) * g_nu_bar[id_nue] + inelastic_kernels_m1.abs[id_anue] * (1. - g_nu_bar[id_anue]);
  ann_term[id_nux] = (pair_kernels_m1.abs_x + brem_kernels_m1.abs[id_nux]) * g_nu_bar[id_anux] + inelastic_kernels_m1.abs[id_nux] * (1. - g_nu_bar[id_nux]);
  ann_term[id_anux] = (pair_kernels_m1.abs_x + brem_kernels_m1.abs[id_anux]) * g_nu_bar[id_nux] + inelastic_kernels_m1.abs[id_anux] * (1. - g_nu_bar[id_anux]);

  double n_integrand_1[total_num_species], e_integrand_1[total_num_species];
  double n_integrand_2[total_num_species], e_integrand_2[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    n_integrand_1[idx] = nu * nu * nu_bar * nu_bar * pro_term[idx];
    e_integrand_1[idx] = nu * n_integrand_1[idx];

    n_integrand_2[idx] = g_nu[idx] * nu * nu * nu_bar * nu_bar * (pro_term[idx] + ann_term[idx]);
    // n_integrand_2[idx] = g_nu[idx] * nu * nu * nu_bar * nu_bar * ann_term[idx];
    e_integrand_2[idx] = nu * n_integrand_2[idx];
  }

  MyQuadratureIntegrand result = {.n = 16};

  result.integrand[0] = n_integrand_1[id_nue];
  result.integrand[1] = n_integrand_1[id_anue];
  result.integrand[2] = n_integrand_1[id_nux];
  result.integrand[3] = n_integrand_1[id_anux];
  result.integrand[4] = n_integrand_2[id_nue];
  result.integrand[5] = n_integrand_2[id_anue];
  result.integrand[6] = n_integrand_2[id_nux];
  result.integrand[7] = n_integrand_2[id_anux];
  result.integrand[8] = e_integrand_1[id_nue];
  result.integrand[9] = e_integrand_1[id_anue];
  result.integrand[10] = e_integrand_1[id_nux];
  result.integrand[11] = e_integrand_1[id_anux];
  result.integrand[12] = e_integrand_2[id_nue];
  result.integrand[13] = e_integrand_2[id_anue];
  result.integrand[14] = e_integrand_2[id_nux];
  result.integrand[15] = e_integrand_2[id_anux];

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
MyQuadratureIntegrand M1CoeffsSingleIntegrand(double *var, void *p) {

  // energy and parameters
  double nu = var[0]; // [MeV]

  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *)p;
  OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
  }

  double n_integrand_1[total_num_species] = {0.}, e_integrand_1[total_num_species] = {0.};
  double n_integrand_2[total_num_species] = {0.}, e_integrand_2[total_num_species] = {0.};
  double e_integrand_3[total_num_species] = {0.};

  if (opacity_flags.use_abs_em == 1) {
    MyOpacity abs_em_beta = StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars); // [s^-1]
    // MyOpacity abs_em_beta = AbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars); // [s^-1]

    for (int idx = 0; idx < total_num_species; idx++) {
      n_integrand_1[idx] = nu * nu * abs_em_beta.em[idx]; // [MeV^-1 cm^-3 s^-1]
      e_integrand_1[idx] = nu * n_integrand_1[idx];       // [cm^-3 s^-1]

      n_integrand_2[idx] = nu * nu * g_nu[idx] * abs_em_beta.abs[idx]; // ab = em + ab (stimulated absorption)
      e_integrand_2[idx] = nu * n_integrand_2[idx];
    }
  }

  if (opacity_flags.use_iso == 1) {
    const double iso_scatt = IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars);

    for (int idx = 0; idx < total_num_species; idx++) {
      e_integrand_3[idx] = 4. * kPi * nu * nu * nu * nu * nu * g_nu[idx] * iso_scatt;
    }
  }

  MyQuadratureIntegrand result = {.n = 12};

  result.integrand[0] = n_integrand_1[id_nue];
  result.integrand[1] = n_integrand_1[id_anue];
  result.integrand[2] = n_integrand_2[id_nue];
  result.integrand[3] = n_integrand_2[id_anue];
  result.integrand[4] = e_integrand_1[id_nue];
  result.integrand[5] = e_integrand_1[id_anue];
  result.integrand[6] = e_integrand_2[id_nue];
  result.integrand[7] = e_integrand_2[id_anue];
  result.integrand[8] = e_integrand_3[id_nue];
  result.integrand[9] = e_integrand_3[id_anue];
  result.integrand[10] = e_integrand_3[id_nux];
  result.integrand[11] = e_integrand_3[id_anux];

  return result;
}

void ComputeM1OpacitiesTest(MyQuadrature *quad_1d, GreyOpacityParams *my_grey_opacity_params, double t, M1Matrix *out_n, M1Matrix *out_j) {

  double nu, nu_bar;

  const int n = quad_1d->nx;

  double pro_term_ij[total_num_species], pro_term_ji[total_num_species];
  double ann_term_ij[total_num_species], ann_term_ji[total_num_species];

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species], g_nu_bar[total_num_species];

  M1Matrix beta, pair, brem, inel;

  BetaOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars, &my_grey_opacity_params->opacity_pars, t, &beta);

  PairOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars, &my_grey_opacity_params->kernel_pars, t, &pair);

  my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
  BremOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars, &my_grey_opacity_params->kernel_pars, t, &brem);

  InelasticOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars, &my_grey_opacity_params->kernel_pars, t, &inel);


  for (int idx = 0; idx < total_num_species; idx++) {
    out_j->m1_mat_ab[idx] = (double **)malloc(sizeof(double *) * 2 * quad_1d->nx);
    out_j->m1_mat_em[idx] = (double **)malloc(sizeof(double *) * 2 * quad_1d->nx);

    out_n->m1_mat_ab[idx] = (double **)malloc(sizeof(double *) * 2 * quad_1d->nx);
    out_n->m1_mat_em[idx] = (double **)malloc(sizeof(double *) * 2 * quad_1d->nx);

    for (int i = 0; i < 2 * n; i++) {
      out_j->m1_mat_ab[idx][i] = (double *)malloc(sizeof(double) * 2 * quad_1d->nx);
      out_j->m1_mat_em[idx][i] = (double *)malloc(sizeof(double) * 2 * quad_1d->nx);

      out_n->m1_mat_ab[idx][i] = (double *)malloc(sizeof(double) * 2 * quad_1d->nx);
      out_n->m1_mat_em[idx][i] = (double *)malloc(sizeof(double) * 2 * quad_1d->nx);
    }
  }

  for (int i = 0; i < n; i++) {

    for (int j = i; j < n; j++) {

      // energies and parameters
      nu = t * quad_1d->points[i];
      nu_bar = t * quad_1d->points[j];

      for (int idx = 0; idx < total_num_species; idx++) {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
        g_nu_bar[idx] = TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
      }

      pro_term_ij[id_nue] = (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) * (1. - g_nu_bar[id_anue]) + inel.m1_mat_em[id_nue][i][j] * g_nu_bar[id_nue];
      pro_term_ij[id_anue] = (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) * (1. - g_nu_bar[id_nue]) + inel.m1_mat_em[id_anue][i][j] * g_nu_bar[id_anue];
      pro_term_ij[id_nux] = (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) * (1. - g_nu_bar[id_anux]) + inel.m1_mat_em[id_nux][i][j] * g_nu_bar[id_nux];
      pro_term_ij[id_anux] = (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) * (1. - g_nu_bar[id_nux]) + inel.m1_mat_em[id_anux][i][j] * g_nu_bar[id_anux];

      pro_term_ji[id_nue] = (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) * (1. - g_nu[id_anue]) + inel.m1_mat_em[id_nue][j][i] * g_nu[id_nue];
      pro_term_ji[id_anue] = (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) * (1. - g_nu[id_nue]) + inel.m1_mat_em[id_anue][j][i] * g_nu[id_anue];
      pro_term_ji[id_nux] = (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) * (1. - g_nu[id_anux]) + inel.m1_mat_em[id_nux][j][i] * g_nu[id_nux];
      pro_term_ji[id_anux] = (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) * (1. - g_nu[id_nux]) + inel.m1_mat_em[id_anux][j][i] * g_nu[id_anux];

      ann_term_ij[id_nue] = (pair.m1_mat_ab[id_nue][i][j] + brem.m1_mat_ab[0][i][j]) * g_nu_bar[id_anue] + inel.m1_mat_ab[id_nue][i][j] * (1. - g_nu_bar[id_nue]);
      ann_term_ij[id_anue] = (pair.m1_mat_ab[id_anue][i][j] + brem.m1_mat_ab[0][i][j]) * g_nu_bar[id_nue] + inel.m1_mat_ab[id_anue][i][j] * (1. - g_nu_bar[id_anue]);
      ann_term_ij[id_nux] = (pair.m1_mat_ab[id_nux][i][j] + brem.m1_mat_ab[0][i][j]) * g_nu_bar[id_anux] + inel.m1_mat_ab[id_nux][i][j] * (1. - g_nu_bar[id_nux]);
      ann_term_ij[id_anux] = (pair.m1_mat_ab[id_anux][i][j] + brem.m1_mat_ab[0][i][j]) * g_nu_bar[id_nux] + inel.m1_mat_ab[id_anux][i][j] * (1. - g_nu_bar[id_anux]);

      ann_term_ji[id_nue] = (pair.m1_mat_ab[id_nue][j][i] + brem.m1_mat_ab[0][j][i]) * g_nu[id_anue] + inel.m1_mat_ab[id_nue][j][i] * (1. - g_nu[id_nue]);
      ann_term_ji[id_anue] = (pair.m1_mat_ab[id_anue][j][i] + brem.m1_mat_ab[0][j][i]) * g_nu[id_nue] + inel.m1_mat_ab[id_anue][j][i] * (1. - g_nu[id_anue]);
      ann_term_ji[id_nux] = (pair.m1_mat_ab[id_nux][j][i] + brem.m1_mat_ab[0][j][i]) * g_nu[id_anux] + inel.m1_mat_ab[id_nux][j][i] * (1. - g_nu[id_nux]);
      ann_term_ji[id_anux] = (pair.m1_mat_ab[id_anux][j][i] + brem.m1_mat_ab[0][j][i]) * g_nu[id_nux] + inel.m1_mat_ab[id_anux][j][i] * (1. - g_nu[id_anux]);

      for (int idx = 0; idx < total_num_species; idx ++) {
        out_n->m1_mat_em[idx][i][j] = nu * nu * nu_bar * nu_bar * pro_term_ij[idx];
        out_n->m1_mat_em[idx][j][i] = nu_bar * nu_bar * nu * nu * pro_term_ji[idx];

        out_n->m1_mat_ab[idx][i][j] = nu * nu * g_nu[idx] * nu_bar * nu_bar * (pro_term_ij[idx] + ann_term_ij[idx]);
        out_n->m1_mat_ab[idx][j][i] = nu_bar * nu_bar * g_nu_bar[idx] * nu * nu * (pro_term_ji[idx] + ann_term_ji[idx]);

        out_j->m1_mat_em[idx][i][j] = nu * out_n->m1_mat_em[idx][i][j];
        out_j->m1_mat_em[idx][j][i] = nu_bar * out_n->m1_mat_em[idx][j][i];

        out_j->m1_mat_ab[idx][i][j] = nu * out_n->m1_mat_ab[idx][i][j];
        out_j->m1_mat_ab[idx][j][i] = nu_bar * out_n->m1_mat_ab[idx][j][i];
      }

     
      // energies and parameters
      nu = t / quad_1d->points[i];
      nu_bar = t / quad_1d->points[j];

      for (int idx = 0; idx < total_num_species; idx++) {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
        g_nu_bar[idx] = TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
      }

      pro_term_ij[id_nue] = (pair.m1_mat_em[id_nue][n+i][n+j] + brem.m1_mat_em[0][n+i][n+j]) * (1. - g_nu_bar[id_anue]) + inel.m1_mat_em[id_nue][n+i][n+j] * g_nu_bar[id_nue];
      pro_term_ij[id_anue] = (pair.m1_mat_em[id_anue][n+i][n+j] + brem.m1_mat_em[0][n+i][n+j]) * (1. - g_nu_bar[id_nue]) + inel.m1_mat_em[id_anue][n+i][n+j] * g_nu_bar[id_anue];
      pro_term_ij[id_nux] = (pair.m1_mat_em[id_nux][n+i][n+j] + brem.m1_mat_em[0][n+i][n+j]) * (1. - g_nu_bar[id_anux]) + inel.m1_mat_em[id_nux][n+i][n+j] * g_nu_bar[id_nux];
      pro_term_ij[id_anux] = (pair.m1_mat_em[id_anux][n+i][n+j] + brem.m1_mat_em[0][n+i][n+j]) * (1. - g_nu_bar[id_nux]) + inel.m1_mat_em[id_anux][n+i][n+j] * g_nu_bar[id_anux];

      pro_term_ji[id_nue] = (pair.m1_mat_em[id_nue][n+j][n+i] + brem.m1_mat_em[0][n+j][n+i]) * (1. - g_nu[id_anue]) + inel.m1_mat_em[id_nue][n+j][n+i] * g_nu[id_nue];
      pro_term_ji[id_anue] = (pair.m1_mat_em[id_anue][n+j][n+i] + brem.m1_mat_em[0][n+j][n+i]) * (1. - g_nu[id_nue]) + inel.m1_mat_em[id_anue][n+j][n+i] * g_nu[id_anue];
      pro_term_ji[id_nux] = (pair.m1_mat_em[id_nux][n+j][n+i] + brem.m1_mat_em[0][n+j][n+i]) * (1. - g_nu[id_anux]) + inel.m1_mat_em[id_nux][n+j][n+i] * g_nu[id_nux];
      pro_term_ji[id_anux] = (pair.m1_mat_em[id_anux][n+j][n+i] + brem.m1_mat_em[0][n+j][n+i]) * (1. - g_nu[id_nux]) + inel.m1_mat_em[id_anux][n+j][n+i] * g_nu[id_anux];

      ann_term_ij[id_nue] = (pair.m1_mat_ab[id_nue][n+i][n+j] + brem.m1_mat_ab[0][n+i][n+j]) * g_nu_bar[id_anue] + inel.m1_mat_ab[id_nue][n+i][n+j] * (1. - g_nu_bar[id_nue]);
      ann_term_ij[id_anue] = (pair.m1_mat_ab[id_anue][n+i][n+j] + brem.m1_mat_ab[0][n+i][n+j]) * g_nu_bar[id_nue] + inel.m1_mat_ab[id_anue][n+i][n+j] * (1. - g_nu_bar[id_anue]);
      ann_term_ij[id_nux] = (pair.m1_mat_ab[id_nux][n+i][n+j] + brem.m1_mat_ab[0][n+i][n+j]) * g_nu_bar[id_anux] + inel.m1_mat_ab[id_nux][n+i][n+j] * (1. - g_nu_bar[id_nux]);
      ann_term_ij[id_anux] = (pair.m1_mat_ab[id_anux][n+i][n+j] + brem.m1_mat_ab[0][n+i][n+j]) * g_nu_bar[id_nux] + inel.m1_mat_ab[id_anux][n+i][n+j] * (1. - g_nu_bar[id_anux]);

      ann_term_ji[id_nue] = (pair.m1_mat_ab[id_nue][n+j][n+i] + brem.m1_mat_ab[0][n+j][n+i]) * g_nu[id_anue] + inel.m1_mat_ab[id_nue][n+j][n+i] * (1. - g_nu[id_nue]);
      ann_term_ji[id_anue] = (pair.m1_mat_ab[id_anue][n+j][n+i] + brem.m1_mat_ab[0][n+j][n+i]) * g_nu[id_nue] + inel.m1_mat_ab[id_anue][n+j][n+i] * (1. - g_nu[id_anue]);
      ann_term_ji[id_nux] = (pair.m1_mat_ab[id_nux][n+j][n+i] + brem.m1_mat_ab[0][n+j][n+i]) * g_nu[id_anux] + inel.m1_mat_ab[id_nux][n+j][n+i] * (1. - g_nu[id_nux]);
      ann_term_ji[id_anux] = (pair.m1_mat_ab[id_anux][n+j][n+i] + brem.m1_mat_ab[0][n+j][n+i]) * g_nu[id_nux] + inel.m1_mat_ab[id_anux][n+j][n+i] * (1. - g_nu[id_anux]);

      for (int idx = 0; idx < total_num_species; idx ++) {
        out_n->m1_mat_em[idx][n+i][n+j] = nu * nu * nu_bar * nu_bar * pro_term_ij[idx];
        out_n->m1_mat_em[idx][n+j][n+i] = nu_bar * nu_bar * nu * nu * pro_term_ji[idx];

        out_n->m1_mat_ab[idx][n+i][n+j] = nu * nu * g_nu[idx] * nu_bar * nu_bar * (pro_term_ij[idx] + ann_term_ij[idx]);
        out_n->m1_mat_ab[idx][n+j][n+i] = nu_bar * nu_bar * g_nu_bar[idx] * nu * nu * (pro_term_ji[idx] + ann_term_ji[idx]);

        out_j->m1_mat_em[idx][n+i][n+j] = nu * out_n->m1_mat_em[idx][n+i][n+j];
        out_j->m1_mat_em[idx][n+j][n+i] = nu_bar * out_n->m1_mat_em[idx][n+j][n+i];

        out_j->m1_mat_ab[idx][n+i][n+j] = nu * out_n->m1_mat_ab[idx][n+i][n+j];
        out_j->m1_mat_ab[idx][n+j][n+i] = nu_bar * out_n->m1_mat_ab[idx][n+j][n+i];
      }
    }

    for (int j = 0; j < n; j++) {

      // energies and parameters
      nu = t * quad_1d->points[i];
      nu_bar = t / quad_1d->points[j];

      for (int idx = 0; idx < total_num_species; idx++) {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
        g_nu_bar[idx] = TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
      }

      pro_term_ij[id_nue] = (pair.m1_mat_em[id_nue][i][n+j] + brem.m1_mat_em[0][i][n+j]) * (1. - g_nu_bar[id_anue]) + inel.m1_mat_em[id_nue][i][n+j] * g_nu_bar[id_nue];
      pro_term_ij[id_anue] = (pair.m1_mat_em[id_anue][i][n+j] + brem.m1_mat_em[0][i][n+j]) * (1. - g_nu_bar[id_nue]) + inel.m1_mat_em[id_anue][i][n+j] * g_nu_bar[id_anue];
      pro_term_ij[id_nux] = (pair.m1_mat_em[id_nux][i][n+j] + brem.m1_mat_em[0][i][n+j]) * (1. - g_nu_bar[id_anux]) + inel.m1_mat_em[id_nux][i][n+j] * g_nu_bar[id_nux];
      pro_term_ij[id_anux] = (pair.m1_mat_em[id_anux][i][n+j] + brem.m1_mat_em[0][i][n+j]) * (1. - g_nu_bar[id_nux]) + inel.m1_mat_em[id_anux][i][n+j] * g_nu_bar[id_anux];

      pro_term_ji[id_nue] = (pair.m1_mat_em[id_nue][n+j][i] + brem.m1_mat_em[0][n+j][i]) * (1. - g_nu[id_anue]) + inel.m1_mat_em[id_nue][n+j][i] * g_nu[id_nue];
      pro_term_ji[id_anue] = (pair.m1_mat_em[id_anue][n+j][i] + brem.m1_mat_em[0][n+j][i]) * (1. - g_nu[id_nue]) + inel.m1_mat_em[id_anue][n+j][i] * g_nu[id_anue];
      pro_term_ji[id_nux] = (pair.m1_mat_em[id_nux][n+j][i] + brem.m1_mat_em[0][n+j][i]) * (1. - g_nu[id_anux]) + inel.m1_mat_em[id_nux][n+j][i] * g_nu[id_nux];
      pro_term_ji[id_anux] = (pair.m1_mat_em[id_anux][n+j][i] + brem.m1_mat_em[0][n+j][i]) * (1. - g_nu[id_nux]) + inel.m1_mat_em[id_anux][n+j][i] * g_nu[id_anux];

      ann_term_ij[id_nue] = (pair.m1_mat_ab[id_nue][i][n+j] + brem.m1_mat_ab[0][i][n+j]) * g_nu_bar[id_anue] + inel.m1_mat_ab[id_nue][i][n+j] * (1. - g_nu_bar[id_nue]);
      ann_term_ij[id_anue] = (pair.m1_mat_ab[id_anue][i][n+j] + brem.m1_mat_ab[0][i][n+j]) * g_nu_bar[id_nue] + inel.m1_mat_ab[id_anue][i][n+j] * (1. - g_nu_bar[id_anue]);
      ann_term_ij[id_nux] = (pair.m1_mat_ab[id_nux][i][n+j] + brem.m1_mat_ab[0][i][n+j]) * g_nu_bar[id_anux] + inel.m1_mat_ab[id_nux][i][n+j] * (1. - g_nu_bar[id_nux]);
      ann_term_ij[id_anux] = (pair.m1_mat_ab[id_anux][i][n+j] + brem.m1_mat_ab[0][i][n+j]) * g_nu_bar[id_nux] + inel.m1_mat_ab[id_anux][i][n+j] * (1. - g_nu_bar[id_anux]);

      ann_term_ji[id_nue] = (pair.m1_mat_ab[id_nue][n+j][i] + brem.m1_mat_ab[0][n+j][i]) * g_nu[id_anue] + inel.m1_mat_ab[id_nue][n+j][i] * (1. - g_nu[id_nue]);
      ann_term_ji[id_anue] = (pair.m1_mat_ab[id_anue][n+j][i] + brem.m1_mat_ab[0][n+j][i]) * g_nu[id_nue] + inel.m1_mat_ab[id_anue][n+j][i] * (1. - g_nu[id_anue]);
      ann_term_ji[id_nux] = (pair.m1_mat_ab[id_nux][n+j][i] + brem.m1_mat_ab[0][n+j][i]) * g_nu[id_anux] + inel.m1_mat_ab[id_nux][n+j][i] * (1. - g_nu[id_nux]);
      ann_term_ji[id_anux] = (pair.m1_mat_ab[id_anux][n+j][i] + brem.m1_mat_ab[0][n+j][i]) * g_nu[id_nux] + inel.m1_mat_ab[id_anux][n+j][i] * (1. - g_nu[id_anux]);

      for (int idx = 0; idx < total_num_species; idx ++) {
        out_n->m1_mat_em[idx][i][n+j] = nu * nu * nu_bar * nu_bar * pro_term_ij[idx];
        out_n->m1_mat_em[idx][n+j][i] = nu_bar * nu_bar * nu * nu * pro_term_ji[idx];

        out_n->m1_mat_ab[idx][i][n+j] = nu * nu * g_nu[idx] * nu_bar * nu_bar * (pro_term_ij[idx] + ann_term_ij[idx]);
        out_n->m1_mat_ab[idx][n+j][i] = nu_bar * nu_bar * g_nu_bar[idx] * nu * nu * (pro_term_ji[idx] + ann_term_ji[idx]);

        out_j->m1_mat_em[idx][i][n+j] = nu * out_n->m1_mat_em[idx][i][n+j];
        out_j->m1_mat_em[idx][n+j][i] = nu_bar * out_n->m1_mat_em[idx][n+j][i];

        out_j->m1_mat_ab[idx][i][n+j] = nu * out_n->m1_mat_ab[idx][i][n+j];
        out_j->m1_mat_ab[idx][n+j][i] = nu_bar * out_n->m1_mat_ab[idx][n+j][i];
      }
    }
  }

  for (int idx = 0; idx < total_num_species; idx++) {
    for (int i = 0; i < 2 * n; i++) {
      free(pair.m1_mat_ab[idx][i]);
      free(pair.m1_mat_em[idx][i]);

      free(inel.m1_mat_ab[idx][i]);
      free(inel.m1_mat_em[idx][i]);
    }

    free(beta.m1_mat_ab[idx][0]);
    free(beta.m1_mat_em[idx][0]);

    free(pair.m1_mat_ab[idx]);
    free(pair.m1_mat_em[idx]);

    free(inel.m1_mat_ab[idx]);
    free(inel.m1_mat_em[idx]);

    free(beta.m1_mat_ab[idx]);
    free(beta.m1_mat_em[idx]);
  }

  for (int i = 0; i < 2 * n; i++) {
    free(brem.m1_mat_ab[0][i]);
    free(brem.m1_mat_em[0][i]);
  }

  free(brem.m1_mat_ab[0]);
  free(brem.m1_mat_em[0]);

  return;
}


/* Computes the opacities for the M1 code
 *
 */
M1Opacities ComputeM1Opacities(MyQuadrature *quad_1d, MyQuadrature *quad_2d, GreyOpacityParams *my_grey_opacity_params) {
  // compute some constants
  static const double four_pi_hc3 = (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
  static const double four_pi_hc3_sqr = four_pi_hc3 * four_pi_hc3;                 // [MeV^-6 cm^-6]

  // set up 1d integration
  MyFunctionMultiD integrand_m1_1d;
  MyQuadratureIntegrand integrand_m1_1d_info = {.n = 12};
  integrand_m1_1d.function = &M1CoeffsSingleIntegrand;
  integrand_m1_1d.dim = 1;
  integrand_m1_1d.params = my_grey_opacity_params;
  integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

  // set up 2d integration
  MyFunctionMultiD integrand_m1_2d;
  MyQuadratureIntegrand integrand_m1_2d_info = {.n = 16};
  integrand_m1_2d.function = &M1CoeffsDoubleIntegrand;
  integrand_m1_2d.dim = 2;
  integrand_m1_2d.params = my_grey_opacity_params;
  integrand_m1_2d.my_quadrature_integrand = integrand_m1_2d_info;

  double n[total_num_species];
  double J[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    n[idx] = my_grey_opacity_params->m1_pars.n[idx];
    J[idx] = my_grey_opacity_params->m1_pars.J[idx];
  }

  // double s_2d_x[12];
  // double s_2d_y[12];

  // for (int i=0; i<11; i++) {
  // s_2d_x[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
  // s_2d_y[i] = s_2d_x[i];
  //}
  // @TODO: choose this appropriately

  // double s_1d[11];
  // for (int i=0; i<11; i++) {
  // s_1d[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
  //}

  const double s = 1.5 * my_grey_opacity_params->eos_pars.temp;
  MyQuadratureIntegrand integrals_1d = GaussLegendreIntegrateFixedSplit1D(quad_1d, &integrand_m1_1d, s);

  //MyQuadratureIntegrand integrals_2d = GaussLegendreIntegrateFixedSplit2D(quad_2d, &integrand_m1_2d, s);
  
  MyQuadratureIntegrand integrals_2d = GaussLegendreIntegrateTest2D(quad_1d, my_grey_opacity_params, s);

  
  // printf("%.5e\n", integrals_2d.integrand[0]);

  M1Opacities m1_opacities;

  m1_opacities.eta_0[id_nue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[0] + integrals_1d.integrand[0]);
  m1_opacities.eta_0[id_anue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[1] + integrals_1d.integrand[1]);
  m1_opacities.eta_0[id_nux] = four_pi_hc3_sqr * integrals_2d.integrand[2];
  m1_opacities.eta_0[id_anux] = four_pi_hc3_sqr * integrals_2d.integrand[3];

  m1_opacities.kappa_0_a[id_nue] = four_pi_hc3 / (kClight * n[id_nue]) * (four_pi_hc3 * integrals_2d.integrand[4] + integrals_1d.integrand[2]);
  m1_opacities.kappa_0_a[id_anue] = four_pi_hc3 / (kClight * n[id_anue]) * (four_pi_hc3 * integrals_2d.integrand[5] + integrals_1d.integrand[3]);
  m1_opacities.kappa_0_a[id_nux] = four_pi_hc3_sqr / (kClight * n[id_nux]) * integrals_2d.integrand[6];
  m1_opacities.kappa_0_a[id_anux] = four_pi_hc3_sqr / (kClight * n[id_anux]) * integrals_2d.integrand[7];

  m1_opacities.eta[id_nue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[8] + integrals_1d.integrand[4]);
  m1_opacities.eta[id_anue] = four_pi_hc3 * (four_pi_hc3 * integrals_2d.integrand[9] + integrals_1d.integrand[5]);
  m1_opacities.eta[id_nux] = four_pi_hc3_sqr * integrals_2d.integrand[10];
  m1_opacities.eta[id_anux] = four_pi_hc3_sqr * integrals_2d.integrand[11];

  m1_opacities.kappa_a[id_nue] = four_pi_hc3 / (kClight * J[id_nue]) * (four_pi_hc3 * integrals_2d.integrand[12] + integrals_1d.integrand[6]);
  m1_opacities.kappa_a[id_anue] = four_pi_hc3 / (kClight * J[id_anue]) * (four_pi_hc3 * integrals_2d.integrand[13] + integrals_1d.integrand[7]);
  m1_opacities.kappa_a[id_nux] = four_pi_hc3_sqr / (kClight * J[id_nux]) * integrals_2d.integrand[14];
  m1_opacities.kappa_a[id_anux] = four_pi_hc3_sqr / (kClight * J[id_anux]) * integrals_2d.integrand[15];

  m1_opacities.kappa_s[id_nue] = four_pi_hc3 / (kClight * J[id_nue]) * integrals_1d.integrand[8];
  m1_opacities.kappa_s[id_anue] = four_pi_hc3 / (kClight * J[id_anue]) * integrals_1d.integrand[9];
  m1_opacities.kappa_s[id_nux] = four_pi_hc3 / (kClight * J[id_nux]) * integrals_1d.integrand[10];
  m1_opacities.kappa_s[id_anux] = four_pi_hc3 / (kClight * J[id_anux]) * integrals_1d.integrand[11];

  return m1_opacities;
}

/* Compute the integrands for the computation of the spectral emissivity and inverse mean free path */
MyQuadratureIntegrand SpectralIntegrand(double *var, void *p) {
  // energies and parameters
  double nu_bar = var[0]; // [MeV]

  GreyOpacityParams *my_grey_opacity_params = (GreyOpacityParams *)p;
  MyEOSParams my_eos_params = my_grey_opacity_params->eos_pars;
  OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
  }

  // compute the pair kernels
  MyKernelQuantity pair_kernels_m1 = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  if (opacity_flags.use_pair) {
    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega_prime = nu_bar;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = 1.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.filter = 0.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax = 0;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.mu = 1.;
    my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime = 1.;
    // pair_kernels_m1 = PairKernelsM1(&my_eos_params, &my_grey_opacity_params->kernel_pars.pair_kernel_params);
    pair_kernels_m1 = PairKernelsM1Optimized(&my_eos_params, &my_grey_opacity_params->kernel_pars.pair_kernel_params);
  }

  // compute the bremsstrahlung kernels
  MyKernelOutput brem_kernels_m1 = {0};
  if (opacity_flags.use_brem) {
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega_prime = nu_bar;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
    brem_kernels_m1 = BremKernelsLegCoeff(&my_grey_opacity_params->kernel_pars.brem_kernel_params, &my_eos_params);
  }

  // compute the inelastic NES/NPS kernels
  MyKernelOutput inelastic_kernels_m1 = {0};
  if (opacity_flags.use_inelastic_scatt) {
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
    inelastic_kernels_m1 = InelasticScattKernels(&my_grey_opacity_params->kernel_pars.inelastic_kernel_params, &my_grey_opacity_params->eos_pars);
  }

  double pro_term[total_num_species];

  pro_term[id_nue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_nue]) * (1. - g_nu[id_anue]) + inelastic_kernels_m1.em[id_nue] * g_nu[id_nue];
  pro_term[id_anue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_anue]) * (1. - g_nu[id_nue]) + inelastic_kernels_m1.em[id_anue] * g_nu[id_anue];
  pro_term[id_nux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_nux]) * (1. - g_nu[id_anux]) + inelastic_kernels_m1.em[id_nux] * g_nu[id_nux];
  pro_term[id_anux] = (pair_kernels_m1.em_x + brem_kernels_m1.em[id_anux]) * (1. - g_nu[id_nux]) + inelastic_kernels_m1.em[id_anux] * g_nu[id_anux];

  double ann_term[total_num_species];

  ann_term[id_nue] = (pair_kernels_m1.abs_e + brem_kernels_m1.abs[id_nue]) * g_nu[id_anue] + inelastic_kernels_m1.abs[id_nue] * (1. - g_nu[id_nue]);
  ann_term[id_anue] = (pair_kernels_m1.abs_e + brem_kernels_m1.abs[id_anue]) * g_nu[id_nue] + inelastic_kernels_m1.abs[id_anue] * (1. - g_nu[id_anue]);
  ann_term[id_nux] = (pair_kernels_m1.abs_x + brem_kernels_m1.abs[id_nux]) * g_nu[id_anux] + inelastic_kernels_m1.abs[id_nux] * (1. - g_nu[id_nux]);
  ann_term[id_anux] = (pair_kernels_m1.abs_x + brem_kernels_m1.abs[id_anux]) * g_nu[id_nux] + inelastic_kernels_m1.abs[id_anux] * (1. - g_nu[id_anux]);

  double integrand_1[total_num_species], integrand_2[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    integrand_1[idx] = nu_bar * nu_bar * pro_term[idx];
    integrand_2[idx] = nu_bar * nu_bar * ann_term[idx];
  }

  MyQuadratureIntegrand result = {.n = 8};

  result.integrand[0] = integrand_1[id_nue];
  result.integrand[1] = integrand_1[id_anue];
  result.integrand[2] = integrand_1[id_nux];
  result.integrand[3] = integrand_1[id_anux];
  result.integrand[4] = integrand_2[id_nue];
  result.integrand[5] = integrand_2[id_anue];
  result.integrand[6] = integrand_2[id_nux];
  result.integrand[7] = integrand_2[id_anux];

  return result;
}

/* Computes the spectral emissivity and inverse mean free path */
SpectralOpacities ComputeSpectralOpacitiesNotStimulated(const double nu, MyQuadrature *quad_1d, GreyOpacityParams *my_grey_opacity_params) {
  // compute some constants
  static const double four_pi_hc3 = (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]

  // set up 1d integration
  MyFunctionMultiD integrand_m1_1d;
  MyQuadratureIntegrand integrand_m1_1d_info = {.n = 8};
  integrand_m1_1d.function = &SpectralIntegrand;
  integrand_m1_1d.dim = 1;
  integrand_m1_1d.params = my_grey_opacity_params;
  integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

  // compute the neutrino & anti-neutrino distribution function
  double g_nu[total_num_species];

  for (int idx = 0; idx < total_num_species; idx++) {
    g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
  }

  double s[8];
  for (int i = 0; i < 8; i++) {
    s[i] = nu;
  }

  my_grey_opacity_params->kernel_pars.pair_kernel_params.omega = nu;
  my_grey_opacity_params->kernel_pars.brem_kernel_params.omega = nu;
  my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;

  MyQuadratureIntegrand integrals_1d = GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s);

  MyOpacity abs_em_beta = {0};
  if (my_grey_opacity_params->opacity_flags.use_abs_em) {
    abs_em_beta = StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars); // [s^-1]
  }

  double iso_scatt = 0.;
  if (my_grey_opacity_params->opacity_flags.use_iso) {
    iso_scatt = IsoScattLegCoeff(nu, &my_grey_opacity_params->opacity_pars, &my_grey_opacity_params->eos_pars, 0);
  }

  SpectralOpacities sp_opacities;

  sp_opacities.j[id_nue] = abs_em_beta.em[id_nue] + four_pi_hc3 * integrals_1d.integrand[0];
  sp_opacities.j[id_anue] = abs_em_beta.em[id_anue] + four_pi_hc3 * integrals_1d.integrand[1];
  sp_opacities.j[id_nux] = abs_em_beta.em[id_nux] + four_pi_hc3 * integrals_1d.integrand[2];
  sp_opacities.j[id_anux] = abs_em_beta.em[id_anux] + four_pi_hc3 * integrals_1d.integrand[3];

  sp_opacities.kappa[id_nue] = (abs_em_beta.abs[id_nue] + four_pi_hc3 * integrals_1d.integrand[4]) / kClight;
  sp_opacities.kappa[id_anue] = (abs_em_beta.abs[id_anue] + four_pi_hc3 * integrals_1d.integrand[5]) / kClight;
  sp_opacities.kappa[id_nux] = (abs_em_beta.abs[id_nux] + four_pi_hc3 * integrals_1d.integrand[6]) / kClight;
  sp_opacities.kappa[id_anux] = (abs_em_beta.abs[id_anux] + four_pi_hc3 * integrals_1d.integrand[7]) / kClight;

  sp_opacities.j_s[id_nue] = 4. * kPi * nu * nu * g_nu[id_nue] * iso_scatt;
  sp_opacities.j_s[id_anue] = 4. * kPi * nu * nu * g_nu[id_anue] * iso_scatt;
  sp_opacities.j_s[id_nux] = 4. * kPi * nu * nu * g_nu[id_nux] * iso_scatt;
  sp_opacities.j_s[id_anux] = 4. * kPi * nu * nu * g_nu[id_anux] * iso_scatt;

  sp_opacities.kappa_s[id_nue] = 4. * kPi * nu * nu * iso_scatt / kClight;
  sp_opacities.kappa_s[id_anue] = sp_opacities.kappa_s[id_nue];
  sp_opacities.kappa_s[id_nux] = sp_opacities.kappa_s[id_nue];
  sp_opacities.kappa_s[id_anux] = sp_opacities.kappa_s[id_nue];

  return sp_opacities;
}
