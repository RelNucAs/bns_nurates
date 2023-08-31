//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_m1_coeffs.c
//  \brief test routines for the grey m1 source coefficients

#include <stdio.h>
#include <math.h>

#include "../src/constants.h"
#include "../src/bns_nurates.h"
#include "../src/opacities/opacities.h"
#include "../src/distribution/distribution.h"

int main() {

  // Define thermodynamics quantities
  double nb = 3.730000e+14 / kMb;
  double temp = 1.204000e+01;
  double lep_mass = kMe;
  double yp = 3.134000e-01;
  double yn = 6.866000e-01;
  double mu_l = 2.495089e+02;
  double mu_hat = 6.015553e+01;
  double omega = 10.08;

  // EOS parameters
  MyEOSParams EOS_pars = {.nb = nb,
                          .temp = temp,
                          .yp = yp,
                          .yn = yn,
                          .mu_e = mu_l,
                          .mu_mu = 0.,
                          .mu_p = 0.,
                          .mu_n = mu_hat + kQ,
                          .dU = 0.};

  // Opacity paraneters
  OpacityParams opacity_pars = {.use_dU = 0, .use_WM_ab = 0, .use_WM_sc = 0};

  // Neutrino distribution function parameters (thick + thin)
  NuDistributionParams distr_pars = {.w_t = 1., .temp_t = 10., .eta_t = 4.,
                                     .w_f = 0., .temp_f = 10., .c_f = 0.};


  double n = NuNumber(&distr_pars);
  double J = NuEnergy(&distr_pars);

  M1Quantities m1_pars = {.n = n, .J = J, .chi = 1./3.};

  GreyOpacityParams grey_pars = {.opacity_pars = opacity_pars,
                                 .eos_pars = EOS_pars,
                                 .distr_pars = distr_pars,
                                 .m1_pars = m1_pars};

  /* Compute emission coefficients */

  // SourceCoeffs em_coeffs = ... 

  // Print results to screen
  //printf("Emission coefficients:\n");
  //printf("R_nue = %.5e, R_anue = %.5e, R_nux = %.5e\n", em_coeffs.R_nue, em_coeffs.R_anue , em_coeffs.R_nux);
  //printf("Q_nue = %.5e, Q_anue = %.5e, Q_nux = %.5e\n", em_coeffs.Q_nue, em_coeffs.Q_anue , em_coeffs.Q_nux);
  //printf("\n");


  /* Compute opacity coefficients */
  
  // SourceCoeffs ab_coeffs = ...

  // Print result to screen
  //printf("Opacity coefficients:\n");
  //printf("R_nue = %.5e, R_anue = %.5e, R_nux = %.5e\n", ab_coeffs.R_nue, ab_coeffs.R_anue , ab_coeffs.R_nux);
  //printf("Q_nue = %.5e, Q_anue = %.5e, Q_nux = %.5e\n", ab_coeffs.Q_nue, ab_coeffs.Q_anue , ab_coeffs.Q_nux);
  //printf("\n");


  /* Compute scattering coefficients */
  
  // SourceCoeffs sc_coeffs = ...

  // Print result to screen
  //printf("Scattering coefficients:\n");
  //printf("R_nue = %.5e, R_anue = %.5e, R_nux = %.5e\n", sc_coeffs.R_nue, sc_coeffs.R_anue , sc_coeffs.R_nux);
  //printf("Q_nue = %.5e, Q_anue = %.5e, Q_nux = %.5e\n", sc_coeffs.Q_nue, sc_coeffs.Q_anue , sc_coeffs.Q_nux);
  //printf("\n");

  return 0;
}