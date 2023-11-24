//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_m1_coeffs.c
//  \brief test routines for the grey m1 source coefficients

#include <stdio.h>
#include <math.h>

#include "../include/constants.h"
#include "../include/bns_nurates.h"
#include "opacities.h"
#include "../include/distribution.h"

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


  //double n = NuNumber(&distr_pars);
  //double J = NuEnergy(&distr_pars);

  //M1Quantities m1_pars = {.n = n, .J = J, .chi = 1./3.};

  //GreyOpacityParams grey_pars = {.opacity_pars = opacity_pars,
  //                               .eos_pars = EOS_pars,
  //                               .distr_pars = distr_pars,
  //                               .m1_pars = m1_pars};


  return 0;
}