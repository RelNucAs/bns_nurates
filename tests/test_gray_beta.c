//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_brem.c
//  \brief test routines for the grey coefficients of the beta processes

#include <stdio.h>
#include <math.h>

#include "../src/constants.h"
#include "../src/bns_nurates.h"
#include "../src/opacities/opacities.h"

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

  MyEOSParams EOS_pars = {.nb = nb,
                          .temp = temp,
                          .yp = yp,
                          .yn = yn,
                          .mu_e = mu_l,
                          .mu_mu = 0.,
                          .mu_p = 0.,
                          .mu_n = mu_hat + kQ,
                          .dU = 0.};

  // Compute emission coefficient (in stimulated absorption framework)
  SourceCoeffs coeffs = NuclAbsEmissionCoeffs(&EOS_pars);

  // Print result to screen
  printf("Emission coefficients output:\n");
  printf("R_nue = %.5e, R_anue = %.5e, R_nux = %.5e\n", coeffs.R_nue, coeffs.R_anue , coeffs.R_nux);
  printf("Q_nue = %.5e, Q_anue = %.5e, Q_nux = %.5e\n", coeffs.Q_nue, coeffs.Q_anue , coeffs.Q_nux);
    
  return 0;
}