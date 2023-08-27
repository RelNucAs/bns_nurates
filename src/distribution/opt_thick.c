// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file opt_thick.h
//  \brief compute neutrino distribution function for optically thick regime

#include <math.h>

#include "distribution.h"
#include "../bns_nurates.h"
#include "../functions/functions.h"

/* Neutrino distribution function in optically thick regime: Fermi-Dirac distribution
 *
 * omega: neutrino energy
 * distr_pars: uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 */
double NuFThick(double omega, NuDistributionParams *distr_pars) {
  double temp = distr_pars->temp_t;
  double mu = temp * distr_pars->eta_t;
  return FermiDistr(omega, temp, mu);
}

/* Recover distribution function parameters for optically thick regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * eos_pars:        uses fluid temperature
 * out_distr_pars:  computes trapped neutrino temperature and degeneracy parameter
 */
void CalculateThickParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, NuDistributionParams *out_distr_pars) {

  // trapped neutrinos taken to be in equilibrium with fluid
  out_distr_pars->eta_t = eos_pars->temp;

  // calculate trapped distribution contribution from Eddington factor
  out_distr_pars->w_t = 0.5 * (3. * M1_pars->chi - 1.);

  // calculate degeneracy parameter
  double y = M1_pars->J / (M1_pars->n * eos_pars->temp); // average neutrino energy over thermal energy
  const double y_0 = 114.;

  // least square fit parameters from Federico's notes
  const double kA = 1.3300219438752203;
  const double kB = 368.1278987362366;
  const double kC = -2701.5865312193273;
  const double kD = 3930.881704339285;
  const double kE = 11356.551670239718;
  const double kF = -33734.77092570755;
  const double kG = 21233.158792661736;
  const double kH = 275.93123409700814;
  const double kI = -2015.9177416521718;
  const double kL = 4403.086279325216;
  const double kM = -1871.7098365259042;
  const double kN = -2166.313342616665;

  // limit y to values greater than 3.05
  y = (y > 3.05) ? y : 3.05;

  out_distr_pars->eta_t = (y < y_0) ? (kA * pow(y, 6) + kB * pow(y, 5) + kC * pow(y, 4) + kD * pow(y, 3) + kE * pow(y, 2) + kF * y + kG)
      / (pow(y, 5) + kH * pow(y, 4) + kI * pow(y, 3) + kL * pow(y, 2) + kM * y + kN) : (4. / 3.) * y;

}