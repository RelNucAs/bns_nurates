// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file opt_thin.h
//  \brief compute neutrino distribution function for optically thin regime

#include <math.h>

#include "distribution.h"
#include "../constants.h"
#include "../functions/functions.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif

#define kExp exp(1.)
#define kSqrt32PiThree sqrt(32. * kPi * kPi * kPi)

/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann distribtion
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 */
double NuFThin(double omega, NuDistributionParams *distr_pars) {
  double temp_f = distr_pars->temp_f;
  double c_f = distr_pars->c_f;
  return pow(omega, c_f) * exp(-omega / temp_f);
}

/* Recover distribution function parameters for optically thin regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * out_distr_pars:  computes free neutrino temperature and c_f
 */
void CalculateThinParamsFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars) {
  double n = M1_pars->n;
  double j = M1_pars->J;

  double j_over_n = j / n;
  double A = j_over_n / kExp;
  double logA = log(A);
  double B = n / kSqrt32PiThree;

  double lambert = W0(-2 * logA / (B * B));

  if (lambert < 0.) {
    DEBUG_PRINT(("WARNING -> opt_thin.c: negative value for Lambert function\n"));
  }

  lambert = (lambert > 0.) ? lambert : 0.;
  double x = -0.5 * lambert / logA;

  out_distr_pars->w_f = 1.5 * (1. - M1_pars->chi);
  out_distr_pars->c_f = x - 3.;
  out_distr_pars->temp_f = j_over_n / x;

}