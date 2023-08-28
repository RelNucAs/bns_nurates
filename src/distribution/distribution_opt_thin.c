// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_opt_thin.c
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

/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann distribtion
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 */
double NuFThin(double omega, NuDistributionParams *distr_pars) {
  double temp_f = distr_pars->temp_f;
  double c_f = distr_pars->c_f;
  return pow(omega, c_f) * exp(- omega / temp_f);
}

/* Recover distribution function parameters for optically thin regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * out_distr_pars:  computes free neutrino temperature and c_f
 */
void CalculateThinParamsFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars) {
  static const double exp3         = 20.085536923187668;
  static const double e            = 2.718281828459045235360287471352;
  static const double ie           = 1. / e;
  static const double isqrt_32_pi3 = 0.031746817967120484;
  
  const double n = M1_pars->n;
  const double j = M1_pars->J;

  const double j_over_n = j / n;
  const double A = j_over_n / e;
  const double B = n * isqrt_32_pi3;
  const double log_A = log(A);

  const double argW0 = -2. * log_A / (B * B);
  double x;
  x = argW0 > -ie ? -W0(argW0) / (2. * log_A) : 3.;

  //if (x < 0.) {
    //DEBUG_PRINT(("WARNING -> opt_thin.c: negative value for Lambert function\n"));
  //}

  // calculate free streaming distribution contribution from Eddington factor
  out_distr_pars->w_f = 1.5 * (1. - M1_pars->chi);

  out_distr_pars->c_f = fmax(x - 3., 0.);
  out_distr_pars->temp_f = j_over_n / x;

  return;
}