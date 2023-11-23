// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_opt_thin.c
//  \brief compute neutrino distribution function & M1 parameters for optically thin regime

#include <math.h>

#include "distribution.h"
#include "constants.h"
#include "functions.h"

/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann distribution
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 * species:     species of neutrino
 */
double NuFThin(double omega, NuDistributionParams *distr_pars, int species) {
  double temp_f = distr_pars->temp_f[species];
  double c_f = distr_pars->c_f[species];
  return pow(omega, c_f) * exp(-omega / temp_f);
}

/* Recover distribution function parameters for optically thin regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * out_distr_pars:  computes free neutrino temperature and c_f
 */
void CalculateThinParamsFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars) {

  static const double exp3 = 20.085536923187668;
  static const double e = 2.718281828459045235360287471352;
  static const double ie = 1. / e;
  static const double isqrt_32_pi3 = 0.031746817967120484;

  for (int species = 0; species < total_num_species; species++) {
    double n = M1_pars->n[species];
    double j = M1_pars->J[species];

    double j_over_n = j / n;
    double A = j_over_n / e;
    double B = n * isqrt_32_pi3 * (kHClight * kHClight * kHClight); //@TODO: warning!
    double log_A = log(A);

    double argW0 = -2. * log_A / (B * B);
    double x = argW0 > -ie ? -W0(argW0) / (2. * log_A) : 3.;

    out_distr_pars->w_f[species] = 0.5 * (3. * M1_pars->chi[species] - 1.);

    out_distr_pars->c_f[species] = fmax(x - 3., 0.);
    out_distr_pars->temp_f[species] = j_over_n / x;

  }
}