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
  double beta_f = distr_pars->beta_f[species];

  return pow(omega / temp_f, c_f) * exp(-omega / temp_f); // @TODO: we need to change this later
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

    out_distr_pars->w_f[species] = 0.5 * (3. * M1_pars->chi[species] - 1.);

    const double MYCF = 0.6;
    out_distr_pars->c_f[species] = MYCF;
    out_distr_pars->temp_f[species] = j / (n * (MYCF + 3.));;
    out_distr_pars->beta_f[species] = n * (kHClight * kHClight * kHClight) / (4. * M_PI * GammaStirling(MYCF + 3.) * pow(out_distr_pars->temp_f[species], MYCF + 3.));

    // alternate reconstruction
    //if(j/n - 3. < 0){
    //  printf("Warning! c_f is going to be negative! c_f: %lf\n", j/n - 3.);
    //}
    //out_distr_pars->c_f[species] = fmax(j/n - 3., 0.);
    //out_distr_pars->temp_f[species] = n * (kHClight * kHClight * kHClight)  / (4. * M_PI * GammaStirling(out_distr_pars->c_f[species] + 3.));
    //out_distr_pars->beta_f[species] = -42.;

  }
}