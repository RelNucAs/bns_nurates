// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_opt_thick.h
//  \brief compute neutrino distribution function & M1 parameters for optically thick regime

#include <math.h>

#include "distribution.h"
#include "bns_nurates.h"
#include "functions.h"

#define kFourThirds 1.333333333333333333333333
#define kFitA  1.3300219438752203
#define kFitB  368.1278987362366
#define kFitC (-2701.5865312193273)
#define kFitD 3930.881704339285
#define kFitE 11356.551670239718
#define kFitF (-33734.77092570755)
#define kFitG 21233.158792661736
#define kFitH 275.93123409700814
#define kFitI (-2015.9177416521718)
#define kFitL 4403.086279325216
#define kFitM (-1871.7098365259042)
#define kFitN (-2166.313342616665)

/* Neutrino distribution function in optically thick regime: Fermi-Dirac distribution
 *
 * omega:       neutrino energy
 * distr_pars:  uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 * species:     species of neutrino
 */
double NuFThick(double omega, NuDistributionParams *distr_pars, int species) {
  double temp = distr_pars->temp_t[species];
  double mu = distr_pars->temp_t[species] * distr_pars->eta_t[species]; // @TODO: simplify if needed
  return FermiDistr(omega, temp, mu);
}

/* Recover distribution function parameters for optically thick regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * eos_pars:        uses fluid temperature
 * out_distr_pars:  computes trapped neutrino temperature and degeneracy parameter
 */
void CalculateThickParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, NuDistributionParams *out_distr_pars) {

  // set degeneracy parameter for different neutrino species
  out_distr_pars->eta_t[0] = (eos_pars->mu_p + eos_pars->mu_e - eos_pars->mu_n) / eos_pars->temp;
  out_distr_pars->eta_t[1] = -(eos_pars->mu_p + eos_pars->mu_e - eos_pars->mu_n) / eos_pars->temp;
  out_distr_pars->eta_t[2] = 0.;

  for (int species = 0; species < total_num_species; species++) {

    out_distr_pars->w_t[species] = 1.5 * (1. - M1_pars->chi[species]);
    out_distr_pars->temp_t[species] = eos_pars->temp;

    // @TODO: Currently disabling alternative. What do we do with this ?
    double y = fmax(M1_pars->J[species] / (M1_pars->n[species] * eos_pars->temp), 3.05); // average neutrino energy over thermal energy
    double y_0 = 114.;

    out_distr_pars->eta_t[species] =
        y < y_0 ? (y * (y * (y * (y * (y * (kFitA * y + kFitB) + kFitC) + kFitD) + kFitE) + kFitF) + kFitG) / (y * (y * (y * (y * (y + kFitH) + kFitI) + kFitL) + kFitM) + kFitN) : kFourThirds * y;

  }
}