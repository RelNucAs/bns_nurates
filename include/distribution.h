// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution.h
//  \brief header file for distribution function reconstruction from M1 parameters
//         supports three different species: electron neutrino, electron anti-neutrino and mu/tau neutrino/antineutrino

#ifndef BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
#define BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_

#include "bns_nurates.h"

/* ===========================================================================
 * Functions for the optically thick distribution function
 * ===========================================================================
 */

// Neutrino distribution function for the optically thick regime
double NuFThick(double omega, NuDistributionParams *distr_pars, int species);

// Recover parameters of thick distribution function from M1 quantities
void CalculateThickParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, NuDistributionParams *out_distr_pars);

// Recover parameters considering equilibrium with fluid
void CalculateEquilibriumParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, NuDistributionParams *out_distr_pars);
/* ===========================================================================
 * Functions for the optically thin distribution function
 * ===========================================================================
 */

// Neutrino distribution function in optically thin regime
double NuFThin(double omega, NuDistributionParams *distr_pars, int species);

// Recover parameters of thin distribution function from M1 quantities
void CalculateThinParamsFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars);

/* ===========================================================================
 * Functions for constructing the total distribution function
 * ===========================================================================
 */

// Total neutrino distribution combining optically thick and thin regimes
double TotalNuF(double omega, NuDistributionParams *distr_pars, int species);

// Recover parameters of thick and thin distribution function from M1 quantities
NuDistributionParams CalculateDistrParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, bool use_equilibrium);

/* ===========================================================================
 * Function for evaluating parameters of neutrino distribution function at equilibrium
 * ===========================================================================
 */

NuDistributionParams NuEquilibriumParams(MyEOSParams *eos_pars);

/* ===========================================================================
 * Functions for computing neutrino number and energy density
 * ===========================================================================
 */

// Neutrino number density
MyQuadratureIntegrand NuNumber(NuDistributionParams *distr_pars);

// Neutrino energy density
MyQuadratureIntegrand NuEnergy(NuDistributionParams *distr_pars);

#endif //BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
