//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
#define BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_

#include "constants.h"

/*===========================================================================*/

// nucfrmfac.h

// Computation of single nucleon form factors
void NucFrmFac(const double E, double* cv, double* ca, double* F2,
               const int reacflag);

/*===========================================================================*/

// weak_magnetism.c

// Correction for neutrino absorption on neutron (nue + n -> e- + p) +
// Correction for electron antineutrino absorption on proton (anue + p -> e+ +
// n)
void WMAbsEm(const double omega, double* R, double* Rbar);

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N)
void WMScatt(const double omega, double* R0, double* R1, const int reacflag);

/*===========================================================================*/

#endif // BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
