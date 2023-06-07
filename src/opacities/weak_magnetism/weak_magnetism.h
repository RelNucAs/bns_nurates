//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
#define BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_

#include "../../constants.h"

/*===========================================================================*/

// nucfrmfac.h

// Computation of single nucleon form factors
void nucfrmfac(const double E,
	       double* cv, double* ca, double* F2,
	       const int reacflag);

/*===========================================================================*/

// weak_magnetism.c

// Correction for electron neutrino absorption on neutron (nue + n -> e- + p)
double WM_nue_abs(const double e_nu);

// Correction for electron antineutrino absorption on proton (anue + p -> e+ + n)
double WM_anue_abs(const double e_nu);

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N)
void WM_scatt(const double e_nu, double* R0, double* R1, const int reacflag);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
