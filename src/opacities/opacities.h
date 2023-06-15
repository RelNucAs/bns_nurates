//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include "../bns_nurates.h"

/*===========================================================================*/

// nu_abs_em.h

// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(const double omega, MyEOSParams *eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(const double omega, MyEOSParams *eos_pars);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
