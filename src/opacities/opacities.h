//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include <stdbool.h>

#include "../bns_nurates.h"

// Opacity parameters
struct OpacityParams {
  double omega;    // (anti)neutrino energy [MeV]
  bool use_dU;     // flag for dU correction
  bool use_WM_ab;  // flag for WM correction (and related) on absorption rates
  bool use_WM_sc;  // flag for WM correction (and related) on scattering rates
};
typedef struct OpacityParams OpacityParams;

/*===========================================================================*/

// nu_abs_em.h

// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(OpacityParams *opacity_pars, MyEOSParams *eos_pars);

/*===========================================================================*/

// pair opacities
MyOpacityQuantity PairEmissivityAbsorptivityIntegrandFermi(double omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
MyOpacityQuantity PairOpacitiesFermi(MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
#endif //BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
