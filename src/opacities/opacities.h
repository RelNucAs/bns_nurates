//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include <stdbool.h>

#include "../bns_nurates.h"


/*===========================================================================*/

// nu_abs_em_beta.c

// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Emission coefficients
SourceCoeffs NuclAbsEmissionCoeffs(MyEOSParams *eos_pars);

// Absoprtion coefficients
SourceCoeffs NuclAbsOpacityCoeffs(GreyOpacityParams *grey_pars);

/*===========================================================================*/

// nu_scatt_iso.c

// Isoenergetic scattering on protons
double IsoScattProton(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Isoenergetic scattering on neutrons
double IsoScattNeutron(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Sum of the two contributions for neutrino-nucleon scattering (protons + neutrons)
double IsoScattTotal(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Scattering coefficients
SourceCoeffs IsoScattCoeffs(GreyOpacityParams *grey_pars);

/*===========================================================================*/

// nu_grey_total.c

// Total emission coefficients
SourceCoeffs EmissionCoeffs(MyEOSParams *eos_pars);

// Total opacity coefficients
SourceCoeffs OpacityCoeffs(GreyOpacityParams *grey_pars);

// Total scattering coefficients
SourceCoeffs ScatteringCoeffs(GreyOpacityParams *grey_pars);

/*===========================================================================*/

// pair opacities
MyOpacityQuantity PairEmissivityAbsorptivityIntegrand(double var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
MyOpacityQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
// bremsstrahlung opacities
MyOpacityQuantity BremEmissivityAbsorptivityIntegrandFermi(double omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
MyOpacityQuantity BremOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);

// M1 opacities
M1Opacities M1OpacitiesIntegrand(var3d x, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params, GreyOpacityParams *my_grey_opacity_params);
M1Opacities ComputeM1Opacities(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params, GreyOpacityParams *my_grey_opacity_params);

#endif //BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
