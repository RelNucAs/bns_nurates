// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file opacities.h
//  \brief header files for opacity functions

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include <stdbool.h>

#include "bns_nurates.h"


/*===========================================================================*/

// nu_abs_em_beta.c

// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

/*===========================================================================*/

// nu_scatt_iso.c

// Isoenergetic scattering on protons
double IsoScattProton(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Isoenergetic scattering on neutrons
double IsoScattNeutron(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Sum of the two contributions for neutrino-nucleon scattering (protons + neutrons)
double IsoScattTotal(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars);

// Legendre coefficient of order l fot the isoenergetic scattering kernel
double IsoScattLegCoeff(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const int l);

/*===========================================================================*/

// pair opacities
MyKernelQuantity PairEmissivityAbsorptivityIntegrand(double var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
MyKernelQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
// bremsstrahlung opacities
MyKernelQuantity BremEmissivityAbsorptivityIntegrandFermi(double omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);
MyKernelQuantity BremOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params);

// M1 opacities
MyQuadratureIntegrand M1CoeffsDoubleIntegrand(double *var, void *p);
MyQuadratureIntegrand M1CoeffsSingleIntegrand(double *var, void *p);
M1Opacities ComputeM1Opacities(MyQuadrature *quad_1d, MyQuadrature *quad_2d, GreyOpacityParams *my_grey_opacity_params);

SpectralOpacities ComputeSpectralOpacities(const double nu, MyQuadrature *quad_1d, GreyOpacityParams *my_grey_opacity_params);


#endif //BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
