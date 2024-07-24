// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file opacities.h
//  \brief header files for opacity functions

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include <stdbool.h>

#include "bns_nurates.hpp"

#define kirchoff_flag false

/*===========================================================================*/

// nu_abs_em_beta.c

// (Anti)neutrino absorption on nucleons
MyOpacity AbsOpacity(const double omega, OpacityParams* opacity_pars,
                     MyEOSParams* eos_pars);

// Stimulated absoption version
MyOpacity StimAbsOpacity(const double omega, OpacityParams* opacity_pars,
                         MyEOSParams* eos_pars);

// Build matrix for integration
void BetaOpacitiesTable(MyQuadrature* quad, MyEOSParams* eos_pars,
                        OpacityParams* opacity_pars, double t, M1Matrix* out);

/*===========================================================================*/

// nu_scatt_iso.c

// Isoenergetic scattering on protons
double IsoScattProton(const double omega, OpacityParams* opacity_pars,
                      MyEOSParams* eos_pars);

// Isoenergetic scattering on neutrons
double IsoScattNeutron(const double omega, OpacityParams* opacity_pars,
                       MyEOSParams* eos_pars);

// Sum of the two contributions for neutrino-nucleon scattering (protons +
// neutrons)
double IsoScattTotal(const double omega, OpacityParams* opacity_pars,
                     MyEOSParams* eos_pars);

// Legendre coefficient of order l fot the isoenergetic scattering kernel
double IsoScattLegCoeff(const double omega, OpacityParams* opacity_pars,
                        MyEOSParams* eos_pars, const int l);

/*===========================================================================*/

// pair opacities
#ifdef GSL_INCLUDES_H_
MyKernelQuantity
PairEmissivityAbsorptivityIntegrandFermi(double var, MyEOSParams* my_eos_params,
                                         MyKernelParams* my_kernel_params);
MyKernelQuantity PairOpacitiesFermi(MyQuadrature* quad,
                                    MyEOSParams* my_eos_params,
                                    MyKernelParams* my_kernel_params);
#endif // GSL_INCLUDES_H

// bremsstrahlung opacities
#ifdef GSL_INCLUDES_H_
MyKernelQuantity
BremEmissivityAbsorptivityIntegrandFermi(double omega_prime,
                                         MyEOSParams* my_eos_params,
                                         MyKernelParams* my_kernel_params);
MyKernelQuantity BremOpacitiesFermi(MyQuadrature* quad,
                                    MyEOSParams* my_eos_params,
                                    MyKernelParams* my_kernel_params);
#endif // GSL_INCLUDES_H

// M1 opacities
MyQuadratureIntegrand M1CoeffsDoubleIntegrand(double* var, void* p);
MyQuadratureIntegrand M1CoeffsSingleIntegrand(double* var, void* p);
M1Opacities
ComputeM1OpacitiesNotStimulated(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params);
M1Opacities ComputeM1Opacities(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                               GreyOpacityParams* my_grey_opacity_params);
void ComputeM1DoubleIntegrandNotStimulated(
    MyQuadrature* quad_1d, GreyOpacityParams* my_grey_opacity_params, double t,
    M1Matrix* out_n, M1Matrix* out_j);
void ComputeM1DoubleIntegrand(MyQuadrature* quad_1d,
                              GreyOpacityParams* my_grey_opacity_params,
                              double t, M1Matrix* out_n, M1Matrix* out_j);
SpectralOpacities ComputeSpectralOpacitiesNotStimulatedAbs(
    const double nu, MyQuadrature* quad_1d,
    GreyOpacityParams* my_grey_opacity_params);
SpectralOpacities ComputeSpectralOpacitiesStimulatedAbs(
    const double nu, MyQuadrature* quad_1d,
    GreyOpacityParams* my_grey_opacity_params);


#endif // BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
