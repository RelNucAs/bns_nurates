//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file kernels.h
//  \brief header file for kernels, kernel parameters and associated functions

#ifndef BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
#define BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_

#include <stdbool.h>

#include "bns_nurates.h"

#define tol_t 1.0E-06

// ===============================================================================
// (1) Bremsstrahlung kernel

// bremsstrahlung helper functions and kernels
double BremKernelS(double x, double y, double eta_star);
double BremKernelG(double y, double eta_star);
MyKernelOutput BremKernelsLegCoeff(BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params);
void BremKernelsTable(const int n, double* nu_array,
                      GreyOpacityParams* grey_pars, M1Matrix* out);
MyKernelOutput BremKernelsBRT06(BremKernelParams* kernel_params,
                                MyEOSParams* eos_pars);
void BremKernelsTableBRT06(const int n, double* nu_array,
                           GreyOpacityParams* grey_pars, M1Matrix* out);

// End of Bremsstrahlung kernel
// ===============================================================================


// ===============================================================================
// (2) Pair process kernel


// pair helper functions and kernels


MyKernelOutput PairKernelsOptimized(MyEOSParams* eos_pars,
                                    PairKernelParams* kernel_pars);
void PairKernelsTable(const int n, double* nu_array,
                      GreyOpacityParams* grey_pars, M1Matrix* out);

#ifdef GSL_INCLUDES_H_
void PairTInterpolated(PairKernelParams* kernel_pars, double alpha,
                       double* out);
MyKernelQuantity PairKernelsPhiMuIntegrated(MyEOSParams* eos_pars,
                                            PairKernelParams* kernel_pars);
void PairKernelsM1Test(MyEOSParams* eos_pars, PairKernelParams* kernel_pars,
                       MyKernelOutput* out_for, MyKernelOutput* out_inv);
double PairPhi(int l, double omega, double omega_prime, double eta, double temp,
               int e_x);
MyKernelQuantity PairKernelsM1(MyEOSParams* eos_pars,
                               PairKernelParams* kernel_pars);
MyKernelQuantity PairKernels(MyEOSParams* eos_pars,
                             PairKernelParams* kernel_pars);
double PairT(int l, double alpha, double tolerance);
void TabulatePairTFunction(double xi, double xf, double x[dim_pair_t],
                           double t[6][dim_pair_t]);
double PairTFitted(int l, double alpha);
double PairF(int k, double eta, double x1);
double PairFOptimized(int k, double eta, double x1);
double PairFBackup(int k, double eta, double x1);
double PairG(int n, double a, double b, double eta, double y, double z);
double PairPsi(int l, double y, double z, double eta);
#endif // GSL_INCLUDES_H_

// End of pair process kernel
// ===============================================================================
/*===========================================================================*/


// ===============================================================================
// (3) Neutrino-electron (positron) scattering kernel

// @TODO: add here functions from src/opacities/kernels/kernel_nes.c
MyKernelOutput NESKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params);
MyKernelOutput NPSKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params);
MyKernelOutput InelasticScattKernels(InelasticScattKernelParams* kernel_params,
                                     MyEOSParams* eos_params);
void InelasticKernelsTable(const int n, double* nu_array,
                           GreyOpacityParams* grey_pars, M1Matrix* out);

// End of pair process kernel
// ===============================================================================
/*===========================================================================*/


#endif // BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
