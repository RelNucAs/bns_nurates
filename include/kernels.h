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
MyKernelOutput BremKernelsLegCoeff(BremKernelParams *kernel_params, MyEOSParams *eos_params);

// End of Bremsstrahlung kernel
// ===============================================================================


// ===============================================================================
// (2) Pair process kernel

void TabulatePairTFunction(double xi, double xf, double x[dim_pair_t], double t[6][dim_pair_t]);

// pair helper functions and kernels
void PairTInterpolated(PairKernelParams *kernel_pars, double alpha, double *out);
double PairTFitted(int l, double alpha);
double PairT(int l, double alpha, double tolerance);
double PairF(int k, double eta, double x1);
double PairG(int n, double a, double b, double eta, double y, double z);
double PairPsi(int l, double y, double z, double eta);
double PairPhi(int l, double omega, double omega_prime, double eta, double temp, int e_x);
MyKernelQuantity PairKernels(MyEOSParams *eos_pars, PairKernelParams *kernel_pars);
MyKernelQuantity PairKernelsM1(MyEOSParams *eos_pars, PairKernelParams *kernel_pars);
MyKernelQuantity PairKernelsM1Optimized(MyEOSParams *eos_pars, PairKernelParams *kernel_pars);
void PairKernelsM1Test(MyEOSParams *eos_pars, PairKernelParams *kernel_pars, MyKernelQuantity *out_for, MyKernelQuantity *out_inv);
MyKernelQuantity PairKernelsPhiMuIntegrated(MyEOSParams *eos_pars, PairKernelParams *kernel_pars);

// End of pair process kernel
// ===============================================================================
/*===========================================================================*/


// ===============================================================================
// (3) Neutrino-electron (positron) scattering kernel

// @TODO: add here functions from src/opacities/kernels/kernel_nes.c
MyKernelOutput NESKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params);
MyKernelOutput NPSKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params);
MyKernelOutput InelasticScattKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params);
void CrossedInelasticScattKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params, MyKernelOutput *out_1, MyKernelOutput *out_2);

// End of pair process kernel
// ===============================================================================
/*===========================================================================*/


#endif //BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
