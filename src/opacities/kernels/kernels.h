//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file kernels.h
//  \brief header file for kernels and associated functions

#ifndef BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
#define BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_

#include "../../bns_nurates.h"

// Parameters which are needed by a kernel (which does not come from an EOS table)
struct PairKernelParams {
  double omega;
  double omega_prime;
  double cos_theta;
  double lmax;
  double filter;
};
typedef struct PairKernelParams PairKernelParams;

// bremsstrahlung helper functions and kernels
double BremKernelS(double x, double y, double eta_star);
double BremKernelG(double y, double eta_star);

MyKernel BremKernels(MyEOSParams *params);

// pair helpher functions and kernels
double PairT(int l, double alpha, double tolerance);
double PairF(int k, double eta, double x1);
double PairG(int n, double a, double b, double eta, double y, double z);
double PairPhi(int l, double omega, double omega_prime, double eta, double temp, int e_x);
MyKernel PairKernels(PairKernelParams *kernel_pars, MyEOSParams *eos_pars);
#endif //BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
