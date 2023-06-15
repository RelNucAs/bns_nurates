//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file kernels.h
//  \brief header file for kernels, kernel parameters and associated functions

#ifndef BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
#define BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_

#include "../../bns_nurates.h"

// @TODO: Rename elastic scattering stuff

// ===============================================================================
// (1) Bremsstrahlung kernel

struct BremKernelParams {
  double omega;
  double omega_prime;
  double m_N;
};
typedef struct BremKernelParams BremKernelParams;

// bremsstrahlung helper functions and kernels
double BremKernelS(double x, double y, double eta_star);
double BremKernelG(double y, double eta_star);

MyKernel BremKernels(BremKernelParams *kernel_params, MyEOSParams *eos_params);

// End of Bremsstrahlung kernel
// ===============================================================================


// ===============================================================================
// (2) Pair process kernel

// pair kernel specific parameters
struct PairKernelParams {
  double omega;
  double omega_prime;
  double cos_theta;
  double mu;
  double mu_prime;
  double lmax;
  double filter;
};
typedef struct PairKernelParams PairKernelParams;

// pair helper functions and kernels
double PairT(int l, double alpha, double tolerance);
double PairF(int k, double eta, double x1);
double PairG(int n, double a, double b, double eta, double y, double z);
double PairPhi(int l, double omega, double omega_prime, double eta, double temp, int e_x);
MyKernel PairKernels(PairKernelParams *kernel_pars, MyEOSParams *eos_pars);
MyKernel PairKernelsPhiIntegrated(PairKernelParams *kernel_pars, MyEOSParams *eos_pars);
// End of pair process kernel
// ===============================================================================
/*===========================================================================*/


// ===============================================================================
// (2) Elastic scattering on nucleons kernel

// elastic scattering specific parameters
struct ElasticScattParams {
  double omega;
  double mu;
  double mu_prime;
};
typedef struct ElasticScattParams ElasticScattParams;

// bremsstrahlung helper functions and kernels

// Scattering kernel for neutrino-neutron scatyytering
MyKernel nu_p_scatt_kern(ElasticScattParams *kernel_params, MyEOSParams *eos_pars);

// Scattering kernel for neutrino-neutron scatyytering
MyKernel nu_n_scatt_kern(ElasticScattParams *kernel_params, MyEOSParams *eos_pars);

// Sum of scattering kernels for neutrino-nucleon scattering (protons + neutrons)
MyKernel nu_N_scatt_kern_tot(ElasticScattParams *kernel_params, MyEOSParams *eos_pars);

// End of elastic scattering on nucleons kernel
// ===============================================================================
/*===========================================================================*/

#endif //BNS_NURATES_SRC_OPACITIES_KERNELS_KERNELS_H_
