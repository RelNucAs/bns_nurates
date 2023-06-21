//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  nu_abs_em_brem.c
//  \brief computes opacities for the bremsstrahlung reaction

#include "../bns_nurates.h"
#include "kernels/kernels.h"
#include "../functions/functions.h"

MyOpacityIntegrand BremOpacityIntegrand(double omega_prime, MyEOSParams my_eos_params, BremKernelParams my_kernel_params) {

  MyKernel BremKernel = BremKernels(&my_kernel_params, &my_eos_params);

  // @TODO: fixme
  double result_production = omega_prime * omega_prime * BremKernel.production * (1. - FermiDistr(omega_prime, my_eos_params.temp, my_eos_params.mu_e));
  double result_absorption = omega_prime * omega_prime * BremKernel.production * (1. - FermiDistr(omega_prime, my_eos_params.temp, my_eos_params.mu_e));

  MyOpacityIntegrand result = {.em = result_production, .ab = result_absorption};

  return result;
}