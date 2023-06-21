//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file nu_abs_em_pair.c
//  \brief calculate emissivity and inverse mean free path for the pair process

#include <math.h>
#include "../bns_nurates.h"
#include "../constants.h"
#include "kernels/kernels.h"
#include "../functions/functions.h"
#include "../integration/integration.h"

MyOpacityQuantity PairEmissivityAbsorptivityIntegrandFermi(double omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {
  my_kernel_params->pair_kernel_params.omega_prime = omega_prime;
  MyOpacityQuantity pair_kernel = PairKernels(my_eos_params, &my_kernel_params->pair_kernel_params);

  double fermi_e = FermiDistr(omega_prime, my_eos_params->temp, my_eos_params->mu_e / my_eos_params->temp);
  double fermi_x = FermiDistr(omega_prime, my_eos_params->temp, 0.);

  MyOpacityQuantity result;
  double prefactor = 4. * kPi / (kClight * pow(kH * kClight, 3.));

  result.em_e = prefactor * omega_prime * omega_prime * (1. - fermi_e) * pair_kernel.em_e;
  result.em_x = prefactor * omega_prime * omega_prime * (1. - fermi_x) * pair_kernel.em_x;
  result.abs_e = prefactor * omega_prime * omega_prime * fermi_e * pair_kernel.abs_e;
  result.abs_x = prefactor * omega_prime * omega_prime * fermi_x * pair_kernel.abs_x;

  return result;
}

MyOpacityQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  double t = 1.5 * my_eos_params->temp;
  MyFunctionSpecial func;
  func.function = &PairEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyOpacityQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}