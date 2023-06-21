//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  nu_abs_em_brem.c
//  \brief computes opacities for the bremsstrahlung reaction

#include <math.h>
#include "../bns_nurates.h"
#include "../constants.h"
#include "kernels/kernels.h"
#include "../functions/functions.h"
#include "../integration/integration.h"

MyOpacityQuantity BremEmissivityAbsorptivityIntegrandFermi(double omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {
  my_kernel_params->brem_kernel_params.omega_prime = omega_prime;
  MyOpacityQuantity brem_kernel = PairKernels(my_eos_params, &my_kernel_params->pair_kernel_params); // @TODO: fixme!

  double fermi_e = FermiDistr(omega_prime, my_eos_params->temp, my_eos_params->mu_e / my_eos_params->temp);
  double fermi_x = FermiDistr(omega_prime, my_eos_params->temp, 0.);

  MyOpacityQuantity result;
  double prefactor = kGSqr * 1. / (pow(2. * kPi, 2.)); // @TODO: fixme!

  result.em_e = prefactor * omega_prime * omega_prime * (1. - fermi_e) * brem_kernel.em_e;
  result.em_x = prefactor * omega_prime * omega_prime * (1. - fermi_x) * brem_kernel.em_x;
  result.abs_e = prefactor * omega_prime * omega_prime * fermi_e * brem_kernel.abs_e;
  result.abs_x = prefactor * omega_prime * omega_prime * fermi_x * brem_kernel.abs_x;

  return result;
}

MyOpacityQuantity BremOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  double t = 1.5 * my_eos_params->temp;
  MyFunctionSpecial func;
  func.function = &BremEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyOpacityQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}