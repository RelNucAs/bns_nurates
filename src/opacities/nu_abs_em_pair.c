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

MyOpacityQuantity PairEmissivityAbsorptivityIntegrandFermi(double *var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  my_kernel_params->pair_kernel_params.omega_prime = var[0];
  MyOpacityQuantity pair_kernel = PairKernelsPhiMuIntegrated(&my_kernel_params->pair_kernel_params, my_eos_params);

  double fermi_e = FermiDistr(var[0], my_eos_params->temp, my_eos_params->mu_e / my_eos_params->temp);
  double fermi_x = FermiDistr(var[0], my_eos_params->temp, 0.);

  MyOpacityQuantity result;

  result.em_e = var[0] * var[0] * (1. - fermi_e) * pair_kernel.em_e;
  result.em_x = var[0] * var[0] * (1. - fermi_x) * pair_kernel.em_x;
  result.abs_e = var[0] * var[0] * fermi_e * pair_kernel.abs_e;
  result.abs_x = var[0] * var[0] * fermi_x * pair_kernel.abs_x;

  return result;
}

MyOpacityQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  //double t = 1.5 * my_eos_params->temp;
  double t = 1.;
  MyFunctionSpecial func;
  func.function = &PairEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyOpacityQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}