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

MyOpacityQuantity BremEmissivityAbsorptivityIntegrandFermi(double *omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  MyOpacityQuantity result = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
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