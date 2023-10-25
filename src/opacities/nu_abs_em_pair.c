//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file nu_abs_em_pair.c
//  \brief calculate emissivity and inverse mean free path for the pair process

#include <math.h>
#include "bns_nurates.h"
#include "constants.h"
#include "kernels.h"
#include "functions.h"
#include "integration.h"

// Calculate the 4 integrands needed for computing pair opacities (assuming no dependence on mu' and phi')
// var[3]: [omega' mu' phi'] can in generally have 3 variables but the kernel we are considering has already
// been integrated over mu' and phi' and therefore only depends on the anti-neutrino energy (nu') and EOS/Kernel
// parameters.
MyKernelQuantity PairEmissivityAbsorptivityIntegrandFermi(double *var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  my_kernel_params->pair_kernel_params.omega_prime = var[0];
  MyKernelQuantity pair_kernel = PairKernelsPhiMuIntegrated(my_eos_params, &my_kernel_params->pair_kernel_params);

  double fermi_e = FermiDistr(var[0], my_eos_params->temp, my_eos_params->mu_e / my_eos_params->temp);
  double fermi_x = FermiDistr(var[0], my_eos_params->temp, 0.);

  MyKernelQuantity result;

  result.em_e = var[0] * var[0] * (1. - fermi_e) * pair_kernel.em_e;
  result.em_x = var[0] * var[0] * (1. - fermi_x) * pair_kernel.em_x;

  result.abs_e = var[0] * var[0] * fermi_e * pair_kernel.abs_e;
  result.abs_x = var[0] * var[0] * fermi_x * pair_kernel.abs_x;

  return result;

}

MyKernelQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  double t = 1.5 * my_eos_params->temp;
  MyFunctionSpecial func;
  func.function = &PairEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyKernelQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}