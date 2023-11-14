//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nes.c
//  \brief contains kernels for inelastic neutrino scattering
//         on electrons and positrons
//
// Computation of inelastic neutrino-electron and neutrino-positron
// scattering using .... from ....
// @TODO: add paper from which kernel is taken

// @TODO: inclusion of headers here
#include "kernels.h"
#include "bns_nurates.h"


// @TODO: define here global constants for the file
// #define kSqrtPi 1.7724538509055159 // sqrt(kPi)

// @TODO: write kernel implementation here ...
MyKernelOutput NESKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params) {
  MyKernelOutput nes_kernel = {0.};

  return nes_kernel;
}

MyKernelOutput NPSKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params) {
  MyKernelOutput nps_kernel = {0.};
  
  return nps_kernel;
}

MyKernelOutput InelasticScattKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params) {
  MyKernelOutput nes_kernel = NESKernels(kernel_params, eos_params);
  MyKernelOutput nps_kernel = NPSKernels(kernel_params, eos_params);

  MyKernelOutput tot_kernel = {0.};

  for (int idx = 0; idx<total_num_species; idx++) {
    tot_kernel.em[idx] = nes_kernel.em[idx] + nps_kernel.em[idx];
    tot_kernel.abs[idx] = nes_kernel.abs[idx] + nps_kernel.abs[idx];
  }

  return tot_kernel;
}