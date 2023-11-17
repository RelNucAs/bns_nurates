//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nes.c
//  \brief contains kernels for inelastic neutrino scattering
//         on electrons and positrons
//
// Computation of inelastic neutrino-electron and neutrino-positron
// scattering using Eq. (43) from Mezzacappa & Bruenn, ApJ v.410, p.740 (1993)
// https://ui.adsabs.harvard.edu/abs/1993ApJ...410..740M/abstract

#include "math.h"

#include "kernels.h"
#include "bns_nurates.h"
#include "functions.h"
#include "functions.h"
#include "constants.h"


#define __SIN2W__  0.2229

static const double kNes = 2. * kGf * kGf * kClight / (3. * kPi * kHbarClight* kHbarClight* kHbarClight* kHbarClight);

double kBPlus = (2. * __SIN2W__ + 1.) * (2. * __SIN2W__ + 1.);
double kBMinus = (2. * __SIN2W__ - 1.) * (2. * __SIN2W__ - 1.);
double kBZero = 4. * __SIN2W__ * __SIN2W__;

void ComputeFDIForInelastic(double w, double wp, double eta, double* fdi_diff_w, double* fdi_diff_abs){
  double abs_val = fabs(w - wp);

  fdi_diff_w[0] = FDI_p1(eta - wp) - FDI_p1(eta - w);
  fdi_diff_w[1] = FDI_p2(eta - wp) - FDI_p2(eta - w);
  fdi_diff_w[2] = FDI_p3(eta - wp) - FDI_p3(eta - w);
  fdi_diff_w[3] = FDI_p4(eta - wp) - FDI_p4(eta - w);
  fdi_diff_w[4] = FDI_p5(eta - wp) - FDI_p5(eta - w);

  fdi_diff_abs[0] = FDI_p3(eta) - FDI_p3(eta - abs_val);
  fdi_diff_abs[1] = FDI_p4(eta) - FDI_p4(eta - abs_val);
  fdi_diff_abs[2] = FDI_p5(eta) - FDI_p5(eta - abs_val);
}

double MezzacappaIntOut(double w_, double wp_, double b1_, double b2_, double* fdi_diff_w, double* fdi_diff_abs)
{
    int sign_ = 2 * (wp_ - w_ < 0.) - 1;
        
    double x_ = (w_ > wp_) ? w_ : wp_,
        y_ = (w_ < wp_) ? w_ : wp_;

    return
        (
            (b1_ + b2_) * (sign_ * fdi_diff_abs[2] - fdi_diff_w[4]) * 0.2 // All G5 terms

            - b1_ * (w_ + wp_) * (fdi_diff_w[3] + 2. * (w_ + wp_) * fdi_diff_w[2]) // G4(eta - wp) + G3(eta - wp) term

            - 6. * b1_ * w_ * wp_ * ( (w_ + wp_) * fdi_diff_w[1] + w_ * wp_ * fdi_diff_w[0] ) // G2(eta - wp) + G1(eta - wp) term

            + sign_ * (
                (b1_ * x_ - b2_ * y_) * fdi_diff_abs[1]
                + 2. * (b1_ * x_ * x_ + b2_ * y_ * y_) * fdi_diff_abs[0]
                ) // G4(eta) + G3(eta) terms
            ) / ((1 - SafeExp(wp_ - w_)) * w_ * w_ * wp_ * wp_ + (wp_ == w_));
}

double MezzacappaIntIn(double w_, double wp_, double b1_, double b2_, double* fdi_diff_w, double* fdi_diff_abs)
{
    int sign_ = 2 * (wp_ - w_ < 0.) - 1;
        
    double x_ = (w_ > wp_) ? w_ : wp_,
        y_ = (w_ < wp_) ? w_ : wp_;

    return
        (
            (b1_ + b2_) * (sign_ * fdi_diff_abs[2] + fdi_diff_w[4]) * 0.2 // All G5 terms

            + b1_ * (w_ + wp_) * (fdi_diff_w[3] + 2. * (w_ + wp_) * fdi_diff_w[2]) // G4(eta - wp) + G3(eta - wp) term

            + 6. * b1_ * w_ * wp_ * ( (w_ + wp_) * fdi_diff_w[1] + w_ * wp_ * fdi_diff_w[0] ) // G2(eta - wp) + G1(eta - wp) term

            + sign_ * (
                (b1_ * x_ - b2_ * y_) * fdi_diff_abs[1]
                + 2. * (b1_ * x_ * x_ + b2_ * y_ * y_) * fdi_diff_abs[0]
                ) // G4(eta) + G3(eta) terms
            ) / ((1 - SafeExp(wp_ - w_)) * w_ * w_ * wp_ * wp_ + (wp_ == w_));
}

MyKernelOutput NESKernels(InelasticScattKernelParams *kernel_params, MyEOSParams *eos_params)
{
  const double t = eos_params -> temp, w = kernel_params -> omega / t,
    wp = kernel_params -> omega_prime / t;

  const double eta_e = eos_params -> mu_e / t;

  MyKernelOutput output;

  double fdi_diff_abs[3], fdi_diff_w[5];

  ComputeFDIForInelastic(w, wp, eta_e, fdi_diff_w, fdi_diff_abs);

  output.abs[id_nue] = kNes * t * t * MezzacappaIntOut(w, wp, kBPlus, kBZero, fdi_diff_w, fdi_diff_abs);
  output.abs[id_anue] = kNes * t * t * MezzacappaIntOut(w, wp, kBZero, kBPlus, fdi_diff_w, fdi_diff_abs);
  output.abs[id_nux] = kNes * t * t * MezzacappaIntOut(w, wp, kBMinus, kBZero, fdi_diff_w, fdi_diff_abs);
  //output.abs[id_anux] = kNes * MezzacappaIntOut(w, wp, kBZero, kBMinus, fdi_diff_w, fdi_diff_abs); @TODO: extend library for including anux

  output.em[id_nue] = kNes * t * t * MezzacappaIntIn(wp, w, kBPlus, kBZero, fdi_diff_w, fdi_diff_abs);
  output.em[id_anue] = kNes * t * t * MezzacappaIntIn(wp, w, kBZero, kBPlus, fdi_diff_w, fdi_diff_abs);
  output.em[id_nux] = kNes * t * t * MezzacappaIntIn(wp, w, kBMinus, kBZero, fdi_diff_w, fdi_diff_abs);
  //output.em[id_anux] = kNes * MezzacappaIntIn(wp, w, kBZero, kBMinus, fdi_diff_w, fdi_diff_abs); @TODO: extend library for including anux

  return output;
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// @TODO: define here global constants for the file
// #define kSqrtPi 1.7724538509055159 // sqrt(kPi)

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
