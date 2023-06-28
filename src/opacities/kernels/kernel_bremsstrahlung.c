//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_bremsstrahlung.c
//  \brief contains bremsstrahlung kernels and associated helper functions
//
// Computation of nucleon-nucleon bremsstrahlung kernel using the analytic
// fitting formula in Hannestad & Raffelt 1998, Apj, 507, 339
// (https://iopscience.iop.org/article/10.1086/306303/pdf)


#include <math.h>
#include "kernels.h"
#include "../../constants.h"
#include "../../bns_nurates.h"
#include "../../functions/functions.h"

// Constants for the bremsstrahlung reaction
static const double kBremXmin   = 1.e-10;
static const double kBremYmin   = 1.e-10;
static const double kBremEtamin = 1.e-10;
static const double kPiSquared = kPi * kPi;
static const double kHSquared  = kH  * kH;
static const double kFiveThirds = 5./3.;
static const double kFiveSixths = 5./6.;
static const double kCa = -1.26 / 2.;
// @TODO: decide how to define the following constant
static const double kGfBrem = 1.1663787e-11; //MeV-2
static const double kMnGrams = kMn * kMeV / (kClight * kClight); // neutron mass in grams
static const double kMpGrams = kMp * kMeV / (kClight * kClight); // proton mass in grams

/* Compute the analytical fit for the s-component of the kernel for
 * neutrino bremsstrahlung and inelastic scattering in a nucleon field
 *
 * Inputs:
 *      x:        rescaled total neutrino energy (w+wp/T)
 *      y:        pion mass parameter defined in Eqn. (38)
 *      eta_star: nucleon degeneracy parameter
 *
 * Output:
 *      s:  a dimensionless quantity as defined in Eqn. (49)
 */
double BremKernelS(double x, double y, double eta_star) {

  // prevent singular behaviour
  if (x >= 0.) {
    x = (x > kBremXmin) ? x : kBremXmin;
  } else {
    x = (x < -kBremXmin) ? x : kBremXmin;
  }

  y = (y > kBremYmin) ? y : kBremYmin;
  eta_star = (eta_star > kBremEtamin) ? eta_star : kBremEtamin;

  // compute non-degenerate approximation, s_nd in Eqn. (45)
  // @TODO: compute kSqrtPi and kPi2OneEighth only once
  const double kSqrtPi = sqrt(kPi);
  const double kPi2OneEighth = pow(kPi,1./8.);
  double s_nd_numerator = 2. * kSqrtPi * pow(x + 2. - exp(-y / 12.), 1.5) * (x * x + 2. * x * y + kFiveThirds * y * y + 1.);
  double s_nd_denominator = kSqrtPi + pow(kPi2OneEighth + x + y, 4.);
  double s_nd = s_nd_numerator / s_nd_denominator;

  // detail balance behaviour in Eqn. (41)
  if (x < 0.) {
    s_nd = s_nd * SafeExp(-x);
  }

  // compute degenerate approximation, s_d in Eqn. (46)
  double u     = sqrt(y / (2. * eta_star)) + 1.e-10;
  double u2    = u * u;
  double u_arg = u2 / (2. * sqrt(2. * u2 + 4.));
  double f_u   = 1. - kFiveSixths * u * atan(2. / u) + u2 / (3. * (u2 + 4.)) + atan(1. / u_arg) * u_arg / 3.;

  // @TODO: compute pow(0.5 * kPi, 2.5) constant only once
  double s_d = 3. * pow(0.5 * kPi, 2.5) * pow(eta_star, -2.5) * (x * x + 4. * kPiSquared) * x * f_u / (4. * kPiSquared * (1. - SafeExp(-x)));
  
  //if (s_d < 0.) {
  //  printf("s_D = %.5e\n"              , s_d);
  //  printf("one minus exp(-x) = %.5e\n", 1.-SafeExp(-x));
  //  printf("x = %.5e\n"                , x);
  //  printf("y = %.5e\n"                , y);
  //  printf("u = %.5e\n"                , u);
  //  printf("eta_star = %.5e\n"         , eta_star);
  //  printf("f_u = %.5e\n"              , f_u);
  //}

  // F, Eqn. (50)
  double f_denominator = (3. + pow(x - 1.2, 2.) + pow(x, -4.)) * (1. + eta_star * eta_star) * (1. + pow(y, 4.));
  double f_brem = 1. + 1. / f_denominator; //Eq.(50)

  // @TODO: compute pow(x, 1.1) only once
  // G, Eqn. (50)
  double g_brem = 1. - 0.0044 * pow(x, 1.1) * y / (0.8 + 0.06 * pow(y, 1.05)) * sqrt(eta_star) / (eta_star + 0.2);

  // h and C, Eqn. (50)
  double h_brem = 0.1 * eta_star / (2.39 + 0.1 * pow(eta_star, 1.1));
  double c_brem = 1.1 * pow(x, 1.1) * h_brem / (2.3 + h_brem * pow(x, 0.93) + 0.0001 * pow(x, 1.2)) * 30. / (30. + 0.005 * pow(x, 2.8));

  // p, Eqn. (50)
  double p_brem = 0.67 + 0.18 * pow(y, 0.4);

  // @TODO: decide what to do here in case of errors
  if (s_nd < 0.) printf("\ns_ND = %.5e\n", s_nd);
  if (s_d  < 0.)  printf(  "s_D  = %.5e\n", s_d);

  // interpolated formula for s in Eqn. (49)
  double s_brem = pow(pow(s_nd, -p_brem) + pow(s_d, -p_brem), -1. / p_brem) * f_brem * (1. + c_brem * g_brem);

  return s_brem;
}

/* Compute the analytic fit of the g-component of the kernel for neutrino
 * bremsstrahlung and inelastic scattering in a nucleon field. Implements
 * Eqn. (52) of Hannestad & Raffelt (1998).
 *
 * Inputs:
 *    y:        pion mass parameter defined in Eqn. (38)
 *    eta_star: nucleon degeneracy parameter
 *
 * Output:
 *    g: a dimensionless quantity as defined in Eqn. (52)
 */
double BremKernelG(double y, double eta_star) {

  // prevent singular behavior
  y = (y > kBremYmin) ? y : kBremYmin;
  eta_star = (eta_star > kBremEtamin) ? eta_star : kBremEtamin;

  // alpha_1, Eqn. (53)
  double y2 = y * y;
  double eta_star_inv = 1. / eta_star;
  double alpha_1_denom = 25. * y2 + 1.;

  double alpha_1 = (0.5 + eta_star_inv) / (1. + eta_star_inv) * (1. / alpha_1_denom) + (0.5 + eta_star / 15.6) * 25. * y2 / alpha_1_denom;

  // alpha_2, Eqn. (53)
  double alpha_2 = (0.63 + 0.04 * pow(eta_star, 1.45)) / (1. + 0.02 * pow(eta_star, 2.5));

  // @TODO: compute pow(eta_star, 1.5) only once
  // alpha_3, Eqn. (53)
  double alpha_3 = 1.2 * SafeExp(0.6 * eta_star - 0.4 * pow(eta_star, 1.5));

  // p_1, Eqn. (53)
  double p_1 = (1.8 + 0.45 * eta_star) / (1. + 0.15 * pow(eta_star, 1.5));

  // p_2, Eqn. (53)
  double p_2 = 2.3 - 0.05 * eta_star / (1. + 0.025 * eta_star);

  // g, Eqn. (52)
  double g = (alpha_1 + alpha_2 * pow(y, p_1)) / (1. + alpha_3 * pow(y, p_2) + alpha_2 * pow(y, p_1 + 2.) / 13.75);

  return g;
}

/* Compute the absorption kernels for a given NN Bremsstrahlung channel */
double BremSingleChannelAbsKernel(double n_nuc, double m_nuc, BremKernelParams *kernel_params, MyEOSParams *eos_params) {

  // EOS parameters
  double temp = eos_params->temp;

  // kernel parameters
  double omega = kernel_params->omega;
  double omega_prime = kernel_params->omega_prime;

  // temperature in units of 10 MeV
  double temp_10 = temp / 10.;

  // @TODO: optimize computation of constants
  // nucleon effective degeneracy parameter, Eqn. (36) using baryon number density instead of matter density
  double eta_star = pow(3. * kPiSquared * n_nuc, 2. / 3.) * (kHSquared * kMeV / (4. * kPiSquared)) / (2. * m_nuc * temp);

  // spin-fluctuation rate gamma, Eqn. (37)
  double gamma = 1.63 * pow(eta_star, 3. / 2.) * temp_10;

  // pion mass parameter y, Eqn. (38)
  double y = 1.94 / temp_10;

  // dimensionless neutrino energy sum
  double x = (omega + omega_prime) / temp;

  // dimensionless fitting parameter s, N.B.: x and y values are not changed by the function
  const double sb = BremKernelS(x, y, eta_star);

  // dimensionless fitting parameter g, N.B.: x and y values are not changed by the function
  const double gb = BremKernelG(y, eta_star);

  // differential absorption kernel, Eqn. (35)
  double s_kernel_abs = gamma / (x * x + pow(0.5 * gamma * gb, 2.)) * sb / temp;

  return s_kernel_abs;
}

/* Compute the production and absorption kernels for theBremsstrahlung reactions by summing
 the contributions of all NN channels */
MyKernel BremKernels(BremKernelParams *kernel_params, MyEOSParams *eos_params) {
  // @TODO: compute the following constant only once
  const double kBremConst = kClight * pow(kHbar * kClight, 5.) * kCa * kCa * kGfBrem * kGfBrem;

  // EOS parameters
  double nb   = eos_params->nb;
  double temp = eos_params->temp;
  double xn   = eos_params->yn;
  double xp   = eos_params->yp;

  double x_mean = sqrt(xn*xp);

  // kernel parameters
  double omega = kernel_params->omega;
  double omega_prime = kernel_params->omega_prime;

  // dimensionless neutrino energy sum
  double x = (omega + omega_prime) / temp;

  // compute single channel kernels
  double s_abs_nn = BremSingleChannelAbsKernel(nb*xn, kMnGrams, kernel_params, eos_params);
  double s_abs_pp = BremSingleChannelAbsKernel(nb*xp, kMpGrams, kernel_params, eos_params);
  // @TODO: compute sqrt(kMnGrams*kMpGrams) only once
  double s_abs_np = BremSingleChannelAbsKernel(nb*x_mean, sqrt(kMnGrams*kMpGrams), kernel_params, eos_params);
 
  // @TODO: compute 28./3. only once
  // total absorption kernel
  double s_abs_tot = kBremConst * nb * (xn * s_abs_nn + xp * s_abs_pp + 28./3. * x_mean * s_abs_np);
  // total production kernel (from detailed balance)
  double s_em_tot  = s_abs_tot * SafeExp(-x);
  double production_e;
  double absorption_e;
  double production_x;
  double absorption_x;
  double production;
  double absorption;
  MyKernel brem_kernel = {.absorption_e = s_abs_tot,
                          .production_e = s_em_tot,
                          .absorption_x = s_abs_tot,
                          .production_x = s_em_tot};

  return brem_kernel;
}