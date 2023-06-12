//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_bremsstrahlung.c
//  \brief contains bremsstrahlung kernels and associated helper functions

#include <math.h>
#include "kernels.h"
#include "../../constants.h"
#include "../../bns_nurates.h"

/* Compute the analytical fit for the s component of the kernel for
 * neutrino bremsstrahlung and inelastic scattering in a nucleon field
 *
 * Inputs:
 *      x:        dimensionless energy omega/PairT
 *      y:        pion mass parameter defined in Eqn. (38)
 *      eta_star: neutron degeneracy parameter
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
  double s_nd_numerator = 2. * pow(kPi, 0.5) * pow(x + 2. - exp(-y / 12.), 1.5) * (x * x + 2. * x * y + 5. / 3. * y * y + 1.);
  double s_nd_denominator = pow(kPi, 0.5) + pow(pow(kPi, 1. / 8.) + x + y, 4.);
  double s_nd = s_nd_numerator / s_nd_denominator;

  // detail balance behaviour in Eqn. (41)
  if (x < 0.) {
    s_nd = s_nd * exp(-x);
  }

  // compute non-degenerate approximation, s_d in Eqn. (46)
  double u = sqrt(y / (2. * eta_star)) + 1.e-10;
  double u_arg = (u * u) / (2. * sqrt(2. * u * u + 4.));
  double f_u = 1. - 5. / 6. * u * atan(2. / u) + (u * u) / (3. * (u * u + 4.)) + atan(1. / u_arg) * u_arg / 3.;

  double s_d = 3. * pow(0.5 * kPi, 2.5) * pow(eta_star, -2.5) * (x * x + 4. * kPi * kPi) * x * f_u / (4. * kPi * kPi * (1. - exp(-x)));

  // PairF, Eqn. (50)
  double f_denominator = (3. + pow(x - 1.2, 2.) + pow(x, -4.)) * (1. + eta_star * eta_star) * (1. + pow(y, 4.));
  double f_brem = 1. + 1. / f_denominator; //Eq.(50)

  // PairG, Eqn. (50)
  double g_brem = 1. - 0.0044 * pow(x, 1.1) * y / (0.8 + 0.06 * pow(y, 1.05)) * sqrt(eta_star) / (eta_star + 0.2);

  // C, Eqn. (50)
  double h_brem = 0.1 * eta_star / (2.39 + 0.1 * pow(eta_star, 1.1));
  double c_brem = 1.1 * pow(x, 1.1) * h_brem / (2.3 + h_brem * pow(x, 0.93) + 0.0001 * pow(x, 1.2)) * 30. / (30. + 0.005 * pow(x, 2.8));

  // p, Eqn. (50)
  double p_brem = 0.67 + 0.18 * pow(y, 0.4);

  // interpolated formula for s in Eqn. (49)
  double s_brem = pow(pow(s_nd, -p_brem) + pow(s_d, -p_brem), -1. / p_brem) * f_brem * (1. + c_brem * g_brem);

  return s_brem;
}

/* Compute the analytic fit of the g component of the kernel for neutrino
 * bremsstrahlung and inelastic scattering in a nucleon field. Implements
 * Eqn. (52) of Hannestadt & Raffelt (1998).
 *
 * Inputs:
 *    y:        pion mass parameter defined in Eqn. (38)
 *    eta_star: neutron degeneracy parameter
 *
 * Output:
 *    g: a dimensionless quantity as defined in Eqn. (52)
 */
double BremKernelG(double y, double eta_star) {

  // prevent singular behavior
  y = (y > kBremYmin) ? y : kBremYmin;
  eta_star = (eta_star > kBremEtamin) ? eta_star : kBremEtamin;

  // alpha_1, Eqn. (53)
  double eta_star_inv = 1. / eta_star;
  double alpha_1_denom = 25. * y * y + 1.;

  double alpha1 = (0.5 + eta_star_inv) / (1. + eta_star_inv) * (1. / alpha_1_denom) + (0.5 + eta_star / 15.6) * 25. * y * y / alpha_1_denom;

  // alpha_2, Eqn. (53)
  double alpha_2 = (0.63 + 0.04 * pow(eta_star, 1.45)) / (1. + 0.02 * pow(eta_star, 2.5));

  // alpha_3, Eqn. (53)
  double alpha_3 = 1.2 * exp(0.6 * eta_star - 0.4 * pow(eta_star, 1.5));

  // p_1, Eqn. (53)
  double p_1 = (1.8 + 0.45 * eta_star) / (1. + 0.15 * pow(eta_star, 1.5));

  // p_2, Eqn. (53)
  double p_2 = 2.3 - 0.05 * eta_star / (1. + 0.025 * eta_star);

  // g, Eqn. (52)
  double g = (alpha1 + alpha_2 * pow(y, p_1)) / (1. + alpha_3 * pow(y, p_2) + alpha_2 * pow(y, p_1 + 2.) / 13.75);

  return g;

}
/* Compute the production and absorption kernels for the Bremsstrahlung process
 */
MyKernel BremKernels(BremKernelParams *kernel_params, MyEOSParams *eos_params) {

  // EOS parameters
  double nb = eos_params->nb;
  double temp = eos_params->temp;

  // kernel parameters
  double omega = kernel_params->omega;
  double omega_prime = kernel_params->omega_prime;
  double m_N = kernel_params->m_N;

  // temperature in units of 10 MeV
  double temp_10 = temp / 10.;

  // neutron degeneracy parameter, Eqn. (36) using baryon number density instead of matter density
  double eta_star = pow(3. * kPi * kPi * nb, 2. / 3.) * (kH * kH * kMeV / (4. * kPi * kPi)) / (2. * m_N * temp);

  // spin-fluctuation rate gamma, Eqn. (37)
  double gamma = 1.63 * pow(eta_star, 3. / 2.) * temp_10;

  // pion mass parameter y, Eqn. (38)
  double y = 1.94 / temp_10;

  // dimensionless neutrino energy sum
  double x = (omega + omega_prime) / temp;

  // dimensionless fitting parameter s, N.B.: x and y values are not changed by the function
  double sb = BremKernelS(x, y, eta_star);

  // dimensionless fitting parameter g, N.B.: x and y values are not changed by the function
  const double gb = BremKernelG(y, eta_star);

  // differential absoption kernel, Eqn. (35)
  double s_kernel_abs = gamma / (x * x + pow(0.5 * gamma * gb, 2.)) * sb / temp;

  // differential emission kernel from detail balance
  double s_kernel_prod = s_kernel_abs * exp(-x);

  MyKernel brem_kernel = {.absorption_e = s_kernel_abs, .production_e = s_kernel_prod};

  return brem_kernel;
}