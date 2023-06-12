//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.c
//  \brief contains pair kernels and associated helper functions

#include <math.h>
#include <assert.h>
#include <stddef.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_zeta.h>
#include "kernels.h"
#include "../../constants.h"

/* Compute T_l(alpha) as defined in Appendix B of Pons et. al. (1998)
 * upto a given accuracy
 *
 * T_l(alpha) = sum_{n=1..inf} (-1)^{n+1} e^{-n alpha} / n^l
 *
 * This is well-defined for alpha > 0 and l >= 1.
 *
 * Inputs:
 *    l:          mode number (>= 1)
 *    alpha:      function argument (> 0)
 *    tolerance:  a measure of accuracy of the computation (1e-6 is a good number)
 *
 * Output:
 *    T_l(alpha) = sum_{n=1..inf} (-1)^{n+1} e^{-n alpha} / n^l
 */
double PairT(int l, double alpha, double tolerance) {

  assert(alpha >= 0 && l >= 1);

  if (alpha == 0 && l != 1) {
    return pow(2., -l) * (pow(2., l) - 2.) * gsl_sf_zeta_int(l);  // Computed in Mathematica
  } else if (alpha == 0 && l == 1) {
    return log(2.0);
  } else {
    double val = fabs(42. * tolerance);
    double result = 0.;
    int n = 1;

    while (tolerance < fabs(val)) {
      val = pow(-1., n + 1) * exp(-n * alpha) / pow(n, l);
      result += val;
      n++;
    }
    return result;
  }
}

/* Compute F_k(eta,x1) as defined in Appendix B of Pons et. al. (1998)
 *
 * Inputs:
 *    k:    an integer
 *    eta:  double
 *    x1:   double
 *
 * Output:
 *  F_k:    as defined below
 *
 * eta < 0:         F_k(eta,x1) = k! [T_{k+1}(-eta) - sum_{l=0..k} T_{k+1-l}(x1-eta) x1^l / l!]
 * 0 <= eta <= x1:  F_k(eta,x1) = eta^{k+1}/(k+1)
 *                              + k! [2 sum_{l=0..int((k-1)/2)} T_{2l+2}(0) eta^{k-1-2l}/(k-1-2l)!
 *                                    + (-1)^k T_{k+1}(eta)
 *                                    - sum_{l=0..k} T_{k+1-l}(x1-eta) x1^l / l!]
 * x1 < eta:        F_k(eta,x1) = x1^{k+1}/(k+1)
 *                              + k! [(-1)^k T_{k+1}(eta)
 *                                     - sum_{l=0..k} (-1)^{k-l} T_{k+1-l}(eta-x1) x1^l / l!]
 */
double PairF(int k, double eta, double x1) {

  double result = 0.;
  double tol = 1e-6;

  if (eta < 0.) {
    double sum = 0.;
    for (int l = 0; l <= k; l++) {
      sum += PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
    }
    result = tgamma(k + 1) * (PairT(k + 1., -eta, tol) - sum);

  } else if (k != 0 && 0. <= eta && eta <= x1) {

    double sum = 0.;
    for (int l = 0; l <= (int) ((k - 1.) / 2.); l++) {
      sum += PairT(2 * l + 2, 0., tol) * pow(eta, k - 1. - 2. * l) / tgamma(k - 1 - 2 * l + 1);
    }

    double sum_2 = 0.;
    for (int l = 0; l <= k; l++) {
      sum_2 += PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
    }

    result = pow(eta, k + 1.) / (k + 1.) + tgamma(k + 1) * (2.0 * sum + pow(-1., k) * PairT(k + 1, eta, tol) - sum_2);

  } else if (k == 0 && 0. <= eta && eta <= x1) {

    double sum = 0.;
    for (int l = 0; l <= k; l++) {
      sum += PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
    }

    result = pow(eta, k + 1.) / (k + 1.) + tgamma(k + 1) * (pow(-1., k) * PairT(k + 1, eta, tol) - sum);

  } else {

    double sum = 0.;
    for (int l = 0; l <= k; l++) {
      sum += pow(-1., k - l) * PairT(k + 1 - l, eta - x1, tol) * pow(x1, l) / tgamma(l + 1);
    }
    result = pow(x1, k + 1.) / (k + 1.) + tgamma(k + 1) * (pow(-1., k) * PairT(k + 1, eta, tol) - sum);

  }

  return result;
}

/* Calculate G_n(a,b) from Eqn. (12) of Pons et. al. (1998)
 *
 * Inputs:
 *    n:            integer
 *    a,b,eta,y,z:  double
 *
 * Output:
 *    G_n(a,b,eta,y,z) = F_n(eta,b) - F_n(eta,a) - F_n(eta+y+z,b) + F_n(eta+y+z,a)
 * */
double PairG(int n, double a, double b, double eta, double y, double z) {

  double result = 0.;
  result = PairF(n, eta, b) - PairF(n, eta, a) - PairF(n, eta + y + z, b) + PairF(n, eta + y + z, a);

  return result;
}

/* Calculate Psi_l(y,z) from Eqn. (11) of Pons et. al.
 *
 * Inputs:
 *      l:    integer (mode number)
 *      y:    dimensionless energy
 *      z:    dimensionless energy
 *      eta:  electron degenracy parameter
 *
 * Output:
 *      Psi_l(y,z) = sum_{n=0..2} [c_{ln} G_n(y,y+z) + d_{ln} G_n(z,y+z)] + sum_{n=3..(2l+5)} a_{ln} [G_n(0,min(y,z)) - G_n(max(y,z),y+z)]
 */
double PairPsi(int l, double y, double z, double eta) {

  assert(0 <= l && l <= 3);

  double a[4][12], c[4][3], d[4][3];

  a[0][3] = 8. / (3. * y * y);
  a[0][4] = -4. / (3. * y * y * z);
  a[0][5] = 4. / (15. * y * y * z * z);
  c[0][0] = (4. * y / (z * z)) * (2. * z * z / 3. + y * z + 2. * y * y / 5.);
  c[0][1] = -(4. * y / (3. * z * z)) * (3. * y + 4. * z);
  c[0][2] = 8. * y / (3. * z * z);
  d[0][0] = 4. * z * z * z / (15. * y * y);
  d[0][1] = -4. * z * z / (3. * y * y);
  d[0][2] = 8. * z / (3. * y * y);

  // Coefficients for l = 1
  a[1][3] = 8. / (3. * y * y);
  a[1][4] = -4. * (4. * y + 3. * z) / (3. * y * y * y * z);
  a[1][5] = 4. * (13. * y + 18. * z) / (15. * y * y * y * z * z);
  a[1][6] = -4. * (y + 3. * z) / (5. * y * y * y * z * z * z);
  a[1][7] = 16. / (35. * y * y * y * z * z * z);
  c[1][0] = -(4. * y / (z * z * z)) *
      (2. * y * y * y / 7. + 4. * y * y * z / 5. + 4. * y * z * z / 5. + z * z * z / 3.);
  c[1][1] = (4. * y / (z * z * z)) * (4. * y * y / 5. + 7. * y * z / 5. + 2. * z * z / 3.);
  c[1][2] = -(4. * y / (z * z * z)) * (3. * y / 5. + z / 3.);
  d[1][0] = -(4. * z * z * z / (105. * y * y * y)) * (14. * y + 9. * z);
  d[1][1] = (4. * z * z / (5. * y * y * y)) * (7. * y / 3. + 2. * z);
  d[1][2] = -(4. * z / (y * y * y)) * (y / 3. + 3. * z / 5.);

  // Coefficients for l = 2
  a[2][3] = 8. / (3. * y * y);
  a[2][4] = -4. * (10. * y + 9. * z) / (3. * y * y * y * z);
  a[2][5] = 4. * (73. * y * y + 126. * y * z + 36. * z * z) / (15. * y * y * y * y * z * z);
  a[2][6] = -12. * (y * y + 3. * y * z + 8. * z * z / 5.) / (y * y * y * y * z * z * z);
  a[2][7] = 48. * (2. * y * y + 13. * y * z + 12. * z * z) / (35. * y * y * y * y * z * z * z * z);
  a[2][8] = -24. * (y + 2. * z) / (7. * y * y * y * y * z * z * z * z);
  a[2][9] = 8. / (7. * y * y * y * y * z * z * z * z);
  c[2][0] = 4. * y * (2. * y * y * y * y / 7. + 6. * y * y * y * z / 7 + 32. * y * y * z * z / 35. +
      2. * y * z * z * z / 5. + z * z * z * z / 15.) / (z * z * z * z);
  c[2][1] =
      -4. * y * (6. * y * y * y / 7. + 12. * y * y * z / 7. + y * z * z + 2. * z * z * z / 15.) / (z * z * z *
          z);
  c[2][2] = 8. * y * (z * z / 6. + 3 * y * z / 2. + 12 * y * y / 7.) / (5. * z * z * z * z);
  d[2][0] = 4. * z * z * z * (16. * y * y + 27. * y * z + 12. * z * z) / (105. * y * y * y * y);
  d[2][1] = -4. * z * z * (18. * z * z / 35. + 6. * y * z / 7. + y * y / 3.) / (y * y * y * y);
  d[2][2] = 8. * z * (y * y / 6. + 3. * y * z / 2. + 12. * z * z / 7.) / (5. * y * y * y * y);

  // Coefficients for l = 3
  a[3][3] = 8. / (3. * y * y);
  a[3][4] = -4. * (19. * y + 18. * z) / (3. * y * y * y * z);
  a[3][5] = 4. * (253. * y * y + 468. * y * z + 180. * z * z) / (15. * y * y * y * y * z * z);
  a[3][6] = -8. * (149. * y * y * y + 447. * y * y * z + 330. * y * z * z + 50. * z * z * z) / (15. * y * y * y * y * y * z * z * z);
  a[3][7] = 8. * (116. * y * y * y + 2916. * y * y * z / 5. + 696. * y * z * z + 200. * z * z * z) / (21. * y * y * y * y * y * z * z * z * z);
  a[3][8] = -40. * (5. * y * y * y + 54. * y * y * z + 108. * y * z * z + 50. * z * z * z) / (21. * y * y * y * y * y * z * z * z * z * z);
  a[3][9] = 40. * (10. * y * y + 43. * y * z + 100. * z * z / 3.) / (21. * y * y * y * y * y * z * z * z * z * z);
  a[3][10] = -40. * (3. * y + 5. * z) / (9. * y * y * y * y * y * z * z * z * z * z);
  a[3][11] = 320. / (99. * y * y * y * y * y * z * z * z * z * z);
  c[3][0] = -4. * y * y * (10. * y * y * y * y / 33. + 20. * y * y * y * z / 21. + 68. * y * y * z * z / 63. + 18. * y * z * z * z / 35. + 3. * z * z * z * z / 35.)
      / (z * z * z * z * z);
  c[3][1] = 4. * y * y * (20. * y * y * y / 21. + 130. * y * y * z / 63. + 48. * y * z * z / 35. + 9. * z * z * z / 35.) / (z * z * z * z * z);
  c[3][2] = -4. * y * y * (50. * y * y / 63. + 6. * y * z / 7. + 6. * z * z / 35.) / (z * z * z * z * z);
  d[3][0] = -4. * z * z * z * (9. * y * y * y + 34. * y * y * z + 40. * y * z * z + 500. * z * z * z / 33.) / (105. * y * y * y * y * y);
  d[3][1] = 4. * z * z * (3. * y * y * y + 24. * y * y * z + 130. * y * z * z / 3. + 200. * z * z * z / 9.) / (35. * y * y * y * y * y);
  d[3][2] = -4. * z * z * (50. * z * z / 63. + 6. * z * y / 7. + 6. * y * y / 35.) / (y * y * y * y * y);

  double result = 0.;

  for (int n = 0; n <= 2; n++) {
    result += c[l][n] * PairG(n, y, y + z, eta, y, z) + d[l][n] * PairG(n, z, y + z, eta, y, z);
  }

  double min_yz = (y > z) ? z : y;
  double max_yz = (y > z) ? y : z;

  for (int n = 3; n <= 2 * l + 5; n++) {
    result += a[l][n] * (PairG(n, 0, min_yz, eta, y, z) - PairG(n, max_yz, y + z, eta, y, z));
  }

  return result;
}

/* Calculate Phi_l(y,z) from Eqn. (10) of Pons et. al. (1998)
 *
 * Inputs:
 *      l:            mode number
 *      omega:        neutrino energy [MeV]
 *      omega_prime:  anti-neutrino energy [MeV]
 *      temp:         temperature [MeV]
 *      e_x:          neutrino species type (0: elentron, 1: mu/tau)
 *
 * Output:
 *      Phi_l(y,z) = (G^2 temp^2)/(pi (1 - e^{y+z})) [alpha1 Psi_l(y,z) + alpha2 Psi_l(z,y)]
 */
double PairPhi(int l, double omega, double omega_prime, double eta, double temp, int e_x) {

  assert(e_x >= 0 && e_x <= 1);

  double y = omega / temp;
  double z = omega_prime / temp;

  double psi_1 = PairPsi(l, y, z, eta);
  double psi_2 = PairPsi(l, z, y, eta);

  double result = kGSqr * temp * temp * (kAlpha1[e_x] * kAlpha1[e_x] * psi_1 + kAlpha2[e_x] * kAlpha2[e_x] * psi_2) / (kPi * (1. - exp(y + z)));

  return result;
}

/* Calculates the production and absorption kernels for the pair process from Eqns. (2) and (3) of Pons et. al.
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega') P_l(cos theta)
 * R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T} R^p_{TP}(omega,omega',cos theta)
 *
 * l is an integer.
 * */
MyKernel PairKernels(PairKernelParams *kernel_pars, MyEOSParams *eos_pars) {

  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;
  double cos_theta = kernel_pars->cos_theta;
  int lmax = kernel_pars->lmax;
  double filterpar = kernel_pars->filter;

  double eta = eos_pars->eta;
  double temp = eos_pars->temp;

  double pair_kernel_production_e = 0.;
  double pair_kernel_absorption_e = 0.;
  double pair_kernel_production_x = 0.;
  double pair_kernel_absorption_x = 0.;
  double pair_phi_e = 0.;
  double pair_phi_x = 0.;

  assert(lmax >= 0 && lmax <= 3);

  for (int l = 0; l <= lmax; l++) {
    double legendre_l = gsl_sf_legendre_Pl(l, cos_theta);
    pair_phi_e = PairPhi(l, omega, omega_prime, eta, temp, 0);
    pair_phi_x = PairPhi(l, omega, omega_prime, eta, temp, 1);

    pair_kernel_production_e += (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) * (2. * l + 1.) * pair_phi_e * legendre_l / 2.;
    pair_kernel_production_x += (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) * (2. * l + 1.) * pair_phi_x * legendre_l / 2.;
  }

  pair_kernel_absorption_e = exp((omega + omega_prime) / temp) * pair_kernel_production_e;
  pair_kernel_absorption_x = exp((omega + omega_prime) / temp) * pair_kernel_production_x;

  MyKernel pair_kernel =
      {.absorption_e = pair_kernel_absorption_e, .production_e = pair_kernel_production_e, .absorption_x = pair_kernel_absorption_x, .production_x = pair_kernel_production_x};

  return pair_kernel;

}