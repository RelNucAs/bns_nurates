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
#include "functions.h"
#include "constants.h"


void TabulatePairTFunction(double xi, double xf, double x[dim_pair_t], double t[6][dim_pair_t]) {
  const double dx = (xf - xi) / ((double) dim_pair_t);
  
  for (int i=0; i<dim_pair_t; i++) {
    x[i] = xi + (double) i * dx;
    for (int j=1; j<7; j++) {
      t[j-1][i] = PairT(j, x[i], tol_t) ;
    }
  }

  return;
}

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
double PairTWithAlpha0(int l) {
  static const double log2 = 0.6931471805599453; //log(2.0);
  if (l == 1) {
    return log2;
  } else {
    return pow(2., -l) * (pow(2., l) - 2.) * gsl_sf_zeta_int(l);  // Computed in Mathematica
  }
}

double PairTFitted(int l, double alpha) {
  static const double PairTFit1[6][2] = {
{-1.0000000000000557, 1.1599894050697828e-11},
{-1.0000000000000557, 1.1599894050697828e-11},
{-1.0000000000000557, 1.1599894050697828e-11},
{-1.0000000000000557, 1.1599894050697828e-11},
{-1.0000000000000557, 1.1599894050697828e-11},
{-1.0000000000000557, 1.1599894050697828e-11} };

  static const double PairTFit2[6][13] = {
{7.764460276911894e-13, -9.531900556136287e-11, 5.130720144162727e-09, -1.6018407662684595e-07, 3.2238390304643073e-06, -4.382673963153738e-05, 0.00040809016119114534, -0.0025591672605759157, 0.010059951727924867, -0.01862035702348546, -0.022454086684801058, 0.19153844068873888, 0.6932349884688382},
{1.2097121671548247e-12, -1.2360258213665597e-10, 5.628852629867178e-09, -1.5051104090763755e-07, 2.6176160114632282e-06, -3.0919945500091964e-05, 0.0002501877313487743, -0.0013457908088575755, 0.00426086049575675, -0.003453495292249957, -0.030657101372762428, 0.1290371349050982, 0.8224794727922443},
{7.17877896667645e-13, -7.111741520288681e-11, 3.1269450241639527e-09, -8.029349298289137e-08, 1.331164063502168e-06, -1.48222315936738e-05, 0.00011078575689764255, -0.0005242569477484338, 0.0011880768678479855, 0.0021057203236187075, -0.025154468942278033, 0.0791164948210426, 0.9015387115221923},
{3.228347864267874e-13, -3.1423838815286814e-11, 1.3519272756431448e-09, -3.3760229096333e-08, 5.389942800937989e-07, -5.677096233745809e-06, 3.856825962408047e-05, -0.00014559133060632152, 2.1281035907218678e-05, 0.0030071993630479362, -0.016997198961210103, 0.045559097999768754, 0.9470282577605772},
{1.223042342844013e-13, -1.1689674516055724e-11, 4.908468308315054e-10, -1.1844514209510994e-08, 1.7941662304610306e-07, -1.7228094292976135e-06, 9.471179226979555e-06, -1.096586220062312e-05, -0.0002672125800539407, 0.0023768431000457074, -0.010335683014480727, 0.025128000871717175, 0.9721171063352136},
{3.862549379641218e-14, -3.584054805259039e-12, 1.4415672453643435e-10, -3.249472337669054e-09, 4.3475802836454606e-08, -3.0974403462775056e-07, 7.929667164380111e-08, 2.257168344772742e-05, -0.00025382574940076774, 0.00154830636662188, -0.005890098985143524, 0.013449009951258391, 0.9855501299373773} }; 
  const int idx = l-1;
  double fit_poly = 1.;
  double fit_exp = exp(PairTFit1[idx][0] * alpha + PairTFit1[idx][1]);


  if (alpha < 15.) { 
    if (alpha == 0) return PairTWithAlpha0(l);
    fit_poly = PairTFit2[idx][12] + alpha *
               (PairTFit2[idx][11] + alpha *
                (PairTFit2[idx][10] + alpha *
                 (PairTFit2[idx][9] + alpha *
                  (PairTFit2[idx][8] + alpha *
                   (PairTFit2[idx][7] + alpha *
              	    (PairTFit2[idx][6] + alpha *
                     (PairTFit2[idx][5] + alpha *
                      (PairTFit2[idx][4] + alpha *
                       (PairTFit2[idx][3] + alpha *
                        (PairTFit2[idx][2] + alpha *
                         (PairTFit2[idx][1] + alpha *
                           PairTFit2[idx][0]) ) ) ) ) ) ) ) ) ) );
  }

  return fit_poly * fit_exp;
}

void PairTInterpolated(PairKernelParams *kernel_pars, double alpha, double *out) {

  assert(alpha >= 0);

  const double dx = kernel_pars->alpha[1] - kernel_pars->alpha[0]; // uniform 1d grid 
  int idx = floor( (alpha - kernel_pars->alpha[0]) / dx );

  idx = (idx < 0) ? 0 : idx;
  idx = (idx > dim_pair_t - 2) ? dim_pair_t : idx;

  assert(alpha >= kernel_pars->alpha[idx]);

  const double r = (alpha - kernel_pars->alpha[idx]) / (kernel_pars->alpha[idx+1] - kernel_pars->alpha[idx]);

  for (int i=0; i<6; i++) {
    out[i] = kernel_pars->pair_t[i][idx] * (1. - r) + kernel_pars->pair_t[i][idx+1] * r;
  }
  
  return;
}

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
      sum += PairTWithAlpha0(2 * l + 2) * pow(eta, k - 1. - 2. * l) / tgamma(k - 1 - 2 * l + 1);
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

double PairFOptimized(int k, double eta, double x1) {
  
  double result = 0.;
  int fact = 1;
  double tmp = 1.;
  
  if (eta < 0.) {
    double sum = 0.;
    sum += PairTFitted(k + 1, x1 - eta); // l=0
    for (int l = 1; l <= k; l++) {
      tmp = tmp * x1;
      fact = fact * l;   
      sum += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
    }
    result = fact * (PairTFitted(k + 1, fabs(eta)) - sum);

  } else if (k != 0 && 0. <= eta && eta <= x1) {

    double eta_to_k_minus_one = pow(eta, k - 1.);
    double eta_to_k_plus_one = eta_to_k_minus_one * eta * eta;

    double sum = 0.;
    tmp = 1.;
    sum += PairTWithAlpha0(2) / tgamma(k); // l=0
    for (int l = 1; l <= (int) ((k - 1.) / 2.); l++) {
      tmp = tmp / (eta * eta);
      sum += PairTWithAlpha0(2 * l + 2) * tmp / tgamma(k - 1 - 2 * l + 1);
    }

    double sum_2 = 0.;
    sum_2 += PairTFitted(k + 1, x1 - eta); // l=0
    fact = 1;
    tmp = 1.;
    for (int l = 1; l <= k; l++) {
      fact = fact * l;
      tmp = tmp * x1;
      sum_2 += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
    }

    result = eta_to_k_plus_one / (k + 1.) + fact * (2.0 * eta_to_k_minus_one * sum + pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum_2);

  } else if (k == 0 && 0. <= eta && eta <= x1) {

    double sum = 0.;
    sum += PairTFitted(k + 1, x1 - eta);
    for (int l = 1; l <= k; l++) {
      fact = fact * l;
      tmp = tmp * x1;
      sum += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
    }

    result = pow(eta, k + 1.) / (k + 1.) + fact * (pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum);

  } else {

    double sum = 0.;
    sum += pow(-1., k) * PairTFitted(k + 1, eta - x1);
    for (int l = 1; l <= k; l++) {
      fact = fact * l;
      tmp = tmp * x1;
      sum += pow(-1., k - l) * PairTFitted(k + 1 - l, eta - x1) * tmp / fact;
    }

    tmp = tmp * x1;
    result = tmp / (k + 1.) + fact * (pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum);

  }

  return result;
}

double PairFInterpolated(PairKernelParams *kernel_pars, int k, double eta, double x1) {
  double pair_t_1[6], pair_t_2[6];

  PairTInterpolated(kernel_pars, fabs(eta), pair_t_1); // PairT(|eta|)
  
  double result = 0.;

  if (eta < 0.) {
    double sum = 0.;
    PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1 - eta)
    for (int l = 0; l <= k; l++) {
      sum += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
    }
    result = tgamma(k + 1) * (pair_t_1[k + 1 - 1] - sum);

  } else if (k != 0 && 0. <= eta && eta <= x1) {

    double sum = 0.;
    PairTInterpolated(kernel_pars, 0., pair_t_2); // PairT(0.)
    for (int l = 0; l <= (int) ((k - 1.) / 2.); l++) {
      sum += pair_t_2[2 * l + 2 - 1] * pow(eta, k - 1. - 2. * l) / tgamma(k - 1 - 2 * l + 1);
    }

    double sum_2 = 0.;
    PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1-eta)
    for (int l = 0; l <= k; l++) {
      sum_2 += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
    }

    result = pow(eta, k + 1.) / (k + 1.) + tgamma(k + 1) * (2.0 * sum + pow(-1., k) * pair_t_1[k + 1 - 1] - sum_2);

  } else if (k == 0 && 0. <= eta && eta <= x1) {

    double sum = 0.;
    PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1-eta)
    for (int l = 0; l <= k; l++) {
      sum += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
    }

    result = pow(eta, k + 1.) / (k + 1.) + tgamma(k + 1) * (pow(-1., k) * pair_t_1[k + 1 - 1] - sum);

  } else {

    double sum = 0.;
    PairTInterpolated(kernel_pars, eta - x1, pair_t_2); // PairT(eta-x1)
    for (int l = 0; l <= k; l++) {
      sum += pow(-1., k - l) * pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
    }
    result = pow(x1, k + 1.) / (k + 1.) + tgamma(k + 1) * (pow(-1., k) * pair_t_1[k + 1 - 1] - sum);

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

void PairGOptimized(PairKernelParams *kernel_pars, int n, double eta, double y, double z, double *g_out) {
  /*
  const double pair_f_y_z_1 = PairF(n, eta, y + z);
  const double pair_f_y_z_2 = PairF(n, eta + y + z, y + z);

  const double pair_f_y_1 = PairF(n, eta, y);
  const double pair_f_y_2 = PairF(n, eta + y + z, y);

  const double pair_f_z_1 = PairF(n, eta, z);
  const double pair_f_z_2 = PairF(n, eta + y + z, z);
  */

  // /*
  const double pair_f_y_z_1 = PairFOptimized(n, eta, y + z);
  const double pair_f_y_z_2 = PairFOptimized(n, eta + y + z, y + z);

  const double pair_f_y_1 = PairFOptimized(n, eta, y);
  const double pair_f_y_2 = PairFOptimized(n, eta + y + z, y);

  const double pair_f_z_1 = PairFOptimized(n, eta, z);
  const double pair_f_z_2 = PairFOptimized(n, eta + y + z, z);
  // */

  /*
  const double pair_f_y_z_1 = PairFInterpolated(kernel_pars, n, eta, y + z);
  const double pair_f_y_z_2 = PairFInterpolated(kernel_pars, eta + y + z, y + z);

  const double pair_f_y_1 = PairFInterpolated(kernel_pars, n, eta, y);
  const double pair_f_y_1 = PairFInterpolated(kernel_pars, eta + y + z, y);

  const double pair_f_z_1 = PairFInterpolated(kernel_pars, n, eta, z);
  const double pair_f_z_2 = PairFInterpolated(kernel_pars, eta + y + z, z);
  */


  double pair_f_min_1 = (y > z) ? pair_f_z_1 : pair_f_y_1;
  double pair_f_max_1 = (y > z) ? pair_f_y_1 : pair_f_z_1;

  double pair_f_min_2 = (y > z) ? pair_f_z_2 : pair_f_y_2;
  double pair_f_max_2 = (y > z) ? pair_f_y_2 : pair_f_z_2;

  g_out[0] = pair_f_y_z_1 - pair_f_y_1 - (pair_f_y_z_2 - pair_f_y_2);
  g_out[1] = pair_f_y_z_1 - pair_f_z_1 - (pair_f_y_z_2 - pair_f_z_2);
  g_out[2] = pair_f_min_1 - pair_f_min_2;
  g_out[3] = pair_f_y_z_1 - pair_f_max_1 - (pair_f_y_z_2 - pair_f_max_2);

  return;
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

void PairPsiOptimized(PairKernelParams *kernel_pars, int l, double y, double z, double eta, double *psi_out) {

  assert(0 <= l && l <= 3);

  double a[4][12], c[4][3], d[4][3];

  const double y2 = y * y;
  const double y3 = y * y2;
  const double y4 = y * y3;
  const double y5 = y * y4;

  const double z2 = z * z;
  const double z3 = z * z2;
  const double z4 = z * z3;
  const double z5 = z * z4; 

  // @TODO: compute ratio between constant numbers (e.g 4/3) only once

  a[0][3] = 8. / (3. * y2);
  a[0][4] = -4. / (3. * y2 * z);
  a[0][5] = 4. / (15. * y2 * z2);
  c[0][0] = (4. * y / z2) * (2. * z2 / 3. + y * z + 2. * y2 / 5.);
  c[0][1] = -(4. * y / (3. * z2)) * (3. * y + 4. * z);
  c[0][2] = 8. * y / (3. * z2);
  d[0][0] = 4. * z3 / (15. * y2);
  d[0][1] = -4. * z2 / (3. * y2);
  d[0][2] = 8. * z / (3. * y2);

  // Coefficients for l = 1
  a[1][3] = 8. / (3. * y2);
  a[1][4] = -4. * (4. * y + 3. * z) / (3. * y3 * z);
  a[1][5] = 4. * (13. * y + 18. * z) / (15. * y3 * z2);
  a[1][6] = -4. * (y + 3. * z) / (5. * y3 * z3);
  a[1][7] = 16. / (35. * y3 * z3);
  c[1][0] = -(4. * y / z3) *
      (2. * y3 / 7. + 4. * y2 * z / 5. + 4. * y * z2 / 5. + z3 / 3.);
  c[1][1] = (4. * y / z3) * (4. * y2 / 5. + 7. * y * z / 5. + 2. * z2 / 3.);
  c[1][2] = -(4. * y / z3) * (3. * y / 5. + z / 3.);
  d[1][0] = -(4. * z3 / (105. * y3)) * (14. * y + 9. * z);
  d[1][1] = (4. * z2 / (5. * y3)) * (7. * y / 3. + 2. * z);
  d[1][2] = -(4. * z / y3) * (y / 3. + 3. * z / 5.);

  // Coefficients for l = 2
  a[2][3] = 8. / (3. * y2);
  a[2][4] = -4. * (10. * y + 9. * z) / (3. * y3 * z);
  a[2][5] = 4. * (73. * y2 + 126. * y * z + 36. * z * z) / (15. * y4 * z * z);
  a[2][6] = -12. * (y2 + 3. * y * z + 8. * z2 / 5.) / (y4 * z3);
  a[2][7] = 48. * (2. * y2 + 13. * y * z + 12. * z2) / (35. * y4 * z4);
  a[2][8] = -24. * (y + 2. * z) / (7. * y4 * z4);
  a[2][9] = 8. / (7. * y4 * z4);
  c[2][0] = 4. * y * (2. * y4 / 7. + 6. * y3 * z / 7 + 32. * y2 * z2 / 35. +
      2. * y * z3 / 5. + z4 / 15.) / z4;
  c[2][1] =
      -4. * y * (6. * y3 / 7. + 12. * y2 * z / 7. + y * z2 + 2. * z3 / 15.) / z4;
  c[2][2] = 8. * y * (z2 / 6. + 3 * y * z / 2. + 12 * y2 / 7.) / (5. * z4);
  d[2][0] = 4. * z3 * (16. * y2 + 27. * y * z + 12. * z2) / (105. * y4);
  d[2][1] = -4. * z2 * (18. * z2 / 35. + 6. * y * z / 7. + y2 / 3.) / y4;
  d[2][2] = 8. * z * (y2 / 6. + 3. * y * z / 2. + 12. * z2 / 7.) / (5. * y4);

  // Coefficients for l = 3
  a[3][3] = 8. / (3. * y2);
  a[3][4] = -4. * (19. * y + 18. * z) / (3. * y3 * z);
  a[3][5] = 4. * (253. * y2 + 468. * y * z + 180. * z2) / (15. * y4 * z2);
  a[3][6] = -8. * (149. * y3 + 447. * y2 * z + 330. * y * z2 + 50. * z3) / (15. * y5 * z3);
  a[3][7] = 8. * (116. * y3 + 2916. * y2 * z / 5. + 696. * y * z2 + 200. * z3) / (21. * y5 * z4);
  a[3][8] = -40. * (5. * y3 + 54. * y2 * z + 108. * y * z2 + 50. * z3) / (21. * y5 * z5);
  a[3][9] = 40. * (10. * y2 + 43. * y * z + 100. * z2 / 3.) / (21. * y5 * z5);
  a[3][10] = -40. * (3. * y + 5. * z) / (9. * y5 * z5);
  a[3][11] = 320. / (99. * y5 * z5);
  c[3][0] = -4. * y2 * (10. * y4 / 33. + 20. * y3 * z / 21. + 68. * y2 * z2 / 63. + 18. * y * z3 / 35. + 3. * z4 / 35.)
      / z5;
  c[3][1] = 4. * y2 * (20. * y3 / 21. + 130. * y2 * z / 63. + 48. * y * z2 / 35. + 9. * z3 / 35.) / z5;
  c[3][2] = -4. * y2 * (50. * y2 / 63. + 6. * y * z / 7. + 6. * z2 / 35.) / z5;
  d[3][0] = -4. * z3 * (9. * y3 + 34. * y2 * z + 40. * y * z2 + 500. * z3 / 33.) / (105. * y5);
  d[3][1] = 4. * z2 * (3. * y3 + 24. * y2 * z + 130. * y * z2 / 3. + 200. * z3 / 9.) / (35. * y5);
  d[3][2] = -4. * z2 * (50. * z2 / 63. + 6. * z * y / 7. + 6. * y2 / 35.) / y5;

  double aux;

  double result_y_z = 0., result_z_y = 0.;

  double pair_g_y, pair_g_z;

  double pair_g[4] = {0.};


  for (int n = 0; n <= 2; n++) {
    PairGOptimized(kernel_pars, n, eta, y, z, pair_g);
    pair_g_y = pair_g[0]; //PairG(n, y, y + z, eta, y, z); // a = y, b = y + z
    pair_g_z = pair_g[1]; //PairG(n, z, y + z, eta, y, z); // a = z, b = y + z
    result_y_z += c[l][n] * pair_g_y + d[l][n] * pair_g_z;
    result_z_y += d[l][n] * pair_g_y + c[l][n] * pair_g_z;
  }

  double min_yz = (y > z) ? z : y;
  double max_yz = (y > z) ? y : z;

  for (int n = 3; n <= 2 * l + 5; n++) {
    PairGOptimized(kernel_pars, n, eta, y, z, pair_g);
    aux = a[l][n] * (pair_g[2] - pair_g[3]);
    //aux = a[l][n] * (PairG(n, 0, min_yz, eta, y, z) - PairG(n, max_yz, y + z, eta, y, z));
    result_y_z += aux;
    result_z_y += aux;
  }
  
  psi_out[0] = result_y_z; // Psi(y,z)
  psi_out[1] = result_z_y; // Psi(z,y)

  return;
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

void PairPhiOptimized(PairKernelParams *kernel_pars, int l, double eta, double temp, double *phi_out) {
  static const double kPairPhi = kGSqr / kPi;

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;

  const double y = omega / temp;
  const double z = omega_prime / temp;

  const double phi_prefactor = kPairPhi * temp * temp;
  const double phi_denom = 1. - exp(y + z);

  double pair_psi[2] = {0.};

  PairPsiOptimized(kernel_pars, l, y, z, eta, pair_psi);
  
  phi_out[0] = phi_prefactor * (kAlpha1[0] * kAlpha1[0] * pair_psi[0] + kAlpha2[0] * kAlpha2[0] * pair_psi[1]) / phi_denom;
  phi_out[1] = phi_prefactor * (kAlpha1[1] * kAlpha1[1] * pair_psi[0] + kAlpha2[1] * kAlpha2[1] * pair_psi[1]) / phi_denom;

  return;
}

/* Calculates the production and absorption kernels for the pair process from Eqns. (2) and (3) of Pons et. al.
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega') P_l(cos theta)
 * R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T} R^p_{TP}(omega,omega',cos theta)
 *
 * l is an integer.
 */
MyKernelQuantity PairKernels(MyEOSParams *eos_pars, PairKernelParams *kernel_pars) {

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;
  double cos_theta = kernel_pars->cos_theta;
  int lmax = kernel_pars->lmax;
  double filterpar = kernel_pars->filter;

  // EOS specific parameters
  double eta = eos_pars->mu_e / eos_pars->temp;
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

  MyKernelQuantity pair_kernel =
      {.abs_e = pair_kernel_absorption_e, .em_e = pair_kernel_production_e, .abs_x = pair_kernel_absorption_x, .em_x = pair_kernel_production_x};

  return pair_kernel;

}

/* Calculates the production and absorption kernels for the pair process for M1 (l = 0)
 *
 * R^p_{TP}(omega,omega',cos theta) = (1/2) Phi_l(omega,omega')
 * R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T} R^p_{TP}(omega,omega',cos theta)
 *
 */
MyKernelQuantity PairKernelsM1(MyEOSParams *eos_pars, PairKernelParams *kernel_pars) {
  // EOS specific parameters
  double eta = eos_pars->mu_e / eos_pars->temp;
  double temp = eos_pars->temp;

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;

  double pair_kernel_production_e = 0.;
  double pair_kernel_absorption_e = 0.;
  double pair_kernel_production_x = 0.;
  double pair_kernel_absorption_x = 0.;

  double pair_phi_e = 0.;
  double pair_phi_x = 0.;

  pair_phi_e = PairPhi(0, omega, omega_prime, eta, temp, 0);
  pair_phi_x = PairPhi(0, omega, omega_prime, eta, temp, 1);

  pair_kernel_production_e = 0.5 * pair_phi_e;
  pair_kernel_production_x = 0.5 * pair_phi_x;

  pair_kernel_absorption_e = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_e;
  pair_kernel_absorption_x = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_x;

  MyKernelQuantity pair_kernel = {.abs_e = pair_kernel_absorption_e, .em_e = pair_kernel_production_e,
      .abs_x = pair_kernel_absorption_x, .em_x = pair_kernel_production_x};

  return pair_kernel;

}

MyKernelQuantity PairKernelsM1Optimized(MyEOSParams *eos_pars, PairKernelParams *kernel_pars) {
  // EOS specific parameters
  double eta = eos_pars->mu_e / eos_pars->temp;
  double temp = eos_pars->temp;

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;

  double pair_kernel_production_e = 0.;
  double pair_kernel_absorption_e = 0.;
  double pair_kernel_production_x = 0.;
  double pair_kernel_absorption_x = 0.;
  double pair_phi[2] = {0.};

  PairPhiOptimized(kernel_pars, 0, eta, temp, pair_phi);

  pair_kernel_production_e = 0.5 * pair_phi[0];
  pair_kernel_production_x = 0.5 * pair_phi[1];

  pair_kernel_absorption_e = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_e;
  pair_kernel_absorption_x = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_x;

  MyKernelQuantity pair_kernel = {.abs_e = pair_kernel_absorption_e, .em_e = pair_kernel_production_e,
      .abs_x = pair_kernel_absorption_x, .em_x = pair_kernel_production_x};

  return pair_kernel;

}

/* Calculates the production and absorption kernels for the pair process from Eqns. (2) and (3) of Pons et. al.
 * integrated over the phi variable
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega') P_l(mu) P_l(mu_prime) 2 pi
 * R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T} R^p_{TP}
 *
 * l is an integer.
 */
MyKernelQuantity PairKernelsPhiIntegrated(MyEOSParams *eos_pars, PairKernelParams *kernel_pars) {

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;
  double mu = kernel_pars->mu;
  double mu_prime = kernel_pars->mu_prime;
  int lmax = kernel_pars->lmax;
  double filterpar = kernel_pars->filter;

  // EOS specific parameters
  double eta = eos_pars->mu_e / eos_pars->temp;
  double temp = eos_pars->temp;

  double pair_kernel_production_e = 0.;
  double pair_kernel_absorption_e = 0.;
  double pair_kernel_production_x = 0.;
  double pair_kernel_absorption_x = 0.;
  double pair_phi_e = 0.;
  double pair_phi_x = 0.;

  assert(lmax >= 0 && lmax <= 3);

  for (int l = 0; l <= lmax; l++) {
    double legendre_l_mu = gsl_sf_legendre_Pl(l, mu);
    double legendre_l_mu_prime = gsl_sf_legendre_Pl(l, mu_prime);
    pair_phi_e = PairPhi(l, omega, omega_prime, eta, temp, 0);
    pair_phi_x = PairPhi(l, omega, omega_prime, eta, temp, 1);

    pair_kernel_production_e += 2. * kPi * (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) * (2. * l + 1.) * pair_phi_e * legendre_l_mu * legendre_l_mu_prime / 2.;
    pair_kernel_production_x += 2. * kPi * (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) * (2. * l + 1.) * pair_phi_x * legendre_l_mu * legendre_l_mu_prime / 2.;
  }

  pair_kernel_absorption_e = exp((omega + omega_prime) / temp) * pair_kernel_production_e;
  pair_kernel_absorption_x = exp((omega + omega_prime) / temp) * pair_kernel_production_x;

  MyKernelQuantity pair_kernel =
      {.abs_e = pair_kernel_absorption_e, .em_e = pair_kernel_production_e, .abs_x = pair_kernel_absorption_x, .em_x = pair_kernel_production_x};

  return pair_kernel;

}

MyKernelQuantity PairKernelsPhiMuIntegrated(MyEOSParams *eos_pars, PairKernelParams *kernel_pars) {

  // kernel specific parameters
  double omega = kernel_pars->omega;
  double omega_prime = kernel_pars->omega_prime;

  // EOS specific parameters
  double eta = eos_pars->mu_e / eos_pars->temp;
  double temp = eos_pars->temp;

  int l = 0;

  double pair_phi_e = PairPhi(l, omega, omega_prime, eta, temp, 0);
  double pair_phi_x = PairPhi(l, omega, omega_prime, eta, temp, 1);

  double pair_kernel_production_e = 2. * kPi * pair_phi_e;
  double pair_kernel_production_x = 2. * kPi * pair_phi_x;

  double pair_kernel_absorption_e = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_e;
  double pair_kernel_absorption_x = SafeExp((omega + omega_prime) / temp) * pair_kernel_production_x;

  const double kPrefactor = 1. / (kClight * pow(kH * kClight, 3.));
  //const double kPrefactor = 1.;

  MyKernelQuantity pair_kernel = {.abs_e = kPrefactor * pair_kernel_absorption_e, .em_e = kPrefactor * pair_kernel_production_e,
      .abs_x = kPrefactor * pair_kernel_absorption_x, .em_x = kPrefactor * pair_kernel_production_x};

  return pair_kernel;
}
