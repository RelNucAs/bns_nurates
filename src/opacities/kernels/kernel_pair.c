//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.c
//  \brief contains pair kernels and associated helper functions

#include <math.h>
#include <assert.h>
#include "kernels.h"
//#include "../../constants.h"
#include "../../bns_nurates.h"

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
double T(int l, double alpha, double tolerance) {

  assert(alpha >= 0 && l >= 1);

  if (alpha == 0 && l != 1) {
    //return pow(2., -l) * (pow(2., l) - 2.) * boost::math::zeta<double>(l);  // Computed in Mathematica
    // @TODO: fix
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