//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_integrate.c
//  \brief routines for integrating functions

#include <math.h>
#include "integration.h"

/* Make a convolution between sampled function and weights
 *
 * Inputs:
 *    n:        number of points in the quadrature scheme
 *    wgt:      array of quadrature weights
 *    fnarray:  array of function values at corresponding quadrature positions
 */
double DoIntegration(const int n, const double *wgt, const double *fnarray) {
  double integral = 0.;

  for (int i = 0; i < n; i++) {
    integral += wgt[i] * fnarray[i];
  }

  return integral;
}

/* Integrate a function from 0 to inf using a Gauss-Legendre quadrature by
 * breaking the integral into two parts, from [0,t] and [t,inf].
 *
 * Note: Uses (4.4.2) of Numerical Recipes 3rd edition (Press et. al.) to compute the integral from [t,inf]
 *       The quad struct should contain a Gauss-Legendre quadrature in [0,1] only!
 *
 * Inputs:
 *    quad: A properly populated Gauss-Legendre quadrature struct from x1 = 0 to x2 = 1.
 *    func:    The function struct to be integrated
 *    t:    The value of x at which to break the integral into two
 */
double GaussLegendreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func, double t) {

  double f1[quad->n], f2[quad->n];

  for (int i = 0; i < quad->n; i++) {
    f1[i] = func->function(t * quad->x[i], func->params);
    f2[i] = func->function(t / quad->x[i], func->params) / (quad->x[i] * quad->x[i]);
  }

  return t * (DoIntegration(quad->n, quad->w, f1) + DoIntegration(quad->n, quad->w, f2));
}

/* Integrate a function from 0 to inf using a Gauss-Laguerre quadrature
 *
 * Inputs:
 *    quad: A properly populated Gauss-Laguerre quadrature struct
 *    func:    The function struct to be integrated
 */
double GaussLaguerreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func) {

  double f[quad->n];

  if (quad->alpha == 0.) {
    for (int i = 0; i < quad->n; i++) {
      f[i] = func->function(quad->x[i], func->params) * exp(quad->x[i]);
    }
  } else {
    for (int i = 0; i < quad->n; i++) f[i] = func->function(quad->x[i], func->params) * exp(quad->x[i]) / pow(quad->x[i], quad->alpha);
  }

  return DoIntegration(quad->n, quad->x, f);

}