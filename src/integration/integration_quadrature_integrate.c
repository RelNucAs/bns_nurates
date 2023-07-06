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
 *    nx:        number of points in the quadrature scheme
 *    wtarray:  array of quadrature weights
 *    fnarray:  array of function values at corresponding quadrature positions
 */
double DoIntegration(const int n, const double *wtarray, const double *fnarray) {
  double integral = 0.;

  for (int i = 0; i < n; i++) {
    integral += wtarray[i] * fnarray[i];
  }

  return integral;
}

/* Integrate a function from 0 to inf using a Gauss-Legendre quadrature by
 * breaking the integral into two parts, from [0,t] and [t,inf]. The other dimensions (if applicable)
 * have finite limits.
 *
 * Note: Uses (4.4.2) of Numerical Recipes 3rd edition (Press et. al.) to compute the integral from [t,inf]
 *       The quad struct should contain a Gauss-Legendre quadrature in [0,1] only!
 *
 * Inputs:
 *    quad: A properly populated Gauss-Legendre quadrature struct from x1 = 0 to x2 = 1.
 *    func: The function struct to be integrated
 *    t:    The value of points at which to break the integral into two
 */
double GaussLegendreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func, double t) {

  double f1_x[quad->nx], f2_x[quad->nx], f_y[quad->ny], f_z[quad->nz];
  double w_y[quad->ny], w_z[quad->nz];
  //var3d var = var3d_default;
  double var[3];

  for (int k = 0; k < quad->nz; k++) {
    for (int j = 0; j < quad->ny; j++) {
      for (int i = 0; i < quad->nx; i++) {
        var[0] = t * quad->points[i];
        var[1] = quad->points[quad->nx + j];
        var[2] = quad->points[quad->nx + quad->ny + k];

        f1_x[i] = func->function(var, func->params);

        var[0] = t / quad->points[i];

        f2_x[i] = func->function(var, func->params) / (quad->points[i] * quad->points[i]);
      }
      f_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_x) + DoIntegration(quad->nx, quad->w, f2_x));
      w_y[j] = quad->w[quad->nx + j];
    }
    f_z[k] = DoIntegration(quad->ny, w_y, f_y);
    w_z[k] = quad->w[quad->nx + quad->ny + k];
  }

  return DoIntegration(quad->nz, f_z, w_z);
}

/* Integrate a function from 0 to inf using a Gauss-Laguerre quadrature
 *
 * Inputs:
 *    quad: A properly populated Gauss-Laguerre quadrature struct
 *    func:    The function struct to be integrated
 */
double GaussLaguerreIntegrateZeroInf(MyQuadrature *quad, MyFunction *func) {

  double f[quad->nx];

  if (quad->alpha == 0.) {
    for (int i = 0; i < quad->nx; i++) {
      f[i] = func->function(&quad->points[i], func->params); // * exp(quad->points[i]);
    }
  } else {
    for (int i = 0; i < quad->nx; i++) f[i] = func->function(&quad->points[i], func->params); // * exp(quad->points[i]) / pow(quad->points[i], quad->alpha);
  }

  return DoIntegration(quad->nx, quad->w, f);

}
