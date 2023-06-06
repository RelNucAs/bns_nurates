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

double GaussLaguerreIntegrate(const int n, MyFunction F, double alpha) {
  MyQuadrature quad;
  double f[n];

  LoadQuadrature(quad, n, 1, kGaulag, 0, 1, -42, -42, alpha);

  if (alpha == 0.) {
    for (int i = 0; i < n; i++) f[i] = F.function(quad.x[i], F.params) * exp(quad.x[i]);
  } else {
    for (int i = 0; i < n; i++) f[i] = F.function(quad.x[i], F.params) * exp(quad.x[i]) / pow(quad.x[i], alpha);
  }

  return DoIntegration(n, quad.x, f);
}

double GaussLegendreIntegrate(const int n, MyFunction F, double t) {
  MyQuadrature quad;
  double f[n];
  double f1[n], f2[n];

  LoadQuadrature(quad, n, 1, kGaulag, 0, 1, -42, -42, -42.);

  for (int i = 0; i < n; i++) {
    f1[i] = F.function(t * quad.x[i], F.params);
    f2[i] = F.function(t / quad.x[i], F.params) / (quad.x[i] * quad.x[i]);
  }

  return t * (DoIntegration(n, quad.w, f1) + DoIntegration(n, quad.w, f2));
}