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
 *
 * n:       number of points in the gradrature scheme
 * wgt:     array of quadrature weights
 * fnarray: array of function values at corresponding quadrature positions
*/
inline double DoIntegration(const int n, const double *wgt, const double *fnarray) {
  double integral = 0.;

  for (int i = 0; i < n; i++) {
    integral += wgt[i] * fnarray[i];
  }

  return integral;
}