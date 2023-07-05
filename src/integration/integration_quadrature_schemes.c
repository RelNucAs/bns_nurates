//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_generate.c
//  \brief generate quadrature for integration

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "../bns_nurates.h"
#include "integration.h"
#include "../functions/functions.h"

/* Generate Gauss-Legendre quadratures in [x1,x2].
 * For this routine to generate data successfully, quad struct
 * must have dim, type, nx, x1, x2 metadata populated.
 *
 * Inputs:
 *      quad: A MyQuadrature structure to hold quadratures.
 *            This already contains metadata for the quadrature, the routine
 *            only populates the quadrature points and weights
 */
void GaussLegendre(MyQuadrature *quad) {

  const double kEps = 1.0e-10; //1.0e-14;

  assert(quad->dim == 1);
  assert(quad->type == kGauleg);

  quad->points = NULL;
  quad->w = NULL;
  quad->points = (double *) malloc(quad->nx * sizeof(double));
  quad->w = (double *) malloc(quad->nx * sizeof(double));

  double z1, z, xm, xl, pp, p3, p2, p1;

  int n = quad->nx;
  int m = (n + 1) / 2;
  xm = 0.5 * (quad->x2 + quad->x1);
  xl = 0.5 * (quad->x2 - quad->x1);

  for (int i = 0; i < m; i++) {
    z = cos(3.141592654 * (i + 0.75) / (n + 0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
      }

      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;

    } while (fabs(z - z1) > kEps);

    quad->points[i] = xm - xl * z;
    quad->points[n - 1 - i] = xm + xl * z;
    quad->w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    quad->w[n - 1 - i] = quad->w[i];

  }

}

/* Generate Gauss-Laguerre quadratures in [0,inf).
 * For this routine to generate data successfully, quad struct
 * must have dim, type, nx, alpha metadata populated.
 *
 * Inputs:
 *      quad:   A MyQuadrature structure to hold quadratures.
 *              This already contains metadata for the quadrature, the routine
 *              only populates the quadrature points and weights
 */
void GaussLaguerre(MyQuadrature *quad) {

  const int kMaxIt = 10;
  const double kEps = 1.0e-14;

  assert(quad->dim == 1);
  assert(quad->type == kGaulag);

  quad->points = NULL;
  quad->w = NULL;
  quad->points = (double *) malloc(quad->nx * sizeof(double));
  quad->w = (double *) malloc(quad->nx * sizeof(double));

  double *x = quad->points;
  double *w = quad->w;

  int i, its, j;
  double ai, p1, p2, p3, pp, z, z1;
  int n = quad->nx;

  for (i = 0; i < n; i++) {
    if (i == 0) {
      z = (1.0 + quad->alpha) * (3.0 + 0.92 * quad->alpha) / (1.0 + 2.4 * n + 1.8 * quad->alpha);
    } else if (i == 1) {
      z += (15.0 + 6.25 * quad->alpha) / (1.0 + 0.9 * quad->alpha + 2.5 * n);
    } else {
      ai = i - 1;
      z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * quad->alpha / (1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * quad->alpha);
    }
    for (its = 0; its < kMaxIt; its++) {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1 + quad->alpha - z) * p2 - (j + quad->alpha) * p3) / (j + 1);
      }
      pp = (n * p1 - (n + quad->alpha) * p2) / z;
      z1 = z;
      z = z1 - p1 / pp;
      if (fabs(z - z1) <= kEps) { break; }
    }
    if (its >= kMaxIt) {
      throw("Too many iterations for Gauss-Laguerre!")
    };

    x[i] = z;
    w[i] = -exp(Gammln(quad->alpha + n) - Gammln(((double) n))) / (pp * n * p2);

  }
}