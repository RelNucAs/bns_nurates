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
 *
 * Note: This routine can only generate quadrature in 1d.
 *
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

/* Generate quadratures in multiple dimensions
 *
 * Note: For almost all practical cases, this routine should be called and not the 1d routine!
 *
 * Generates 1d/2d/3d quadratures. Provide a quadrature structure with properly populated metadata.
 * This routine populates the quadrature points & weights arrays with (nx + ny + nz) elements each.
 *
 * The first nx elements contain quadrature information for the x variable, the next ny for the y variable
 * and so on.
 *
 * For the 1d case, the quadrature arrays are of the length (nx + 2) each, with the weight and point information
 * in the y and z direction being populated with 1.
 *
 * Similarly, for the 2d case, the quadrature arrays are of length (nx + ny + 1) each.
 *
 * Inputs:
 *      quad: A MyQuadrature structure to hold quadratures.
 *            This already contains metadata for the quadrature, the routine
 *            only populates the quadrature points and weights
 */
void GaussLegendreMultiD(MyQuadrature *quad) {

  //perform checks
  assert(quad->type == kGauleg);
  assert(quad->dim == 1 || quad->dim == 2 || quad->dim == 3);

  // reallocate memory for the weights and points
  // Note: first nx elements contain quadrature for x variable, then ny elements for y variable and then nx for z
  quad->points = NULL;
  quad->w = NULL;
  quad->points = (double *) malloc(quad->nx + quad->ny + quad->nz * sizeof(double));
  quad->w = (double *) malloc(quad->nx + quad->ny + quad->nz * sizeof(double));

  // local quadrature which contain 1d quadrature along each variable
  MyQuadrature quad_local = quadrature_default;

  if (quad->dim == 1) {
    quad_local.nx = quad->nx;
    quad_local.x1 = quad->x1;
    quad_local.x2 = quad->x2;

    GaussLegendre(&quad_local);

    for (int i = 0; i < quad->nx; i++) {
      quad->points[i] = quad_local.points[i];
      quad->w[i] = quad_local.w[i];
    }
    quad->points[quad->nx] = 1.;
    quad->w[quad->nx] = 1.;
    quad->points[quad->nx + 1] = 1.;
    quad->w[quad->nx + 1] = 1.;

  }

  if (quad->dim == 2) {
    quad_local.nx = quad->ny;
    quad_local.x1 = quad->y1;
    quad_local.x2 = quad->y2;

    GaussLegendre(&quad_local);

    for (int i = 0; i < quad->ny; i++) {
      quad->points[quad->nx + i] = quad_local.points[i];
      quad->w[quad->nx + i] = quad_local.w[i];
    }
    quad->points[quad->nx + quad->ny] = 1.;
    quad->w[quad->nx + quad->ny] = 1.;
    quad->points[quad->nx + quad->ny + 1] = 1.;
    quad->w[quad->nx + quad->ny + 1] = 1.;
  }

  if (quad->dim == 3) {
    quad_local.nx = quad->nz;
    quad_local.x1 = quad->z1;
    quad_local.x2 = quad->z2;

    GaussLegendre(&quad_local);

    for (int i = 0; i < quad->nz; i++) {
      quad->points[quad->nx + quad->ny + i] = quad_local.points[i];
      quad->w[quad->nx + quad->ny + quad->ny + i] = quad_local.w[i];
    }
  }

  free(quad_local.points);
  free(quad_local.w);

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