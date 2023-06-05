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

void gauss_legendre(MyQuadrature *quad) {

  assert(quad->type == kGaulag);
  assert(quad->dim == 1);

  double *x = quad->x;
  double *w = quad->w;

  x = (double *) malloc(quad->n * sizeof(double));
  w = (double *) malloc(quad->n * sizeof(double));

  const double EPS = 1.0e-10; //1.0e-14;
  double z1, z, xm, xl, pp, p3, p2, p1;

  int n = quad->n;
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

    } while (fabs(z - z1) > EPS);

    x[i] = xm - xl * z;
    x[n - 1 - i] = xm + xl * z;
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];

  }
}