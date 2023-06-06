//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_io.c
//  \brief write and read quadratures for integration

#include <stdio.h>
#include <string.h>
#include "../bns_nurates.h"
#include "integration.h"
#include "../functions/functions.h"

// Generate and save quadratures

/*
 * Inputs:
 *
 * n:   number of quadrature points
 * dim: dimension of quadrature (1d/2d)
 * type:  type of quadrature
 * x1:    lower limit of x
 * x2:    upper limit of x
 * y1:    lower limit of y (optional)
 * y2:    upper limit of y (optional)
 * alpha: Gauss-Laguerre parameter (optional)
 */
void SaveQuadrature(const int n, const int dim, enum Quadrature type, const double x1, const double x2, const double y1, const double y2, const double alpha) {

  char fileHeader[100];
  char fileline[50];

  MyQuadrature quad;
  quad.type = type;
  quad.dim = dim;
  quad.n = n;
  quad.x1 = x1;
  quad.x2 = x2;
  quad.y1 = y1;
  quad.y2 = y2;

  char outname[50];
  char *abs_path = "/var/home/maitraya/Documents/bns_nurates/src/";

  FILE *fptr;

  if (type == kGauleg && dim == 1) {
    sprintf(fileHeader, "# Abscissas and weights for Gauss-Legendre integration from x1 = %.3lf to x2 = %.3lf\n", x1, x2);
    GaussLegendre(&quad);
    sprintf(outname, "../quadratures/gausslegquad_n_%d.txt", n);
  } else if (type == kGaulag && dim == 1) {
    GaussLaguerre(quad, alpha);
    sprintf(outname, "../quadratures/gausslagquad_n_%d_alf_%d.txt", n, alpha);
  } else {
    printf("The requested quadrature is not implement!\n");
    exit(1);
  }

  fptr = fopen(strcat(abs_path, outname), "w");
  if (fptr == NULL) {
    printf("The directory %s does not exist!\n", abs_path);
    exit(1);
  }

  fprintf(fptr, fileHeader);

  for (int i = 0; i < n; i++) {
    sprintf(fileline, "%d \t %d\n", quad.x[i], quad.w[i]);
    fprintf(fptr, fileline);
  }

  fclose(fptr);

  free(quad.x);
  free(quad.y);
  free(quad.w);

}

void LoadQuadrature(MyQuadrature quad, const int n, const int dim, enum Quadrature type, const double x1, const double x2, const double y1, const double y2, const double alpha) {
  char fileHeader[100];
  char fileline[50];
  char outname[50];
  char *abs_path = "/var/home/maitraya/Documents/bns_nurates/src/";

  FILE *fptr;

  fptr = fopen(strcat(abs_path, outname), "r");

  if (fptr == NULL) {
    printf("The file %s does not exist!\n", strcat(abs_path, outname));
    exit(1);
  }

  quad.type = type;
  quad.dim = dim;
  quad.n = n;
  quad.x1 = x1;
  quad.x2 = x2;
  quad.y1 = y1;
  quad.y2 = y2;

  quad.x = (double *) malloc(quad.n * sizeof(double));
  quad.w = (double *) malloc(quad.n * sizeof(double));

  fscanf(fptr, fileHeader);

  for (int i = 0; i < n; i++) {
    fscanf(fptr, "%d %d\n", quad.x[i], quad.w[i]);
  }

  fclose(fptr);
}