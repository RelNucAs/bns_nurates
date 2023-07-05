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

/* Generate and save quadratures for integration. Saves files to the
 * filedir and loads the quadrature data to a quad struct provided that
 * the quad struct has properly filled metadata
 *
 * Inputs:
 *    filedir:  location of save directory
 *    quad:     quadrature struct with proper metadata populated
 */
void SaveQuadrature(char *filedir, MyQuadrature *quad) {

  char fileHeader[100];
  char fileline[50];
  char outname[50];
  char filepath[200] = {'\0'};

  FILE *fptr;

  // To add a nre quadrature, add another if/else case here
  if (quad->type == kGauleg && quad->dim == 1) {
    sprintf(fileHeader, "# Abscissas and weights for Gauss-Legendre integration from x1 = %.3lf to x2 = %.3lf\n", quad->x1, quad->x2);
    GaussLegendre(quad);
    sprintf(outname, "/GaussLegQuad_n_%d_x_%f_%f.txt", quad->nx, quad->x1, quad->x2);
  } else if (quad->type == kGaulag && quad->dim == 1) {
    sprintf(fileHeader, "# Abscissas and weights for Gauss-Laguerre integration\n");
    GaussLaguerre(quad);
    sprintf(outname, "/GaussLagQuad_n_%d_alf_%f.txt", quad->nx, quad->alpha);
  } else {
    printf("The requested quadrature is not implemented!\n");
    exit(1);
  }

  strcat(filepath, filedir);
  strcat(filepath, outname);

  fptr = fopen(filepath, "w");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE_NAME__, filepath);
    exit(1);
  }

  fprintf(fptr, fileHeader);

  // @TODO: only works for the 1d case, generalize for 2d case later
  for (int i = 0; i < quad->nx; i++) {
    sprintf(fileline, "%0.16e \t %0.16e\n", quad->points[i], quad->w[i]);
    fprintf(fptr, fileline);
  }

  fclose(fptr);

}

/* Load quadrature from disk to struct quad. If not present on disk,
 * generate and save to disk as well.
 *
 * Inputs:
 *    filedir:  location of save directory
 *    quad:     quadrature struct with proper metadata populated
 */
void LoadQuadrature(char *filedir, MyQuadrature *quad) {

  char fileHeader[100];
  char outname[50];
  char filepath[200] = {'\0'};

  FILE *fptr;

  strcat(filepath, filedir);

  if (quad->type == kGauleg && quad->dim == 1) {
    sprintf(fileHeader, "# Abscissas and weights for Gauss-Legendre integration from x1 = %.3lf to x2 = %.3lf\n", quad->x1, quad->x2);
    sprintf(outname, "/GaussLegQuad_n_%d_x_%f_%f.txt", quad->nx, quad->x1, quad->x2);
  } else if (quad->type == kGaulag && quad->dim == 1) {
    sprintf(fileHeader, "# Abscissas and weights for Gauss-Laguerre integration\n");
    sprintf(outname, "/GaussLagQuad_n_%d_alf_%f.txt", quad->nx, quad->alpha);
  } else {
    printf("The requested quadrature is not implemented!\n");
    exit(1);
  }

  strcat(filepath, outname);

  fptr = fopen(filepath, "r");

  if (fptr == NULL) {
    printf("The file %s does not exist! Creating ...\n", filepath);
    SaveQuadrature(filedir, quad);
  } else {
    quad->points = (double *) realloc(quad->points, quad->nx * sizeof(double));
    quad->w = (double *) realloc(quad->w, quad->nx * sizeof(double));

    fscanf(fptr, fileHeader);

    // @TODO: currently works only for 1d quadratures, generalize for 2d later
    for (int i = 0; i < quad->nx; i++) {
      fscanf(fptr, "%.16e %.16e\n", quad->points[i], quad->w[i]);
    }
  }

  if (fptr != NULL) {
    fclose(fptr);
  }

}