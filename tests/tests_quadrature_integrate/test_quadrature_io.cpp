//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_quadrature_io.c
//  \brief Test quadrature read/write features

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

int main() {

  double max_error = -42.;

  printf("========================== \n");
  printf("Testing quadrature I/O ... \n");
  printf("========================== \n");

  char filedir[300];
  if (getcwd(filedir, sizeof(filedir)) == NULL) {
    perror("Failed to get save directory for quadrature!");
  } else {
    printf("Quadrature save directory: %s\n", filedir);
  }

  // cheking quadrature creating and saving to disk (Gauss-Legendre)
  MyQuadrature quad = {.dim = 1, .type=kGauleg, .nx = 9, .x1 = -1., .x2 = 1.};
  MyQuadrature load_quad = {.dim = 1, .type=kGauleg, .nx = 9, .x1 = -1., .x2 = 1.};
  printf("Saving quadrature information for Gauss-Legendre...\n");
  SaveQuadrature(filedir, &quad);
  printf("Loading saved quadrature ...\n");
  LoadQuadrature(filedir, &load_quad);
  printf("\n");

  printf("nx = %d \t x1 = %0.16e \t x2 = %0.16e\n", quad.nx, quad.x1, quad.x2);
  for (int i = 0; i < quad.nx; i++) {
    printf("%0.16e \t %0.16e \t %0.16e \t %0.16e\n", quad.points[i], load_quad.points[i], quad.w[i], load_quad.w[i]);
    if(fabs(quad.points[i]-load_quad.points[i]) > max_error) {
      max_error = fabs(quad.points[i]-load_quad.points[i]);
    }
    if(fabs(quad.w[i]-load_quad.w[i]) > max_error) {
      max_error = fabs(quad.w[i]-load_quad.w[i]);
    }
  }

  printf("\n");
  printf("Max error: %0.16e\n", max_error);

  if(max_error < 1e-12) {
    return 0;
  } else {
    return 1;
  }
}