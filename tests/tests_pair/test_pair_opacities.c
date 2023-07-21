//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_opacities.c
//  \brief Generate emissivity and absorptivity tables considering thermal equilibrium and compare with Albino's result

#include <stdio.h>
#include "../../src/bns_nurates.h"
#include "../../src/integration/integration.h"
#include "../../src/opacities/opacities.h"

int main() {

  int max_error = -42.;

  // Data from Albino
  double temp = 30.;
  double mu_e = 4. * 30;

  double data_albino[12][3] = { //omega(MeV) emissivity (km^-1) mean free path (km) [x neutrinos]
      {2.0000000000000000, 1208.8403883182957, 235.74209542285573},
      {3.0398221659058673, 1788.6544282931936, 153.90515091871933},
      {4.6202594001663186, 2610.0591334007413, 100.06698725757252},
      {7.0223834684302613, 3729.4256545265321, 64.652907643784999},
      {10.673398462412617, 5162.2821607198175, 41.365021291439454},
      {16.222616615793736, 6812.1178410364482, 26.062527399996728},
      {24.656934788841312, 8368.0149227064740, 16.026205440085686},
      {37.476348457207664, 9238.7406326351138, 9.4772012428358074},
      {56.960717368716004, 8698.7685847840530, 5.2659134502558507},
      {86.575225621661119, 6403.0283673993918, 2.6727774748811468},
      {131.58664493151352, 3059.4607054080661, 1.2518405737406311},
      {199.99999999999986, 637.21742237179706, 0.61609418279642836}
  };


  int lmax = 0;
  int filter = 0.;
  MyQuadrature quad = quadrature_default;
  quad.dim = 1;
  quad.type = kGauleg;
  quad.x1 = 0.;
  quad.x2 = 1.;
  quad.nx = 45;

  GaussLegendreMultiD(&quad);

  printf("Quadrature information:\n");
  printf("points \t\t weights\n");
  for (int i = 0; i < (quad.nx + quad.ny + quad.nz); i++) {
    if (i == quad.nx || i == quad.nx + quad.ny) {
      printf("\n");
    }
    printf("%f \t %f\n", quad.points[i], quad.w[i]);
  }
  printf("\n");

  MyKernelParams my_kernel_params;
  my_kernel_params.pair_kernel_params.filter = filter;
  my_kernel_params.pair_kernel_params.lmax = lmax;
  my_kernel_params.pair_kernel_params.mu = 1.;
  my_kernel_params.pair_kernel_params.mu_prime = 1.;
  my_kernel_params.pair_kernel_params.cos_theta = 1.;
  my_kernel_params.pair_kernel_params.omega = 0.;
  my_kernel_params.pair_kernel_params.omega_prime = 0.;

  printf("Pair kernel parameters:\n");
  printf("Filter parameter = %f\n", my_kernel_params.pair_kernel_params.filter);
  printf("lmax = %f\n", my_kernel_params.pair_kernel_params.lmax);

  MyEOSParams my_eos_params;
  my_eos_params.mu_e = mu_e;
  my_eos_params.temp = temp;

  for (int i = 0; i < 12; i++) {
    my_kernel_params.pair_kernel_params.omega = data_albino[0][i];
    MyOpacityQuantity result = PairOpacitiesFermi(&quad, &my_eos_params, &my_kernel_params);

    printf("%0.16e %.16f %0.16f\n", data_albino[i][0], result.em_x, 1. / result.abs_x);
  }

  return 1;
}