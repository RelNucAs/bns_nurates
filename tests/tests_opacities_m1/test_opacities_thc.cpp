//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_all_reactions.c
//  \brief Generate a table for all reactions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "../tests.hpp"
#include "../include/bns_nurates.hpp"
#include "../include/distribution.hpp"
#include "../include/integration.hpp"
#include "../include/opacities.hpp"
#include "m1_opacities.hpp"

int main() {

  GreyOpacityParams my_grey_opacity_params = {0};

  MyQuadrature my_quadrature_1d;
  MyQuadrature my_quadrature_2d;

  int quad_nx;
  int quad_ny;

  bool use_abs_em;
  bool use_brem;
  bool use_pair;
  bool use_iso_scat;
  bool use_WM_ab;
  bool use_WM_sc;
  bool use_dU;

  // dummy values for now
  use_abs_em = true;
  use_brem = true;
  use_pair = true;
  use_iso_scat = true;
  use_WM_ab = false;
  use_WM_sc = false;
  use_dU = false;

  quad_nx = 60;
  quad_ny = quad_nx;

  int num_data = 4;
  double nb_array[num_data];
  double temp_array[num_data];
  double ye_array[num_data];
  double mu_n_array[num_data];
  double mu_p_array[num_data];
  double mu_e_array[num_data];

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_distribution/data/thc_test_data.txt";
  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data_directory: %s\n", filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  // read in the data file
  int i = 0;
  char line[1000];
  printf("\n");
  printf("#nb temp ye mu_n mu_p mu_e\n");
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[0] == '#') {
      continue;
    }

    sscanf(line, "%lf %lf %lf %lf %lf %lf\n",
           &nb_array[i], &temp_array[i], &ye_array[i], &mu_n_array[i], &mu_p_array[i], &mu_e_array[i]);
    printf("%e %lf %lf %lf %lf %lf\n", nb_array[i], temp_array[i], ye_array[i], mu_n_array[i], mu_p_array[i], mu_e_array[i]);

    i++;
  }
  printf("\n");
  // nb = 1.0512162546485660E+039;
  // temp = 55.800010783866689;
  // ye = 0.23432446891781622;
  // mu_n = 795.93707619905956;
  // mu_p = 439.37056823706115;
  // mu_e = 356.01566068372557;

  // compute quadratures
  my_quadrature_1d.nx = quad_nx;
  my_quadrature_1d.dim = 1;
  my_quadrature_1d.type = kGauleg;
  my_quadrature_1d.x1 = 0.;
  my_quadrature_1d.x2 = 1.;
  GaussLegendreMultiD(&my_quadrature_1d);

  my_quadrature_2d.nx = quad_nx;
  my_quadrature_2d.ny = quad_ny;
  my_quadrature_2d.dim = 2;
  my_quadrature_2d.type = kGauleg;
  my_quadrature_2d.x1 = 0.;
  my_quadrature_2d.x2 = 1.;
  my_quadrature_2d.y1 = 0.;
  my_quadrature_2d.y2 = 1.;
  GaussLegendreMultiD(&my_quadrature_2d);

  // reaction flags
  my_grey_opacity_params.opacity_flags = opacity_flags_default_none;
  my_grey_opacity_params.opacity_flags.use_abs_em = use_abs_em;
  my_grey_opacity_params.opacity_flags.use_brem = use_brem;
  my_grey_opacity_params.opacity_flags.use_pair = use_pair;
  my_grey_opacity_params.opacity_flags.use_iso = use_iso_scat;

  // other flags
  my_grey_opacity_params.opacity_pars.use_WM_ab = use_WM_ab;
  my_grey_opacity_params.opacity_pars.use_WM_sc = use_WM_sc;
  my_grey_opacity_params.opacity_pars.use_dU = use_dU;


  printf(
        "#R_nue R_nua R_nux Q_nue Q_nua Q_nux sigma_0_nue sigma_0_nua sigma_0_nux sigma_1_nue sigma_1_nua sigma_1_nux scat_0_nue scat_0_nua scat_0_nux scat_1_nue scat_1_nua scat_1_nux\n");
  for (int idx = 0; idx < num_data; idx++) {
    // populate EOS quantities
    my_grey_opacity_params.eos_pars.mu_e = mu_e_array[idx];
    my_grey_opacity_params.eos_pars.mu_p = mu_p_array[idx];
    my_grey_opacity_params.eos_pars.mu_n = mu_n_array[idx];
    my_grey_opacity_params.eos_pars.temp = temp_array[idx];
    my_grey_opacity_params.eos_pars.yp = ye_array[idx];
    my_grey_opacity_params.eos_pars.yn = 1 - ye_array[idx];
    my_grey_opacity_params.eos_pars.nb = nb_array[idx];

    // reconstruct distribution function
    my_grey_opacity_params.distr_pars = NuEquilibriumParams(&my_grey_opacity_params.eos_pars);

    // compute n andj
    ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars, &my_grey_opacity_params.distr_pars, &my_grey_opacity_params.m1_pars);

    // populate M1 quantities
    my_grey_opacity_params.m1_pars.chi[id_nue] = 1.;
    my_grey_opacity_params.m1_pars.chi[id_anue] = 1.;
    my_grey_opacity_params.m1_pars.chi[id_nux] = 1.;
    
    // compute opacities
    M1Opacities opacities = ComputeM1Opacities(&my_quadrature_1d, &my_quadrature_2d, &my_grey_opacity_params);

    // extract emissivities
    double R_nue = opacities.eta_0[id_nue];
    double R_nua = opacities.eta_0[id_anue];
    double R_nux = 4.0 * opacities.eta_0[id_nux];
    double Q_nue = opacities.eta[id_nue];
    double Q_nua = opacities.eta[id_anue];
    double Q_nux = 4.0 * opacities.eta[id_nux];

    // extract absorption inverse mean-free path
    double sigma_0_nue = opacities.kappa_0_a[id_nue];
    double sigma_0_nua = opacities.kappa_0_a[id_anue];
    double sigma_0_nux = opacities.kappa_0_a[id_nux];
    double sigma_1_nue = opacities.kappa_a[id_nue];
    double sigma_1_nua = opacities.kappa_a[id_anue];
    double sigma_1_nux = opacities.kappa_a[id_nux];

    // extract scattering inverse mean-free path
    double scat_0_nue = 0;
    double scat_0_nua = 0;
    double scat_0_nux = 0;
    double scat_1_nue = opacities.kappa_s[id_nue];
    double scat_1_nua = opacities.kappa_s[id_anue];
    double scat_1_nux = opacities.kappa_s[id_nux];


    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", R_nue, R_nua, R_nux, Q_nue, Q_nua, Q_nux, sigma_0_nue, sigma_0_nua, sigma_0_nux,
           sigma_1_nue, sigma_1_nua, sigma_1_nux, scat_0_nue, scat_0_nua, scat_0_nux, scat_1_nue, scat_1_nua, scat_1_nux);
  }
  return 0;
}