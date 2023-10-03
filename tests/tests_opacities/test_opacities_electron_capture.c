//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_electron_capture.c
//  \brief Generate a table for electron positron capture

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "../../src/opacities/opacities.h"
#include "../../src/integration/integration.h"

int main() {

  printf("=================================================== \n");
  printf("Testing opacities for electron/positron capture ... \n");
  printf("=================================================== \n");

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_opacities/nurates_CCSN/nurates_1.008E+01.txt";

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data_directory: %s\n", filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE_NAME__, filepath);
    exit(1);
  }

  // store columns here
  int num_data = 102;
  double e_nu;
  int zone[num_data];
  double r[num_data];
  double rho[num_data];
  double T[num_data];
  double Ye[num_data];
  double mu_e[num_data];
  double mu_hat[num_data];
  double Yh[num_data];
  double Ya[num_data];
  double Yp[num_data];
  double Yn[num_data];
  double em_nue[num_data];
  double l_nue_inv[num_data];
  double em_anue[num_data];
  double l_anue_inv[num_data];

  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[1] == '#' && i == 0) {
      sscanf(line + 14, "%lf\n", &e_nu);
      continue;
    } else if (line[1] == '#' && i != 0) {
      continue;
    }

    sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
             &zone[i], &r[i], &rho[i], &T[i], &Ye[i], &mu_e[i], &mu_hat[i], &Yh[i], &Ya[i], &Yp[i], &Yn[i], &em_nue[i], &l_nue_inv[i], &em_anue[i], &l_anue_inv[i]);
    i++;
  }

  printf("\n");
  printf("Printing out %s\n", outname);
  printf("\n");

  printf("E_nu: %lf\n", e_nu);
  for (int i = 0; i < num_data; i++) {
      printf("%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             zone[i], r[i], rho[i], T[i], Ye[i], mu_e[i], mu_hat[i], Yh[i], Ya[i], Yp[i], Yn[i], em_nue[i], l_nue_inv[i], em_anue[i], l_anue_inv[i]);
    }

  printf("\n");
  printf("Generating quadratures ...\n");
  MyQuadrature my_quadrature_1d = {.nx = 60, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendre(&my_quadrature_1d);
  MyQuadrature my_quadrature_2d = {.nx = 60, .ny = 60, .dim = 2, .type = kGauleg, .x1 = 0., .x2 = 1., .y1 = 0., .y2 = 1.};
  GaussLegendre(&my_quadrature_2d);

  GreyOpacityParams my_grey_opacity_params;

  for (int i = 0; i < num_data; i++) {

    // set EOS params from table
    my_grey_opacity_params.eos_pars.mu_e = mu_e[i];
    my_grey_opacity_params.eos_pars.temp = T[i];
    my_grey_opacity_params.eos_pars.yp = Yp[i];
    my_grey_opacity_params.eos_pars.yn = Yn[i];

    my_grey_opacity_params.opacity_pars.use_WM_ab = false;
    my_grey_opacity_params.opacity_pars.use_WM_sc = false;
    my_grey_opacity_params.opacity_pars.use_dU = false;

    // only enable Fermi distribution function for the three neutrino species
    my_grey_opacity_params.distr_pars.temp_t[0] = T[i];
    my_grey_opacity_params.distr_pars.temp_t[1] = T[i];
    my_grey_opacity_params.distr_pars.temp_t[2] = T[i];
    my_grey_opacity_params.distr_pars.w_t[0] = 1.;
    my_grey_opacity_params.distr_pars.w_t[1] = 1.;
    my_grey_opacity_params.distr_pars.w_t[2] = 1.;
    my_grey_opacity_params.distr_pars.eta_t[0] = mu_e[i] / T[i];
    my_grey_opacity_params.distr_pars.eta_t[1] = - mu_e[i] / T[i];
    my_grey_opacity_params.distr_pars.eta_t[2] = 0.;

    // disable thin distribution function for the time being
    my_grey_opacity_params.distr_pars.w_f[0] = 0.;
    my_grey_opacity_params.distr_pars.w_f[1] = 0.;
    my_grey_opacity_params.distr_pars.w_f[2] = 0.;
    my_grey_opacity_params.distr_pars.c_f[0] = 1.;
    my_grey_opacity_params.distr_pars.c_f[1] = 1.;
    my_grey_opacity_params.distr_pars.c_f[2] = 1.;
    my_grey_opacity_params.distr_pars.temp_f[0] = T[i];
    my_grey_opacity_params.distr_pars.temp_f[1] = T[i];
    my_grey_opacity_params.distr_pars.temp_f[2] = T[i];

    my_grey_opacity_params.opacity_flags = opacity_flags_default_none;
    my_grey_opacity_params.opacity_flags.use_abs_em = 1;

  }
  return 0;
}

