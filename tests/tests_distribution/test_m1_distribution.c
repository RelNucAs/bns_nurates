//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_bremsstrahlung.c
//  \brief Generate a table for bremsstrahlung

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "opacities.h"
#include "../../include/integration.h"
#include "../../include/distribution.h"
#include "../../include/constants.h"
#include "../../include/functions.h"

int main() {

  printf("=================================================== \n");
  printf("Testing distribution function reconstruction for M1 \n");
  printf("=================================================== \n");

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_distribution/data/distribution_function_m1_data.txt"; // file containing table generated from Boltzran

  char data_fileout[200] = "/tests/tests_distribution/data/distribution_function_params_from_m1.txt"; // file where data is written
  char data_filepath[300] = {'\0'};

  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data file path: %s\n", data_filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  double energy_vals[] = {3.36092278, 4.28273982, 5.45738823, 6.9542133, 8.86158006, 11.2920898, 14.38922757, 18.33583276,
                          23.36489302, 29.77329872, 37.9393698, 48.34518991, 61.60506618, 78.50179483, 100.03287349, 127.46938843,
                          162.43105312, 206.98182789, 263.75176578, 336.09227757};
  int num_data = 102;
  double T[num_data];
  double mu_e[num_data];
  double mu_p[num_data];
  double mu_n[num_data];
  double n_nue[num_data];
  double j_nue[num_data];
  double chi_nue[num_data];
  double n_anue[num_data];
  double j_anue[num_data];
  double chi_anue[num_data];
  double n_muon[num_data];
  double j_muon[num_data];
  double chi_muon[num_data];
  double n_amuon[num_data];
  double j_amuon[num_data];
  double chi_amuon[num_data];

  // read in the data file
  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[0] == '#' && i == 0) {
      continue;
    }

    sscanf(line,
           "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &T[i], &mu_e[i], &mu_p[i], &mu_n[i],
           &n_nue[i], &j_nue[i], &chi_nue[i],
           &n_anue[i], &j_anue[i], &chi_anue[i],
           &n_muon[i], &j_muon[i], &chi_muon[i],
           &n_amuon[i], &j_amuon[i], &chi_amuon[i]);
    printf("%lf %lf %lf %lf\n", T[i], mu_p[i], mu_e[i], mu_n[i]);
    i++;
  }

  GreyOpacityParams my_grey_opacity_params;
  my_grey_opacity_params.opacity_flags = opacity_flags_default_none;

  char fileline[1000];
  fptr = fopen(data_filepath, "w");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }
  fputs(
      "#w_t_nue temp_t_nue eta_t_nue w_f_nue temp_f_nue c_f_nue w_t_anue temp_t_anue eta_t_anue w_f_anue temp_f_anue c_f_anue w_t_nux temp_t_nux eta_t_nux w_f_nux temp_f_nux c_f_nux\n",
      fptr);

  for (i = 0; i < num_data; i++) {
    my_grey_opacity_params.eos_pars.mu_e = mu_e[i];
    my_grey_opacity_params.eos_pars.mu_p = mu_p[i];
    my_grey_opacity_params.eos_pars.mu_n = mu_n[i];
    my_grey_opacity_params.eos_pars.temp = T[i];
    my_grey_opacity_params.m1_pars.n[id_nue] = n_nue[i];
    my_grey_opacity_params.m1_pars.J[id_nue] = j_nue[i];
    my_grey_opacity_params.m1_pars.chi[id_nue] = chi_nue[i];
    my_grey_opacity_params.m1_pars.n[id_anue] = n_anue[i];
    my_grey_opacity_params.m1_pars.J[id_anue] = j_anue[i];
    my_grey_opacity_params.m1_pars.chi[id_anue] = chi_anue[i];
    my_grey_opacity_params.m1_pars.n[id_nux] = n_muon[i];
    my_grey_opacity_params.m1_pars.J[id_nux] = j_muon[i];
    my_grey_opacity_params.m1_pars.chi[id_nux] = chi_muon[i];
    my_grey_opacity_params.m1_pars.n[id_anux] = n_muon[i];
    my_grey_opacity_params.m1_pars.J[id_anux] = j_muon[i];
    my_grey_opacity_params.m1_pars.chi[id_anux] = chi_muon[i];

    NuDistributionParams my_nudistributionparams = CalculateDistrParamsFromM1(&my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars, true);

    sprintf(fileline, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            my_nudistributionparams.w_t[id_nue], my_nudistributionparams.temp_t[id_nue], my_nudistributionparams.eta_t[id_nue], my_nudistributionparams.w_f[id_nue],
            my_nudistributionparams.temp_f[id_nue], my_nudistributionparams.c_f[id_nue],
            my_nudistributionparams.w_t[id_anue], my_nudistributionparams.temp_t[id_anue], my_nudistributionparams.eta_t[id_anue], my_nudistributionparams.w_f[id_anue],
            my_nudistributionparams.temp_f[id_anue], my_nudistributionparams.c_f[id_anue],
            my_nudistributionparams.w_t[id_nux], my_nudistributionparams.temp_t[id_nux], my_nudistributionparams.eta_t[id_nux], my_nudistributionparams.w_f[id_nux],
            my_nudistributionparams.temp_f[id_nux], my_nudistributionparams.c_f[id_nux]);
    fputs(fileline, fptr);

  }

  fclose(fptr);
}