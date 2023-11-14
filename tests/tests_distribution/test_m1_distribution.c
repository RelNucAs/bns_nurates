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
#include "../../src/opacities/opacities.h"
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
  char outname[200] = "/tests/tests_distribution/data/m1_reconstruction.txt";

  char data_fileout[200] = "/tests/tests_distribution/data/m1_distribution.txt";
  char data_filepath[300] = {'\0'};

  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data_directory: %s\n", filepath);

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
  double n_nue[num_data];
  double j_nue[num_data];
  double chi_nue[num_data];
  double mu_e[num_data];
  double mu_p[num_data];
  double mu_n[num_data];

  // read in the data file
  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[0] == '#' && i == 0) {
      continue;
    }

    sscanf(line, "%lf %lf %lf %lf %lf %lf %lf\n",
           &n_nue[i], &j_nue[i], &chi_nue[i], &T[i], &mu_e[i], &mu_p[i], &mu_n[i]);
    i++;
  }

  GreyOpacityParams my_grey_opacity_params;
  my_grey_opacity_params.opacity_flags = opacity_flags_default_none;

  for (i = 0; i < num_data; i++) {
    my_grey_opacity_params.eos_pars.mu_e = mu_e[i];
    my_grey_opacity_params.eos_pars.mu_p = mu_p[i];
    my_grey_opacity_params.eos_pars.mu_n = mu_n[i];
    my_grey_opacity_params.eos_pars.temp = T[i];
    my_grey_opacity_params.m1_pars.n[0] = n_nue[i];
    my_grey_opacity_params.m1_pars.J[0] = j_nue[i];
    my_grey_opacity_params.m1_pars.chi[0] = chi_nue[i];

    NuDistributionParams my_nudistributionparams = CalculateDistrParamsFromM1(&my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

    /*
    printf("cf = %lf, wf = %lf, temp_f = %e, eta_t = %lf, w_t = %lf, temp_t = %e\n",
           my_nudistributionparams.c_f[0],
           my_nudistributionparams.w_f[0],
           my_nudistributionparams.temp_f[0],
           my_nudistributionparams.eta_t[0],
           my_nudistributionparams.w_t[0],
           my_nudistributionparams.temp_t[0]);
    */

    my_nudistributionparams.w_t[0] = 0.;
    double f_nue = TotalNuF(energy_vals[16], &my_nudistributionparams, 0); // @TODO: fix
    printf("%e\n", 4. * M_PI * f_nue);

  }

}