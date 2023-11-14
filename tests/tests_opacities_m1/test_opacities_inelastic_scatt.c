//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_inelastic_scattering.c
//  \brief Generate a table for inelastic scattering on electrons/positrons

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "../../src/opacities/opacities.h"
#include "../../include/integration.h"
#include "../../include/distribution.h"
#include "../../include/constants.h"
#include "../../include/functions.h"

int main() {

  printf("============================================== \n");
  printf("Testing opacities for inelastic scattering ... \n");
  printf("============================================== \n");

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_opacities_m1/nurates_CCSN/nurates_1.008E+01.txt";

  char data_fileout[200] = "/tests/tests_opacities_m1/output/m1_opacities_inelastic_scatt.txt";
  char data_filepath[300] = {'\0'};

  strcat(filepath, filedir);
  strcat(filepath, outname);
  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);

  printf("Data_directory: %s\n", filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
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

  // read in the data file
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

  // print input data
  if (false) {
    printf("\n");
    printf("Printing out input table ... %s\n", outname);
    printf("\n");

    printf("Input data:\n");
    printf("E_nu: %lf\n", e_nu);
    for (int i = 0; i < num_data; i++) {
      printf("%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             zone[i], r[i], rho[i], T[i], Ye[i], mu_e[i], mu_hat[i], Yh[i], Ya[i], Yp[i], Yn[i], em_nue[i], l_nue_inv[i], em_anue[i], l_anue_inv[i]);
    }

    printf("\n");
    printf("End of input table.\n");
    printf("\n");
  }

  printf("Test for distribution function implementation:\n");

  printf("Generating quadratures ...\n");
  MyQuadrature my_quadrature_1d = {.nx = 60, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&my_quadrature_1d);
  MyQuadrature my_quadrature_2d = {.nx = 60, .ny = 60, .dim = 2, .type = kGauleg, .x1 = 0., .x2 = 1., .y1 = 0., .y2 = 1.};
  GaussLegendreMultiD(&my_quadrature_2d);
  printf("Quadratures generated.\n");

  // activate only inelastic scatt 
  GreyOpacityParams my_grey_opacity_params;
  my_grey_opacity_params.opacity_flags = opacity_flags_default_none;
  my_grey_opacity_params.opacity_flags.use_inelastic_scatt = 1;

  FILE *file;
  file = fopen(data_filepath, "w+");

  printf("\n");

  printf("Generated tables:\n");
  printf("r diff_distr j0-nue j0-anue j0-nux j-nue j-anue j-nux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa-a-nue kappa-a-anue kappa-a-nux kappa-s-nue kappa-s-anue kappa-s-nux\n");

  fprintf(file, "# r diff_distr j0-nue j0-anue j0-nux j-nue j-anue j-nux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa-a-nue kappa-a-anue kappa-a-nux kappa-s-nue kappa-s-anue kappa-s-nux\n");
  // printf("r: [cm], diff_distr [should be zero, checks Fermi Dirac with NuFTotal], j-nue, ");

  for (int i = 0; i < 102; i++) {

    // populate EOS parameters from table
    my_grey_opacity_params.eos_pars.mu_e = mu_e[i];
    my_grey_opacity_params.eos_pars.mu_p = 0.;
    my_grey_opacity_params.eos_pars.mu_n = mu_hat[i] + kQ;
    my_grey_opacity_params.eos_pars.temp = T[i];
    my_grey_opacity_params.eos_pars.yp = Yp[i];
    my_grey_opacity_params.eos_pars.yn = Yn[i];
    my_grey_opacity_params.eos_pars.nb = rho[i] / kMu;

    // Opacity parameters (corrections all switched off)
    my_grey_opacity_params.opacity_pars = opacity_params_default_none;

    // Distribution parameters
    my_grey_opacity_params.distr_pars = NuEquilibriumParams(&my_grey_opacity_params.eos_pars); // consider neutrino distribution function at equilibrium

    // M1 parameters
    MyQuadratureIntegrand Jvals = NuEnergy(&my_grey_opacity_params.distr_pars);
    MyQuadratureIntegrand nvals = NuNumber(&my_grey_opacity_params.distr_pars);
    my_grey_opacity_params.m1_pars.J[0] = Jvals.integrand[0];
    my_grey_opacity_params.m1_pars.J[1] = Jvals.integrand[1];
    my_grey_opacity_params.m1_pars.J[2] = Jvals.integrand[2];
    my_grey_opacity_params.m1_pars.n[0] = nvals.integrand[0];
    my_grey_opacity_params.m1_pars.n[1] = nvals.integrand[1];
    my_grey_opacity_params.m1_pars.n[2] = nvals.integrand[2];

    double distr_fermi =
        FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp, my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p);
    double distr_nuftot = TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);
    double diff_distr = fabs(distr_fermi - distr_nuftot);

    M1Opacities inelastic_scatt_opacities = ComputeM1Opacities(&my_quadrature_1d, &my_quadrature_2d, &my_grey_opacity_params);
    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r[i], diff_distr,
            inelastic_scatt_opacities.eta_0[id_nue], inelastic_scatt_opacities.eta_0[id_anue], inelastic_scatt_opacities.eta_0[id_nux],
            inelastic_scatt_opacities.eta[id_nue], inelastic_scatt_opacities.eta[id_anue], inelastic_scatt_opacities.eta[id_nux],
            inelastic_scatt_opacities.kappa_0_a[id_nue], inelastic_scatt_opacities.kappa_0_a[id_anue], inelastic_scatt_opacities.kappa_0_a[id_nux],
            inelastic_scatt_opacities.kappa_a[id_nue], inelastic_scatt_opacities.kappa_a[id_anue], inelastic_scatt_opacities.kappa_a[id_nux],
            inelastic_scatt_opacities.kappa_s[id_nue], inelastic_scatt_opacities.kappa_s[id_anue], inelastic_scatt_opacities.kappa_s[id_nux]);
    fprintf(file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r[i], diff_distr,
            inelastic_scatt_opacities.eta_0[id_nue], inelastic_scatt_opacities.eta_0[id_anue], inelastic_scatt_opacities.eta_0[id_nux],
            inelastic_scatt_opacities.eta[id_nue], inelastic_scatt_opacities.eta[id_anue], inelastic_scatt_opacities.eta[id_nux],
            inelastic_scatt_opacities.kappa_0_a[id_nue], inelastic_scatt_opacities.kappa_0_a[id_anue], inelastic_scatt_opacities.kappa_0_a[id_nux],
            inelastic_scatt_opacities.kappa_a[id_nue], inelastic_scatt_opacities.kappa_a[id_anue], inelastic_scatt_opacities.kappa_a[id_nux],
            inelastic_scatt_opacities.kappa_s[id_nue], inelastic_scatt_opacities.kappa_s[id_anue], inelastic_scatt_opacities.kappa_s[id_nux]);
    }

  fclose(file);
  return 0;
}

