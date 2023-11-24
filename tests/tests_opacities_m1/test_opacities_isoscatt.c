//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_isoscatt.c
//  \brief Generate a table for iso-energetic neutrino scattering on nucleons

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

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_opacities_m1/nurates_CCSN/nurates_1.008E+01.txt";

  char data_fileout[200] = "/tests/tests_opacities_m1/output/m1_opacities_isoscatt.txt";
  char data_filepath[300] = {'\0'};

  strcat(filepath, filedir);
  strcat(filepath, outname);
  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);

  printf("=================================================== \n");
  printf("Testing opacities for iso-energetic scattering ...  \n");
  printf("=================================================== \n");

  printf("Data load directory: %s\n", filepath);
  printf("Data save directory: %s\n", data_filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  // -------------------------------------------------------------
  // Load data table

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

  // print input data table
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

  // -------------------------------------------------------------

  printf("Test started:\n");
  printf("Generating quadratures ...\n");
  MyQuadrature my_quadrature_1d = {.nx = 60, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&my_quadrature_1d);
  MyQuadrature my_quadrature_2d = {.nx = 60, .ny = 60, .dim = 2, .type = kGauleg, .x1 = 0., .x2 = 1., .y1 = 0., .y2 = 1.};
  GaussLegendreMultiD(&my_quadrature_2d);
  printf("Quadratures generated.\n");

  // activate only iso-energetic scattering
  GreyOpacityParams my_grey_opacity_params;
  my_grey_opacity_params.opacity_flags = opacity_flags_default_none;
  my_grey_opacity_params.opacity_flags.use_iso = 1;

  FILE *file;
  file = fopen(data_filepath, "w+");

  printf("\n");

  printf("Generated tables:\n");
  printf("r diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");

  fprintf(file, "#  diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");
  // printf("r: [cm], diff_distr [should be zero, checks Fermi Dirac with NuFTotal], j-nue, ");

  double diff_distr = -42.;
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
    for (int idx = 0; idx < total_num_species; idx++) {
      my_grey_opacity_params.m1_pars.n[idx] = nvals.integrand[idx];
      my_grey_opacity_params.m1_pars.J[idx] = Jvals.integrand[idx];
    }


    double distr_fermi =
        FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp, my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p);
    double distr_nuftot = TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);


    if (fabs(distr_fermi - distr_nuftot) > diff_distr) {
      diff_distr = fabs(distr_fermi - distr_nuftot);
    }


    M1Opacities coeffs = ComputeM1Opacities(&my_quadrature_1d, &my_quadrature_2d, &my_grey_opacity_params);

    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r[i], diff_distr,
            coeffs.eta_0[id_nue], coeffs.eta_0[id_anue], coeffs.eta_0[id_nux], coeffs.eta_0[id_anux],
            coeffs.eta[id_nue], coeffs.eta[id_anue], coeffs.eta[id_nux], coeffs.eta[id_anux],
            coeffs.kappa_0_a[id_nue], coeffs.kappa_0_a[id_anue], coeffs.kappa_0_a[id_nux], coeffs.kappa_0_a[id_anux],
            coeffs.kappa_a[id_nue], coeffs.kappa_a[id_anue], coeffs.kappa_a[id_nux], coeffs.kappa_a[id_anux],
            coeffs.kappa_s[id_nue], coeffs.kappa_s[id_anue], coeffs.kappa_s[id_nux], coeffs.kappa_s[id_anux]);
    fprintf(file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r[i], diff_distr,
            coeffs.eta_0[id_nue], coeffs.eta_0[id_anue], coeffs.eta_0[id_nux], coeffs.eta_0[id_anux],
            coeffs.eta[id_nue], coeffs.eta[id_anue], coeffs.eta[id_nux], coeffs.eta[id_anux],
            coeffs.kappa_0_a[id_nue], coeffs.kappa_0_a[id_anue], coeffs.kappa_0_a[id_nux], coeffs.kappa_0_a[id_anux],
            coeffs.kappa_a[id_nue], coeffs.kappa_a[id_anue], coeffs.kappa_a[id_nux], coeffs.kappa_a[id_anux],
            coeffs.kappa_s[id_nue], coeffs.kappa_s[id_anue], coeffs.kappa_s[id_nux], coeffs.kappa_s[id_anux]);
  }

  fclose(file);


  if (diff_distr > 1e-12) {
    printf("\n");
    printf("Error in distribution function: %lf", diff_distr);
    return 1;
  }

  return 0;
}
