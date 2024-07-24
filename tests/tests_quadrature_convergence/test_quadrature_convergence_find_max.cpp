//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_quadrature_convergence_integrated.c
//  \brief Test convergence of integrated rates with number of quadrature points

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "functions.hpp"

void ComputeSpectralEmissivities(const int n_leg, const char *reac_type, OpacityFlags *opacity_flags) {

  double s, eta;

  char read_path[300] = {'\0'};
  char read_file[200] = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

  strcpy(read_path, SOURCE_DIR);
  strcat(read_path, read_file);

  printf("Data_directory: %s\n", read_path);

  FILE *fr;
  fr = fopen(read_path, "r");
  if (fr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, read_path);
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
  while (fgets(line, sizeof(line), fr) != NULL) {
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

  fclose(fr);



  char write_path[300] = {'\0'};
  char write_file[200];

  sprintf(write_file, "/tests/tests_quadrature_convergence/output/find_max/bns_nurates_data_%s_quad_%d_3.txt", reac_type, 2 * n_leg);
  printf("%s\n", write_file);

  strcpy(write_path, SOURCE_DIR); 
  strcat(write_path, write_file);

  FILE *fw = fopen(write_path, "w");
  if (fw == NULL) {
    printf("Error opening the file %s\n", write_file);
    return;
  }

  MyQuadrature quad_1d = {.nx = n_leg, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&quad_1d);

  fprintf(fw, "# Computing spectral rates for reaction %s with number of quadrature points: n = 2 * %d\n", reac_type, n_leg);
  fprintf(fw, "#\n");

  // Print opacity flags
  fprintf(fw, "# List of included reactions:\n");
  fprintf(fw, "# Neutrino abs on nucleons - e+/e- capture    : %d\n", opacity_flags->use_abs_em);
  fprintf(fw, "# Elastic scattering on nucleons              : %d\n", opacity_flags->use_iso);
  fprintf(fw, "# Nucleon-nucleon bremsstrahlung (and inverse): %d\n", opacity_flags->use_brem);
  fprintf(fw, "# e+ e- annihilation (and inverse)            : %d\n", opacity_flags->use_pair);
  fprintf(fw, "# Inelastic scattering on e+/e-               : %d\n", opacity_flags->use_inelastic_scatt);
  fprintf(fw, "#\n");

  fprintf(fw, "# radius [cm], density [gcc], temperature [MeV], mu_hat [MeV], mu_e [MeV], j_nue/anue/nux/anux [s-1] as a function of quadrature energies\n");

  MyEOSParams eos_pars;

  for (int i = 0; i < 102; i++) {

    printf("i = %d\n", i);

    // populate EOS parameters from table
    eos_pars.mu_e = mu_e[i];
    eos_pars.mu_p = 0.;
    eos_pars.mu_n = mu_hat[i] + kQ;
    eos_pars.temp = T[i];
    eos_pars.yp = Yp[i];
    eos_pars.yn = Yn[i];
    eos_pars.nb = rho[i] / kMu;

    eta = eos_pars.mu_e / eos_pars.temp;

    s = 0.5 * eos_pars.temp * (FDI_p4(eta) / FDI_p3(eta) + FDI_p4(-eta) / FDI_p3(-eta));

    // Neutrino distribution function parameters
    NuDistributionParams distr_pars = NuEquilibriumParams(&eos_pars); // Fermi-Dirac equilibrium distribution function

    // M1 quantities
    M1Quantities m1_pars;
    MyQuadratureIntegrand Jvals = NuEnergy(&distr_pars);
    MyQuadratureIntegrand nvals = NuNumber(&distr_pars);
    for (int idx = 0; idx < total_num_species; idx++) {
      m1_pars.n[idx] = nvals.integrand[idx];
      m1_pars.J[idx] = Jvals.integrand[idx];
    }

    // Grey opacity parameters
    GreyOpacityParams grey_pars = {.eos_pars = eos_pars, .opacity_flags = *opacity_flags, .opacity_pars = opacity_params_default_none, .distr_pars = distr_pars, .m1_pars = m1_pars};

    fprintf(fw, "%.10e ", r[i]);
    fprintf(fw, "%.10e ", rho[i]);
    fprintf(fw, "%.10e ", T[i]);
    fprintf(fw, "%.10e ", mu_hat[i]);
    fprintf(fw, "%.10e ", mu_e[i]);
    
    for (int j = 0; j < n_leg; j ++) {
        e_nu = s * quad_1d.points[j];
        SpectralOpacities out = ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, &quad_1d, &grey_pars);
 
        for (int idx = 0; idx < total_num_species; idx++) {
          fprintf(fw, "%.10e ", out.j[idx]);
        }
        
    }

    for (int j = 0; j < n_leg; j ++) {
        e_nu = s / quad_1d.points[n_leg-1-j];
        SpectralOpacities out = ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, &quad_1d, &grey_pars);
 
        for (int idx = 0; idx < total_num_species; idx++) {
          fprintf(fw, "%.10e ", out.j[idx]);
        }
        
    }
  
    fprintf(fw, "\n");

  }

  fclose(fw);
   
  printf("Done!\n");
  printf("\n");

  return;
}


void ComputeDoubleIntegrand(const int n_leg, const char *reac_type, OpacityFlags *opacity_flags) {

  double s, eta;
  M1Matrix mat_n, mat_j;

  InitializeM1Matrix(&mat_n, n_leg);
  InitializeM1Matrix(&mat_j, n_leg);

  char read_path[300] = {'\0'};
  char read_file[200] = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

  strcpy(read_path, SOURCE_DIR);
  strcat(read_path, read_file);

  printf("Data_directory: %s\n", read_path);

  FILE *fr;
  fr = fopen(read_path, "r");
  if (fr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, read_path);
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
  while (fgets(line, sizeof(line), fr) != NULL) {
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

  fclose(fr);



  char write_path[300] = {'\0'};
  char write_file[200];

  sprintf(write_file, "/tests/tests_quadrature_convergence/output/find_max/bns_nurates_data_%s_quad_%d_J_2.txt", reac_type, 2 * n_leg);
  printf("%s\n", write_file);

  strcpy(write_path, SOURCE_DIR); 
  strcat(write_path, write_file);

  FILE *fw = fopen(write_path, "w");
  if (fw == NULL) {
    printf("Error opening the file %s\n", write_file);
    return;
  }

  MyQuadrature quad_1d = {.nx = n_leg, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&quad_1d);

  fprintf(fw, "# Computing 2d integrand for reaction %s with number of quadrature points: n = 2 * %d\n", reac_type, n_leg);
  fprintf(fw, "#\n");

  // Print opacity flags
  fprintf(fw, "# List of included reactions:\n");
  fprintf(fw, "# Neutrino abs on nucleons - e+/e- capture    : %d\n", opacity_flags->use_abs_em);
  fprintf(fw, "# Elastic scattering on nucleons              : %d\n", opacity_flags->use_iso);
  fprintf(fw, "# Nucleon-nucleon bremsstrahlung (and inverse): %d\n", opacity_flags->use_brem);
  fprintf(fw, "# e+ e- annihilation (and inverse)            : %d\n", opacity_flags->use_pair);
  fprintf(fw, "# Inelastic scattering on e+/e-               : %d\n", opacity_flags->use_inelastic_scatt);
  fprintf(fw, "#\n");

  MyEOSParams eos_pars;

  fprintf(fw, "# Emissivity integrand as a function of quadrature energies (rows: neutrino energies, columns: antineutrino energies)\n");

  for (int i = 0; i < 102; i++) {

    printf("i = %d\n", i);
    fprintf(fw, "# zone: %d\n", zone[i]);
  
    // populate EOS parameters from table
    eos_pars.mu_e = mu_e[i];
    eos_pars.mu_p = 0.;
    eos_pars.mu_n = mu_hat[i] + kQ;
    eos_pars.temp = T[i];
    eos_pars.yp = Yp[i];
    eos_pars.yn = Yn[i];
    eos_pars.nb = rho[i] / kMu;

    eta = eos_pars.mu_e / eos_pars.temp;

    s = 0.5 * eos_pars.temp * (FDI_p4(eta) / FDI_p3(eta) + FDI_p4(-eta) / FDI_p3(-eta));

    // Neutrino distribution function parameters
    NuDistributionParams distr_pars = NuEquilibriumParams(&eos_pars); // Fermi-Dirac equilibrium distribution function

    // M1 quantities
    M1Quantities m1_pars;
    MyQuadratureIntegrand Jvals = NuEnergy(&distr_pars);
    MyQuadratureIntegrand nvals = NuNumber(&distr_pars);
    for (int idx = 0; idx < total_num_species; idx++) {
      m1_pars.n[idx] = nvals.integrand[idx];
      m1_pars.J[idx] = Jvals.integrand[idx];
    }

    // Grey opacity parameters
    GreyOpacityParams grey_pars = {.eos_pars = eos_pars, .opacity_flags = *opacity_flags, .opacity_pars = opacity_params_default_none, .distr_pars = distr_pars, .m1_pars = m1_pars};
    
    ComputeM1DoubleIntegrand(&quad_1d, &grey_pars, s, &mat_n, &mat_j);

    for (int j = 0; j < total_num_species; j++) {
      fprintf(fw, "# neutrino species: %d\n", j);

      for (int k = 0; k < n_leg; k++) {
        for (int l = 0; l < n_leg; l++) {
          fprintf(fw, "%.10e ", mat_j.m1_mat_em[j][k][l]);
        }
        for (int l = 0; l < n_leg; l++) {
          fprintf(fw, "%.10e ", mat_j.m1_mat_em[j][k][n_leg-1-l]);
        }
        fprintf(fw, "\n");
      }

      for (int k = 0; k < n_leg; k++) {
        for (int l = 0; l < n_leg; l++) {
          fprintf(fw, "%.10e ", mat_j.m1_mat_em[j][n_leg-1-k][l]);
        }
        for (int l = 0; l < n_leg; l++) {
          fprintf(fw, "%.10e ", mat_j.m1_mat_em[j][n_leg-1-k][n_leg-1-l]);
        }
        fprintf(fw, "\n");
      }
      fprintf(fw, "\n");
    }
 
  }

  fclose(fw);
   
  FreeM1Matrix(&mat_n, n_leg);
  FreeM1Matrix(&mat_j, n_leg);

  printf("Done!\n");
  printf("\n");

  return;
}

int main() {

  // Opacity flags
  OpacityFlags opacity_flags = {.use_abs_em = 0, .use_iso = 0, .use_brem = 0, .use_pair = 1, .use_inelastic_scatt = 0};
  const char* reac_type = "pair";
 
  //ComputeSpectralEmissivities(50, reac_type, &opacity_flags);
  ComputeDoubleIntegrand(50, reac_type, &opacity_flags);

  return 0;
}