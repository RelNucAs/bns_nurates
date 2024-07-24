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

void ComputeM1CoeffsGivenQuadrature(const int n_leg, const char *reac_type, OpacityFlags *opacity_flags) {

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

  sprintf(write_file, "/tests/tests_quadrature_convergence/output/integrated/bns_nurates_data_%s_quad_%d_new.txt", reac_type, 2 * n_leg);
  printf("%s\n", write_file);

  strcpy(write_path, SOURCE_DIR); 
  strcat(write_path, write_file);

  FILE *fw = fopen(write_path, "w");
  if (fw == NULL) {
    printf("Error opening the file %s\n", write_file);
    return;
  }

  printf("Computing M1 source coefficients for n_leg = 2 * %d....\n", n_leg);

  MyQuadrature quad_1d = {.nx = n_leg, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&quad_1d);

  fprintf(fw, "# Convergence test of M1 source coefficients with number of quadrature points: n = 2 * %d\n", n_leg);
  fprintf(fw, "#\n");

  // Print opacity flags
  fprintf(fw, "# List of included reactions:\n");
  fprintf(fw, "# Neutrino abs on nucleons - e+/e- capture    : %d\n", opacity_flags->use_abs_em);
  fprintf(fw, "# Elastic scattering on nucleons              : %d\n", opacity_flags->use_iso);
  fprintf(fw, "# Nucleon-nucleon bremsstrahlung (and inverse): %d\n", opacity_flags->use_brem);
  fprintf(fw, "# e+ e- annihilation (and inverse)            : %d\n", opacity_flags->use_pair);
  fprintf(fw, "# Inelastic scattering on e+/e-               : %d\n", opacity_flags->use_inelastic_scatt);
  fprintf(fw, "#\n");

  fprintf(fw, "# radius [cm], density [gcc], eta0_nue/anue/nux/anux [cm-3 s-1], eta_nue/anue/nux/anux [MeV cm-3 s-1], kappa_a0_nue/anue/nux/anux [cm-1], kappa_a_nue/anue/nux/anux [cm-1], kappa_s_nue/anue/nux/anux [cm-1]\n");

  MyEOSParams eos_pars;

  for (int i = 0; i < 102; i++) {

    // populate EOS parameters from table
    eos_pars.mu_e = mu_e[i];
    eos_pars.mu_p = 0.;
    eos_pars.mu_n = mu_hat[i] + kQ;
    eos_pars.temp = T[i];
    eos_pars.yp = Yp[i];
    eos_pars.yn = Yn[i];
    eos_pars.nb = rho[i] / kMu;

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

    M1Opacities m1_coeffs = ComputeM1Opacities(&quad_1d, &quad_1d, &grey_pars);

    fprintf(fw, "%.10e ", r[i]);
    fprintf(fw, "%.10e ", rho[i]);
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fw, "%.10e ", m1_coeffs.eta_0[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fw, "%.10e ", m1_coeffs.eta[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fw, "%.10e ", m1_coeffs.kappa_0_a[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fw, "%.10e ", m1_coeffs.kappa_a[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fw, "%.10e ", m1_coeffs.kappa_s[idx]);
    }
    fprintf(fw, "\n");
  }
  fprintf(fw, "\n");

  fclose(fw);
   
  printf("Done!\n");
  printf("\n");

  return;
}

  

int main() {

  // Number of quadrature points
  int n_leg[7] = {5, 10, 20, 30, 40, 50, 60};

  // Opacity flags (to switch on/off reactions)
  OpacityFlags opacity_flags;
  char *reac_type;

  /* Select the reaction */

  // Beta processes 
  opacity_flags = opacity_flags_default_none;
  opacity_flags.use_abs_em = 1;
  reac_type = "beta";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }
  
  // Pair process 
  opacity_flags = opacity_flags_default_none;
  opacity_flags.use_pair = 1;
  reac_type = "pair";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }
  
  // NN bremsstrahlung
  opacity_flags = opacity_flags_default_none;
  opacity_flags.use_brem = 1;
  reac_type = "brem";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }

  // Inelastic scattering 
  opacity_flags = opacity_flags_default_none;
  opacity_flags.use_inelastic_scatt = 1;
  reac_type = "inel";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }
  
  // Isoenergetic scattering 
  opacity_flags = opacity_flags_default_none;
  opacity_flags.use_iso = 1;
  reac_type = "iso";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }
  
  // All reactions except inelastic scattering 
  opacity_flags = opacity_flags_default_all;
  opacity_flags.use_inelastic_scatt = 0;
  reac_type = "noinel";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }
  
  // All reactions 
  opacity_flags = opacity_flags_default_all;
  reac_type = "all";
 
  for (int i = 0; i < 6; i++) {
    ComputeM1CoeffsGivenQuadrature(n_leg[i], reac_type, &opacity_flags);
  }

  return 0;
}
