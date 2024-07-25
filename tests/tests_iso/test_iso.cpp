//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_iso.c
//  \brief test routines for isoenergetic scattering kernel

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/constants.hpp"
#include "../../include/weak_magnetism.hpp"
#include "opacities.hpp"

int main () {
  printf("Kernel is compared with the result of a tested Fortran code\n");
  printf("\n");
 
  printf("========================================================== \n");
  printf("Testing kernel for isoenergetic scattering on nucleons ... \n");
  printf("========================================================== \n");
  printf("\n");

  // Define path of input file
  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[300];
  strcpy(outname, "/inputs/nurates_CCSN/nurates_1.008E+01.txt");

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Input data file: %s\n", filepath);
  printf("\n");

  // Number of entries (must be compatible with # of rows in the file)
  const int n_data = 102;

  // Neutrino energy (the comparison is done at fixed energy)
  double omega; // [MeV]

  // Weak magnetism correction (code output)
  double w0_n_out, w1_n_out, w0_p_out, w1_p_out; 

  // Data arrays (rates are not corrected for the weak magnetism, WM corrections are provided separately)
  int zone[n_data];        // zone (counter)
  double r[n_data];        // radius [cm]
  double rho[n_data];      // mass density [g/cm3]
  double temp[n_data];     // temperature [MeV]
  double ye[n_data];       // electron fraction [#/baryon]
  double mu_e[n_data];     // electron chemical potential (with rest mass) [MeV]
  double mu_hat[n_data];   // neutron - proton chemical potential (w/o rest mass) [MeV]
  double yh[n_data];       // heavy abundance [#/baryon]
  double ya[n_data];       // alpha abundance [#/baryon]
  double yp[n_data];       // proton abundance [#/baryon]
  double yn[n_data];       // neuton abundance [#/baryon]
  double du[n_data];       // correction to nucleon chemical potential due to interactions [MeV]
  double em_nue[n_data];   // electron neutrino emissivity [1/s]
  double ab_nue[n_data];   // electron neutrino inverse mean free path [1/cm]
  double em_anue[n_data];  // electron antineutrino emissivity [1/s]
  double ab_anue[n_data];  // electron antineutrino inverse mean free path [1/cm]
  double w_nue[n_data], w_anue[n_data];    // weak magnetism correction to beta rates (electron-type neutrino and antineutrino, respectvely)
  double r_is[n_data];     // source term for quasi-elastic scattering (Eq.(A41) in Bruenn1985?)
  double w0_n[n_data], w0_p[n_data];  // weak magnetism correction to zeroth Legendre term of scattering kernel (on neutrons and protons, respectively)
  double w1_n[n_data], w1_p[n_data];  // weak magnetism correction to first Legendre term of scattering kernel (on neutrons and protons, respectively)

  // Opacity parameters (correction to rates)
  OpacityParams opacity_pars = {.use_WM_sc = true}; // WM correction turned on
  
  // EOS parameters
  MyEOSParams eos_pars;

  // Open input file
  FILE *fptr_in;
  fptr_in = fopen(filepath, "r");
  if (fptr_in == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  int i = 0;
  char line[1500];

  // Read in the neutrino energy
  auto val = fgets(line, sizeof(line), fptr_in);
  sscanf(line + 14, "%lf\n", &omega);
  
  printf("Neutrino energy: %.5e MeV\n", omega);
  printf("\n");

  // Define path of output file
  filepath[0] = '\0';
  sprintf(outname, "/tests/tests_iso/test_iso_%.3e.txt", omega);
 
  strcat(filepath, filedir);
  strcat(filepath, outname);

  // Open output file
  FILE *fptr_out;
  fptr_out = fopen(filepath, "w");
  if (fptr_out == NULL) {
    printf("%s: Error in opening in write mode file %s!\n", __FILE__, filepath);
    exit(1);
  }

  printf("Output data file: %s\n", filepath);
  printf("\n");

  // Print neutrino energy in output file
  fprintf(fptr_out, "# E_nu = %.5e MeV\n", omega);
  fprintf(fptr_out, "\n");

  // Print legend as output file header
  fprintf(fptr_out, "#  1: zone\n");
  fprintf(fptr_out, "#  2: r [cm]\n");
  fprintf(fptr_out, "#  3: rho [g/cm^3]\n");
  fprintf(fptr_out, "#  4: T [MeV]\n");
  fprintf(fptr_out, "#  5: Ye\n");
  fprintf(fptr_out, "#  6: mu_e [MeV]\n");
  fprintf(fptr_out, "#  7: mu_hat [MeV]\n");
  fprintf(fptr_out, "#  8: Yh\n");
  fprintf(fptr_out, "#  9: Ya\n");
  fprintf(fptr_out, "# 10: Yp\n");
  fprintf(fptr_out, "# 11: Yn\n");
  fprintf(fptr_out, "# 12: Fortran R_iso with WM [cm-1]\n");
  fprintf(fptr_out, "# 13: C R_iso with WM [cm-1]\n");
  fprintf(fptr_out, "\n");
  
  // Read in data from profile, compute rates and store results in output file
  while (fgets(line, sizeof(line), fptr_in) != NULL) {
    if (line[1] == '#') {
      continue;
    }

    sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &zone[i], &r[i], &rho[i], &temp[i], &ye[i], &mu_e[i], &mu_hat[i],
           &yh[i], &ya[i], &yp[i], &yn[i], &du[i],
           &em_nue[i], &ab_nue[i], &em_anue[i], &ab_anue[i], &w_nue[i], &w_anue[i],
           &r_is[i], &w0_n[i], &w1_n[i], &w0_p[i], &w1_p[i]
           );

    if (i == 0) {
      WMScatt(omega, &w0_p_out, &w1_p_out, 1);
      WMScatt(omega, &w0_n_out, &w1_n_out, 2);

      // Print weak magnetism correction in output file
      fprintf(fptr_out, "# Weak magnetism correction (Fortran):\n");
      fprintf(fptr_out, "# W0_n = %lf\n", w0_n[i]);
      fprintf(fptr_out, "# W1_n = %lf\n", w1_n[i]);
      fprintf(fptr_out, "# W0_p = %lf\n", w0_p[i]);
      fprintf(fptr_out, "# W1_p = %lf\n", w1_p[i]);
      fprintf(fptr_out, "# Weak magnetism correction (C):\n");
      fprintf(fptr_out, "# W0_n = %lf\n", w0_n_out);
      fprintf(fptr_out, "# W1_n = %lf\n", w1_n_out);
      fprintf(fptr_out, "# W0_p = %lf\n", w0_p_out);
      fprintf(fptr_out, "# W1_p = %lf\n", w1_p_out);
      fprintf(fptr_out, "\n");
    }

    // Populate EOS params structure with data read from the profile
    eos_pars.nb   = rho[i] / kMu;   // baryon number density [cm-3]
    eos_pars.temp = temp[i];        
    eos_pars.yp   = yp[i];          
    eos_pars.yn   = yn[i];          

    // Compute isoenergetic scatterig kernel
    double out = IsoScattTotal(omega, &opacity_pars, &eos_pars);

    // Adjust to compare the same quantity
    out = out * omega * omega * 4. * kPi / kClight;

    fprintf(fptr_out, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
             zone[i], r[i], rho[i], temp[i], ye[i], mu_e[i], mu_hat[i],
             yh[i], ya[i], yp[i], yn[i],
             r_is[i], out);
   
    i++;
  }

  // Close files
  fclose(fptr_in);
  fclose(fptr_out);

  // Check if file was entirely read
  if (n_data != zone[i-1]) {
    printf("Warning: n_data (%d) different for number of lines in the input file (%d)",
           n_data, zone[i-1]);
  }

  return 0;
}