//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_beta.c
//  \brief test routines for beta processes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../src/constants.h"
#include "../../src/opacities/opacities.h"

void generate_comparison_data(const bool use_dU) {

  // Define path of input file
  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = {'\0'};
  
  if (use_dU == false) {
    strcpy(outname, "/inputs/nurates_CCSN/nurates_1.008E+01.txt");
  } else {
    strcpy(outname, "/inputs/nurates_CCSN/nurates_1.008E+01_dU.txt");
  }

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Input data file: %s\n", filepath);
  printf("\n");

  // Number of entries (must be compatible with # of rows in the file)
  const int n_data = 102;

  // Neutrino energy (the comparison is done at fixed energy)
  double omega; // [MeV]

  // Data arrays (rates are not corrected for the weak magnetism, WM corrections are provided separately)
  int zone[n_data];        // zone (counter)
  double rad[n_data];        // radius [cm]
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
  double r[n_data], r_bar[n_data];    // weak magnetism correction to beta rates (electron-type neutrino and antineutrino, respectvely)
  double b_is[n_data];     // source term for quasi-elastic scattering (Eq.(A41) in Bruenn1985?)
  double r0_n[n_data], r0_p[n_data];  // weak magnetism correction to zeroth Legendre term of scattering kernel (on neutrons and protons, respectively)
  double r1_n[n_data], r1_p[n_data];  // weak magnetism correction to first Legendre term of scattering kernel (on neutrons and protons, respectively)

  // Opacity parameters (correction to rates)
  OpacityParams opacity_pars = {.use_dU = false}; // dU correction turned off
  
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
  fgets(line, sizeof(line), fptr_in);
  sscanf(line + 14, "%lf\n", &omega);
  
  printf("Neutrino energy: %.5e MeV\n", omega);
  printf("\n");

  // Define path of output file
  filepath[0] = '\0';
  if (use_dU == false) {
    sprintf(outname, "/tests/tests_beta/test_beta_%.3e.txt", omega);
  } else {
    sprintf(outname, "/tests/tests_beta/test_beta_%.3e_dU.txt", omega);
  }

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

  // Print legend as output file header
  fprintf(fptr_out, "# 1: zone\n");
  fprintf(fptr_out, "# 2: r [cm]\n");
  fprintf(fptr_out, "# 3: rho [g/cm^3]\n");
  fprintf(fptr_out, "# 4: T [MeV]\n");
  fprintf(fptr_out, "# 5: Ye\n");
  fprintf(fptr_out, "# 6: mu_e [MeV]\n");
  fprintf(fptr_out, "# 7: mu_hat [MeV]\n");
  fprintf(fptr_out, "# 8: Yh\n");
  fprintf(fptr_out, "# 9: Ya\n");
  fprintf(fptr_out, "# 10: Yp\n");
  fprintf(fptr_out, "# 11: Yn\n");
  fprintf(fptr_out, "# 12: dU [MeV]\n");
  fprintf(fptr_out, "# 13-16: reference rates w/o WM (em_nue [s-1], ab_nue [cm-1], em_anue [s-1], ab_anue [cm-1])\n");
  fprintf(fptr_out, "# 17-18: WM correction to reference rates (R_nue, R_anue)\n");
  fprintf(fptr_out, "# 19-22: reference rates w/o WM (em_nue [s-1], ab_nue [cm-1], em_anue [s-1], ab_anue [cm-1])\n");
  fprintf(fptr_out, "# 23-26: reference rates with WM (em_nue [s-1], ab_nue [cm-1], em_anue [s-1], ab_anue [cm-1])\n");
  
  // Read in data from profile, compute rates and store results in output file
  while (fgets(line, sizeof(line), fptr_in) != NULL) {
    if (line[1] == '#') {
      continue;
    }

    sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &zone[i], &rad[i], &rho[i], &temp[i], &ye[i], &mu_e[i], &mu_hat[i],
           &yh[i], &ya[i], &yp[i], &yn[i], &du[i],
           &em_nue[i], &ab_nue[i], &em_anue[i], &ab_anue[i], &r[i], &r_bar[i],
           &b_is[i], &r0_n[i], &r0_p[i], &r1_n[i], &r1_p[i]
           );

    // Populate EOS params structure with data read from the profile
    eos_pars.nb   = rho[i] / kMu;   // baryon number density [cm-3]
    eos_pars.temp = temp[i];        
    eos_pars.yp   = yp[i];          
    eos_pars.yn   = yn[i];          
    eos_pars.mu_p = 0.;            
    eos_pars.mu_n = mu_hat[i] + kQ; // N.B.: corrected for the nucleon mass difference
    eos_pars.mu_e = mu_e[i];        
    eos_pars.dU   = du[i];

    // Compute emissivity and inverse mean free path for electron-type (anti)neutrinos
    opacity_pars.use_WM_ab = false; // WM correction turned off
    MyOpacity out    = AbsOpacity(omega, &opacity_pars, &eos_pars); 
    
    opacity_pars.use_WM_ab = true;  // WM correction turned on
    MyOpacity out_WM = AbsOpacity(omega, &opacity_pars, &eos_pars);


    fprintf(fptr_out, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
             zone[i], rad[i], rho[i], temp[i], ye[i], mu_e[i], mu_hat[i],
             yh[i], ya[i], yp[i], yn[i], du[i],
             em_nue[i], ab_nue[i], em_anue[i], ab_anue[i], r[i], r_bar[i],
             out.em_nue, out.ab_nue / kClight, out.em_anue, out.ab_anue / kClight,
             out_WM.em_nue, out_WM.ab_nue / kClight, out_WM.em_anue, out_WM.ab_anue / kClight);
   
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

  return;
}





int main() {
      printf("Spectral rates are compared with the results of a tested Fortran code\n");
  printf("\n");

  bool use_dU = false;
  
  printf("============================================================== \n");
  printf("Testing opacities for beta processes without dU correction ... \n");
  printf("============================================================== \n");
  printf("\n");

  generate_comparison_data(use_dU); // dU correction is off

  use_dU = true;
  
  printf("=========================================================== \n");
  printf("Testing opacities for beta processes with dU correction ... \n");
  printf("=========================================================== \n");
  printf("\n");

  generate_comparison_data(use_dU); // dU correction is on

  return 0;
}