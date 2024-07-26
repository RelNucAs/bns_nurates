//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_brem.c
//  \brief test routines for the bremsstrahlung process

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.hpp"
//#include "tests.h"
#include "bns_nurates.hpp"
#include "../../include/kernels.hpp"
#include "integration.hpp"
#include "opacities.hpp"
#include "m1_opacities.hpp"

// Write data for plotting Fig. 3 of Hannestadt and Raffelt
void TestBremKernelS(char *filedir) {

  FILE *fptr;

  char fileline[200];
  char outname[50];
  char filepath[200] = {'\0'};

  sprintf(outname, "/brem_kernel_s.txt");
  strcat(filepath, filedir);
  strcat(filepath, outname);

  int n = 100;
  double x[n];

  x[0] = 0.;
  x[n - 1] = 16.;
  double h = (x[n - 1] - x[0]) / ((double) n - 1);

  for (int i = 1; i < n - 1; i++) {
    x[i] = x[0] + i * h;
  }

  fptr = fopen(filepath, "w");
  if (fptr == NULL) {
    printf("The file %s does not exist!\n", filepath);
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelS(x[i], 0., 0.31), BremKernelS(x[i], 0., 1.01), BremKernelS(x[i], 0., 2.42),
           BremKernelS(x[i], 2., 0.31), BremKernelS(x[i], 4., 0.31));
    sprintf(fileline, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelS(x[i], 0., 0.31), BremKernelS(x[i], 0., 1.01), BremKernelS(x[i], 0., 2.42),
            BremKernelS(x[i], 2., 0.31), BremKernelS(x[i], 4., 0.31));
    fputs(fileline, fptr);
  }

  fclose(fptr);

}

// Write data for plotting Fig. 4 of Hannestadt and Raffelt
void TestBremKernelG(char *filedir) {

  FILE *fptr;

  char fileline[200];
  char outname[50];
  char filepath[200] = {'\0'};

  sprintf(outname, "/brem_kernel_g.txt");
  strcat(filepath, filedir);
  strcat(filepath, outname);

  int n = 100;
  double x[n];

  x[0] = 0.;
  x[n - 1] = 8.;
  double h = (x[n - 1] - x[0]) / ((double) n - 1);

  for (int i = 1; i < n - 1; i++) {
    x[i] = x[0] + i * h;
  }

  fptr = fopen(filepath, "w");
  if (fptr == NULL) {
    printf("The file %s does not exist!\n", filepath);
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelG(x[i], 0.31), BremKernelG(x[i], 1.01), BremKernelG(x[i], 2.42),
           BremKernelG(x[i], 4.22), BremKernelG(x[i], 6.14), BremKernelG(x[i], 8.11));
    sprintf(fileline, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelG(x[i], 0.31), BremKernelG(x[i], 1.01), BremKernelG(x[i], 2.42),
            BremKernelG(x[i], 4.22), BremKernelG(x[i], 6.14), BremKernelG(x[i], 8.11));
    fputs(fileline, fptr);
  }

  fclose(fptr);
}

int main() {
  // @TODO: put here values for comparison
  // Comparison values
  const double ref_em_ker  = 6.29472e-30; // [cm^3 s^-1]
  const double ref_abs_ker = 7.60512e-29; // [cm^3 s^-1]
  
  // Kernel parameters
  BremKernelParams kernel_pars;

  kernel_pars.omega = 10.; // [MeV]
  kernel_pars.omega_prime = 20.; // [MeV]
  kernel_pars.l = 0;
  
  // EOS parameters
  MyEOSParams eos_pars;

  // @TODO: decide a global baryon mass to convert number to mass density and vice versa
  eos_pars.nb   = 3.73e+14 / kMb; // [cm-3]
  eos_pars.temp = 12.04;  // [MeV]
  eos_pars.yn   = 0.6866; 
  eos_pars.yp   = 0.3134;
 
  // Compute kernels
  MyKernelOutput out = BremKernelsLegCoeff(&kernel_pars, &eos_pars);
  
  // Print kernels
  printf("Output emission   kernel = %.5e cm^3/s\n", out.em[id_nue]);
  printf("Output absorption kernel = %.5e cm^3/s\n", out.abs[id_nue]);
  printf("\n");

  // Compare kernels
  const double tol = 1.0E-06; // tolerance on relative difference
  const double em_diff  = fabs(out.em[id_nue]  - ref_em_ker ) / ref_em_ker;  // emission kernel difference
  const double abs_diff = fabs(out.abs[id_nue] - ref_abs_ker) / ref_abs_ker; // absorption kernel difference

  if ( (em_diff > tol) || (abs_diff > tol) ) {
    printf("Test failed!\n");
    printf("Comparison values are:\n");
    printf("Emission   kernel = %.5e cm^3/s\n", ref_em_ker);
    printf("Absorption kernel = %.5e cm^3/s\n", ref_abs_ker);
  } else {
    printf("Test passed!\n");
  } 

  return 0;
}
