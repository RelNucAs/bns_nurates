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
//#include "tests.h"
#include "bns_nurates.h"
#include "kernels.h"
#include "integration.h"
#include "../src/opacities/opacities.h"

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
    printf("%s: The file %s does not exist!\n", filepath);
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelS(x[i], 0., 0.31), BremKernelS(x[i], 0., 1.01), BremKernelS(x[i], 0., 2.42),
           BremKernelS(x[i], 2., 0.31), BremKernelS(x[i], 4., 0.31));
    sprintf(fileline, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelS(x[i], 0., 0.31), BremKernelS(x[i], 0., 1.01), BremKernelS(x[i], 0., 2.42),
            BremKernelS(x[i], 2., 0.31), BremKernelS(x[i], 4., 0.31));
    fprintf(fptr, fileline);
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
    printf("%s: The file %s does not exist!\n", filepath);
    exit(1);
  }

  for (int i = 0; i < n; i++) {
    printf("%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelG(x[i], 0.31), BremKernelG(x[i], 1.01), BremKernelG(x[i], 2.42),
           BremKernelG(x[i], 4.22), BremKernelG(x[i], 6.14), BremKernelG(x[i], 8.11));
    sprintf(fileline, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", x[i], BremKernelG(x[i], 0.31), BremKernelG(x[i], 1.01), BremKernelG(x[i], 2.42),
            BremKernelG(x[i], 4.22), BremKernelG(x[i], 6.14), BremKernelG(x[i], 8.11));
    fprintf(fptr, fileline);
  }

  fclose(fptr);
}