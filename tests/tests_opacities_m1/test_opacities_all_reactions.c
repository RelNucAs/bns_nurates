//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_all_reactions.c
//  \brief Generate a table for all reactions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include "../tests.h"

int main() {
  clock_t start, end;
  time_t start_t, end_t;
  double cpu_time_used, cpu_time_used_t;


  printf("=================================================== \n");
  printf("Testing opacities for all reactions ... \n");
  printf("=================================================== \n");
  
  char filename[200] = "m1_opacities_all_reactions.txt";

  // Opacity flags (activate all reactions)
  OpacityFlags opacity_flags = opacity_flags_default_all;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = opacity_params_default_none;

  start = clock();
  start_t = time(NULL);

  TestM1Opacities(filename, &opacity_flags, &opacity_pars);
  
  end_t = time(NULL);
  end = clock();

  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cpu_time_used_t = ((double) (end_t - start_t));

  printf("Elapsed time: %.3lf sec\n", cpu_time_used);
  printf("Elapsed time: %.3lf sec\n", cpu_time_used_t);

  return 0;
}