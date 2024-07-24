//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_pair.c
//  \brief Generate a table for pair opacities

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "../tests.hpp"

int main() {

  printf("=================================================== \n");
  printf("Testing opacities for pair ... \n");
  printf("=================================================== \n");
  
  char filename[200] = "m1_opacities_pair.txt";

  // Opacity flags (activate only pair)
  OpacityFlags opacity_flags = opacity_flags_default_none;
  opacity_flags.use_pair = 1;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = opacity_params_default_none;

  TestM1Opacities(filename, &opacity_flags, &opacity_pars);

  return 0;
}
