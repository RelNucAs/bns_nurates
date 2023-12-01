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

#include "../tests.h"

int main() {

  printf("=================================================== \n");
  printf("Testing opacities for inelastic scattering ... \n");
  printf("=================================================== \n");
  
  char filename[200] = "m1_opacities_inelastic_scatt.txt";

  // Opacity flags (activate only pair)
  OpacityFlags opacity_flags = opacity_flags_default_none;
  opacity_flags.use_inelastic_scatt = 1;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = opacity_params_default_none;

  TestM1Opacities(filename, &opacity_flags, &opacity_pars);
    
  return 0;
}

