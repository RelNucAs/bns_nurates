//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_bremsstrahlung.c
//  \brief Generate a table for bremsstrahlung

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../tests.hpp"

int main()
{

    Kokkos::initialize();

    printf("=================================================== \n");
    printf("Testing opacities for bremmstrahlung ... \n");
    printf("=================================================== \n");

    char filename[200] = "m1_opacities_brem.txt";

    // Opacity flags (activate only pair)
    OpacityFlags opacity_flags = opacity_flags_default_none;
    opacity_flags.use_brem     = 1;

    // Opacity parameters (corrections all switched off)
    OpacityParams opacity_pars = opacity_params_default_none;

    TestM1Opacities(filename, &opacity_flags, &opacity_pars);

    Kokkos::finalize();

    return 0;
}
