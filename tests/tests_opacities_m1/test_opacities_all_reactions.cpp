//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_all_reactions.c
//  \brief Generate a table for all reactions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <Kokkos_Core.hpp>

#include "../tests.hpp"

int main()
{

    Kokkos::initialize();

    printf("=================================================== \n");
    printf("Testing opacities for all reactions ... \n");
    printf("=================================================== \n");

    char filename[200] = "m1_opacities_all_reactions.txt";

    // Opacity flags (activate all reactions)
    OpacityFlags opacity_flags = opacity_flags_default_all;

    // Opacity parameters (corrections all switched off)
    OpacityParams opacity_pars = opacity_params_default_none;

    TestM1Opacities(filename, &opacity_flags, &opacity_pars);

    Kokkos::finalize();

    return 0;
}
