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
#include <iostream>
#include <random>
#include <chrono>

#include <Kokkos_Core.hpp>

#include "../tests.hpp"
#include "m1_opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

using DevExeSpace   = Kokkos::DefaultExecutionSpace;
using DevMemSpace   = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace  = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;

void TestM1OpacitiesBenchmarks(int nx, int mb_nx)
{

    // Load data from file
    char filepath[300] = {'\0'};
    char filedir[300]  = SOURCE_DIR;
    char outname[200]  = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

    strcat(filepath, filedir);
    strcat(filepath, outname);

    printf("Data_directory: %s\n", filepath);

    FILE* fptr;
    fptr = fopen(filepath, "r");
    if (fptr == NULL)
    {
        printf("%s: The file %s does not exist!\n", __FILE__, filepath);
        exit(1);
    }

    // store columns here
    int num_data = 102;
    BS_REAL e_nu;
    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_zone("zone", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_r("r", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_rho("rho", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_T("T", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Ye("Ye", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_e("mu_e",
                                                               num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_hat("mu_hat",
                                                                 num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yh("Yh", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Ya("Ya", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yp("Yp", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yn("Yn", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_em_nue("em_nue",
                                                                 num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_l_nue_inv("l_nue_inv",
                                                                    num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_em_anue("em_anue",
                                                                  num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_l_anue_inv(
        "l_anue_inv", num_data);

    // read in the data file
    int i = 0;
    char line[1000];
    while (fgets(line, sizeof(line), fptr) != NULL)
    {
        if (line[1] == '#' && i == 0)
        {
            if constexpr (std::is_same_v<BS_REAL, float>)
            {
	        sscanf(line + 14, "%f\n", &e_nu);
            }
	    else
	    {
                sscanf(line + 14, "%lf\n", &e_nu);
            }
	    continue;
        }
        else if (line[1] == '#' && i != 0)
        {
            continue;
        }

        if constexpr (std::is_same_v<BS_REAL, float>)
        {
            sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                   &h_zone(i), &h_r(i), &h_rho(i), &h_T(i), &h_Ye(i), &h_mu_e(i),
                   &h_mu_hat(i), &h_Yh[i], &h_Ya[i], &h_Yp(i), &h_Yn(i),
                   &h_em_nue(i), &h_l_nue_inv(i), &h_em_anue(i), &h_l_anue_inv(i));
        }
	else
        {
	    sscanf(line,
                   "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                   &h_zone(i), &h_r(i), &h_rho(i), &h_T(i), &h_Ye(i), &h_mu_e(i),
                   &h_mu_hat(i), &h_Yh[i], &h_Ya[i], &h_Yp(i), &h_Yn(i),
                   &h_em_nue(i), &h_l_nue_inv(i), &h_em_anue(i), &h_l_anue_inv(i));
        }
	i++;
    }

    fclose(fptr);

    // construct arrays based on meshblock dimensions
    int npts = mb_nx * mb_nx * mb_nx;
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_rho("mb_h_rho",
                                                                 npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_T("mb_h_T", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_Ye("mb_h_Ye",
                                                                npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_mu_e("mb_h_mu_e",
                                                                  npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_mu_hat(
        "mb_h_mu_hat", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_Yp("mb_h_Yp",
                                                                npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> mb_h_Yn("mb_h_Yn",
                                                                npts);

    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_rho("mb_d_rho",
                                                                npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_T("mb_d_T", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_Ye("mb_d_Ye", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_mu_e("mb_d_mu_e",
                                                                 npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_mu_hat(
        "mb_d_mu_hat", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_Yp("mb_d_Yp", npts);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> mb_d_Yn("mb_d_Yn", npts);

    printf("Building data arrays %d points ...\n", npts);

    std::random_device random_device;
    std::mt19937 generate(random_device());
    std::uniform_int_distribution<> uniform_int_distribution(0, num_data - 1);
    for (int i = 0; i < npts; i++)
    {
        int index = uniform_int_distribution(generate);
        // printf("i: %d, random index: %d\n", i, index);

        mb_h_rho(i)    = h_rho(index);
        mb_h_T(i)      = h_T(index);
        mb_h_Ye(i)     = h_Ye(index);
        mb_h_mu_e(i)   = h_mu_e(index);
        mb_h_mu_hat(i) = h_mu_hat(index);
        mb_h_Yp(i)     = h_Yp(index);
        mb_h_Yn(i)     = h_Yn(index);
    }

    Kokkos::deep_copy(mb_d_rho, mb_h_rho);
    Kokkos::deep_copy(mb_d_T, mb_h_T);
    Kokkos::deep_copy(mb_d_Ye, mb_h_Ye);
    Kokkos::deep_copy(mb_d_mu_e, mb_h_mu_e);
    Kokkos::deep_copy(mb_d_mu_hat, mb_h_mu_hat);
    Kokkos::deep_copy(mb_d_Yp, mb_h_Yp);
    Kokkos::deep_copy(mb_d_Yn, mb_h_Yn);

    printf("Data arrays generated.\n");

    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_diff_distribution(
        "diff_distribution", num_data);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_eta_0(
        "coeffs_eta_0", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_eta(
        "coeffs_eta", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_0_a(
        "coeffs_kappa_0_a", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_a(
        "coeffs_kappa_a", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_s(
        "coeffs_kappa_s", num_data, 4);

    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_diff_distribution(
        "h_diff_distribution", num_data);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_eta_0(
        "h_coeffs_eta_0", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_eta(
        "h_coeffs_eta", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_0_a(
        "h_coeffs_kappa_0_a", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_a(
        "h_coeffs_kappa_a", num_data, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_s(
        "h_coeffs_kappa_s", num_data, 4);

    printf("Generating quadratures with %d points ...\n", nx);
    MyQuadrature my_quad = {.type   = kGauleg,
                            .alpha  = -42.,
                            .dim    = 1,
                            .nx     = nx,
                            .ny     = 1,
                            .nz     = 1,
                            .x1     = 0.,
                            .x2     = 1.,
                            .y1     = -42.,
                            .y2     = -42.,
                            .z1     = -42.,
                            .z2     = -42.,
                            .points = {0},
                            .w      = {0}};
    GaussLegendre(&my_quad);
    printf("Quadratures generated.\n");

    constexpr BS_REAL delta_m    = kBS_Q;
    constexpr BS_REAL mu         = kBS_Mu;
    constexpr BS_REAL nm3_to_cm3 = 1e-21;

    auto start = std::chrono::high_resolution_clock::now();
    Kokkos::parallel_for(
        "test_opacities_benchmarks_compute_opacities_loop",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, npts),
        KOKKOS_LAMBDA(const int& i) {
            // printf("Entered in Kokkos parallel_for loop\n");

            GreyOpacityParams my_grey_opacity_params;
            my_grey_opacity_params.opacity_flags = {.use_abs_em          = 1,
                                                    .use_pair            = 1,
                                                    .use_brem            = 1,
                                                    .use_inelastic_scatt = 1,
                                                    .use_iso             = 1};

            // populate EOS parameters from table
            my_grey_opacity_params.eos_pars.mu_e = mb_d_mu_e(i);
            my_grey_opacity_params.eos_pars.mu_p = 0.;
            my_grey_opacity_params.eos_pars.mu_n = mb_d_mu_hat(i) + delta_m;
            my_grey_opacity_params.eos_pars.temp = mb_d_T(i);
            my_grey_opacity_params.eos_pars.yp   = mb_d_Yp(i);
            my_grey_opacity_params.eos_pars.yn   = mb_d_Yn(i);
            my_grey_opacity_params.eos_pars.nb = mb_d_rho(i) / mu * nm3_to_cm3;


            // Opacity parameters (corrections all switched off)
            my_grey_opacity_params.opacity_pars = {.use_dU             = false,
                                                   .use_dm_eff         = false,
                                                   .use_WM_ab          = false,
                                                   .use_WM_sc          = false,
                                                   .use_decay          = false,
                                                   .use_BRT_brem       = false,
                                                   .use_NN_medium_corr = false,
                                                   .neglect_blocking   = false};


            // Distribution parameters
            my_grey_opacity_params.distr_pars = NuEquilibriumParams(
                &my_grey_opacity_params
                     .eos_pars); // consider neutrino distribution function at
                                 // equilibrium

            // M1 parameters
            ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                                 &my_grey_opacity_params.distr_pars,
                                 &my_grey_opacity_params.m1_pars);

            for (int idx = 0; idx < total_num_species; idx++)
            {
                my_grey_opacity_params.m1_pars.chi[idx] = 1. / 3.;
                my_grey_opacity_params.m1_pars.J[idx] =
                    my_grey_opacity_params.m1_pars.J[idx] * kBS_MeV;
            }

            // printf("Computing M1 coefficients\n");

            M1Opacities coeffs =
                ComputeM1Opacities(&my_quad, &my_quad, &my_grey_opacity_params);
            auto testval = coeffs.eta[id_nue];
	    printf("%.15e\n", testval);
        });

    Kokkos::fence();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double time_taken_seconds              = duration.count();

    printf("Opacities computed.\n");
    printf("Number of quadrature points = %d\n", nx);
    printf("Meshblock Nx = %d\n", mb_nx);
    printf("[Nx x Nx x Nx] = %d\n", npts);
    printf("Total time taken = %lf\n", time_taken_seconds);
    printf("zone-cycles/second = %e\n", double(npts) / time_taken_seconds);
}

int main(int argc, char *argv[])
{

    Kokkos::initialize();

    printf("=================================================== \n");
    printf("Benchmark tests for bns_nurates ... \n");
    printf("=================================================== \n");

    int nx, mb_nx;

    if (argc != 3)
    {
        printf("Usage: %s <n_quad> <mesh_size>\n", argv[0]);
	return 1;
    }

    nx = atoi(argv[1]);
    mb_nx = atoi(argv[2]);

    if (nx == 0)
    {
         std::cerr << "Invalid input!" << std::endl;
         return 1;
    }

    if (nx * 2 > n_max)
    {
         std::cerr << "Number of quadrature points exceeds n_max!" << std::endl;
         std::cerr << "2 * nx = " << 2 * nx << std::endl;
         std::cerr << "n_max = "  <<  n_max << std::endl;
         return 1;
    }

    if (mb_nx == 0)
    {
         std::cerr << "Invalid input!" << std::endl;
         return 1;
    }
    
    std::cout << "Number of quadrature points: ";
    std::cout << nx;
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout << "Nx for meshblock with [Nx x Nx x Nx] points: ";
    std::cout << mb_nx;
    std::cout << std::endl;
    std::cout << std::endl;

    TestM1OpacitiesBenchmarks(nx, mb_nx);

    Kokkos::finalize();

    return 0;
}
