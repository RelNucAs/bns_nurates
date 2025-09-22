//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_opacities_benchmarks.cpp
//  \brief Test of bns_nurates performance

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <random>
#include <chrono>

#include <Kokkos_Core.hpp>

#include "m1_opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

#define num_data 102 // number of data points in input profile

using DevExeSpace   = Kokkos::DefaultExecutionSpace;
using DevMemSpace   = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace  = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;

void TestM1OpacitiesBenchmarks(int nx, int mb_nx, OpacityFlags* opacity_flags,
                               OpacityParams* opacity_pars)
{
    constexpr BS_REAL amu        = kBS_Mu;
    constexpr BS_REAL dm         = kBS_Q;
    constexpr BS_REAL mev_to_erg = kBS_MeV;
    constexpr BS_REAL nm3_to_cm3 = 1e-21;

    // Load data from file
    char filepath[300] = {'\0'};
    char filedir[300]  = SOURCE_DIR;
    char outname[200]  = "/inputs/CCSN/ccsn_1d.txt";

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

    // Store columns here
    int zone[num_data];
    BS_REAL r[num_data];
    BS_REAL rho[num_data];
    BS_REAL temp[num_data];
    BS_REAL Ye[num_data];
    BS_REAL mu_e[num_data];
    BS_REAL mu_hat[num_data];
    BS_REAL Yh[num_data];
    BS_REAL Ya[num_data];
    BS_REAL Yp[num_data];
    BS_REAL Yn[num_data];

    // Read in the data file
    int i = 0;
    char line[1000];
    while (fgets(line, sizeof(line), fptr) != NULL)
    {
        if (line[0] == '#')
        {
            continue;
        }

        if constexpr (std::is_same_v<BS_REAL, float>)
        {
            sscanf(line, "%d %f %f %f %f %f %f %f %f %f %f\n", &zone[i], &r[i],
                   &rho[i], &temp[i], &Ye[i], &mu_e[i], &mu_hat[i], &Yh[i],
                   &Ya[i], &Yp[i], &Yn[i]);
            ,
        }
        else
        {
            sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                   &zone[i], &r[i], &rho[i], &temp[i], &Ye[i], &mu_e[i],
                   &mu_hat[i], &Yh[i], &Ya[i], &Yp[i], &Yn[i]);
            ,
        }
        i++;
    }

    fclose(fptr);

    // Construct arrays based on meshblock dimensions
    const int npts = mb_nx * mb_nx * mb_nx;
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

    std::random_device random_device;
    std::mt19937 generate(random_device());
    std::uniform_int_distribution<> uniform_int_distribution(0, num_data - 1);
    for (int i = 0; i < npts; i++)
    {
        int index   = uniform_int_distribution(generate);
        mb_h_idx(i) = index;

        mb_h_rho(i)    = rho[index];
        mb_h_T(i)      = temp[index];
        mb_h_Ye(i)     = Ye[index];
        mb_h_mu_e(i)   = mu_e[index];
        mb_h_mu_hat(i) = mu_hat[index];
        mb_h_Yp(i)     = Yp[index];
        mb_h_Yn(i)     = Yn[index];
    }

    Kokkos::deep_copy(mb_d_rho, mb_h_rho);
    Kokkos::deep_copy(mb_d_T, mb_h_T);
    Kokkos::deep_copy(mb_d_Ye, mb_h_Ye);
    Kokkos::deep_copy(mb_d_mu_e, mb_h_mu_e);
    Kokkos::deep_copy(mb_d_mu_hat, mb_h_mu_hat);
    Kokkos::deep_copy(mb_d_Yp, mb_h_Yp);
    Kokkos::deep_copy(mb_d_Yn, mb_h_Yn);

    // Kokkos views for OpacityFlags and OpacityParams values
    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_opacity_flags(
        "opacity_flags", 5);
    Kokkos::View<bool*, LayoutWrapper, HostMemSpace> h_opacity_pars(
        "opacity_pars", 8);

    h_opacity_flags(0) = opacity_flags->use_abs_em;
    h_opacity_flags(1) = opacity_flags->use_pair;
    h_opacity_flags(2) = opacity_flags->use_brem;
    h_opacity_flags(3) = opacity_flags->use_inelastic_scatt;
    h_opacity_flags(4) = opacity_flags->use_iso;

    h_opacity_pars(0) = opacity_pars->use_dU;
    h_opacity_pars(1) = opacity_pars->use_dm_eff;
    h_opacity_pars(2) = opacity_pars->use_WM_ab;
    h_opacity_pars(3) = opacity_pars->use_WM_sc;
    h_opacity_pars(4) = opacity_pars->use_decay;
    h_opacity_pars(5) = opacity_pars->use_BRT_brem;
    h_opacity_pars(6) = opacity_pars->use_NN_medium_corr;
    h_opacity_pars(7) = opacity_pars->neglect_blocking;

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_opacity_flags(
        "opacity_flags", 5);
    Kokkos::View<bool*, LayoutWrapper, DevMemSpace> d_opacity_pars(
        "opacity_pars", 8);

    Kokkos::deep_copy(d_opacity_flags, h_opacity_flags);
    Kokkos::deep_copy(d_opacity_pars, h_opacity_pars);

    // Generate quadrature data and store them in Kokkos views
    BS_ASSERT(2 * nx <= n_max, "2*nx is larger than n_max");

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

    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_weights("h_weights",
                                                                  my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_points("h_points",
                                                                 my_quad.nx);

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_weights("d_weights",
                                                                 my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_points("d_points",
                                                                my_quad.nx);

    h_n_quad(0) = my_quad.nx;

    for (int i = 0; i < my_quad.nx; i++)
    {
        h_weights(i) = my_quad.w[i];
        h_points(i)  = my_quad.points[i];
    }

    Kokkos::deep_copy(d_n_quad, h_n_quad);
    Kokkos::deep_copy(d_weights, h_weights);
    Kokkos::deep_copy(d_points, h_points);

    // Kokkos views to store results
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_eta_0(
        "coeffs_eta_0", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_eta(
        "coeffs_eta", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_0_a(
        "coeffs_kappa_0_a", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_a(
        "coeffs_kappa_a", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_s(
        "coeffs_kappa_s", npts, 4);


    auto start = std::chrono::high_resolution_clock::now();
    // Kokkos parallel_for loop
    Kokkos::parallel_for(
        "test_opacities_benchmarks",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, npts),
        KOKKOS_LAMBDA(const int& j) {
            // Populate OpacityFlags data
            GreyOpacityParams my_grey_opacity_params;
            my_grey_opacity_params.opacity_flags = {
                .use_abs_em          = d_opacity_flags(0),
                .use_pair            = d_opacity_flags(1),
                .use_brem            = d_opacity_flags(2),
                .use_inelastic_scatt = d_opacity_flags(3),
                .use_iso             = d_opacity_flags(4)};

            // Populate OpacityParams data
            my_grey_opacity_params.opacity_pars = {
                .use_dU             = d_opacity_pars(0),
                .use_dm_eff         = d_opacity_pars(1),
                .use_WM_ab          = d_opacity_pars(2),
                .use_WM_sc          = d_opacity_pars(3),
                .use_decay          = d_opacity_pars(4),
                .use_BRT_brem       = d_opacity_pars(5),
                .use_NN_medium_corr = d_opacity_pars(6),
                .neglect_blocking   = d_opacity_pars(7)};

            // Populate EOS parameters
            my_grey_opacity_params.eos_pars.mu_e = mb_d_mu_e(j);
            my_grey_opacity_params.eos_pars.mu_p = 0.;
            my_grey_opacity_params.eos_pars.mu_n = mb_d_mu_hat(j) + dm;
            my_grey_opacity_params.eos_pars.temp = mb_d_T(j);
            my_grey_opacity_params.eos_pars.yp   = mb_d_Yp(j);
            my_grey_opacity_params.eos_pars.yn   = mb_d_Yn(j);
            my_grey_opacity_params.eos_pars.nb = mb_d_rho(j) / amu * nm3_to_cm3;

            // Populate distribution parameters
            my_grey_opacity_params.distr_pars = NuEquilibriumParams(
                &my_grey_opacity_params
                     .eos_pars); // consider neutrino distribution function at
                                 // equilibrium

            // Populate M1 parameters
            ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                                 &my_grey_opacity_params.distr_pars,
                                 &my_grey_opacity_params.m1_pars);

            for (int idx = 0; idx < total_num_species; idx++)
            {
                my_grey_opacity_params.m1_pars.chi[idx] = 1. / 3.;
                my_grey_opacity_params.m1_pars.J[idx] =
                    my_grey_opacity_params.m1_pars.J[idx] * mev_to_erg;
            }

            // Populate quadrature structure
            MyQuadrature gpu_quad;
            gpu_quad.nx = d_n_quad(0);
            for (int idx = 0; idx < gpu_quad.nx; idx++)
            {
                gpu_quad.w[idx]      = d_weights(idx);
                gpu_quad.points[idx] = d_points(idx);
            }

            // Compute gray emissivities and opacities
            M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad,
                                                    &my_grey_opacity_params);

            // // Save results to Kokkos views
            // d_coeffs_eta_0(j, id_nue)  = coeffs.eta_0[id_nue];
            // d_coeffs_eta_0(j, id_anue) = coeffs.eta_0[id_anue];
            // d_coeffs_eta_0(j, id_nux)  = coeffs.eta_0[id_nux];
            // d_coeffs_eta_0(j, id_anux) = coeffs.eta_0[id_anux];

            // d_coeffs_eta(j, id_nue)  = coeffs.eta[id_nue];
            // d_coeffs_eta(j, id_anue) = coeffs.eta[id_anue];
            // d_coeffs_eta(j, id_nux)  = coeffs.eta[id_nux];
            // d_coeffs_eta(j, id_anux) = coeffs.eta[id_anux];

            // d_coeffs_kappa_0_a(j, id_nue)  = coeffs.kappa_0_a[id_nue];
            // d_coeffs_kappa_0_a(j, id_anue) = coeffs.kappa_0_a[id_anue];
            // d_coeffs_kappa_0_a(j, id_nux)  = coeffs.kappa_0_a[id_nux];
            // d_coeffs_kappa_0_a(j, id_anux) = coeffs.kappa_0_a[id_anux];

            // d_coeffs_kappa_a(j, id_nue)  = coeffs.kappa_a[id_nue];
            // d_coeffs_kappa_a(j, id_anue) = coeffs.kappa_a[id_anue];
            // d_coeffs_kappa_a(j, id_nux)  = coeffs.kappa_a[id_nux];
            // d_coeffs_kappa_a(j, id_anux) = coeffs.kappa_a[id_anux];

            // d_coeffs_kappa_s(j, id_nue)  = coeffs.kappa_s[id_nue];
            // d_coeffs_kappa_s(j, id_anue) = coeffs.kappa_s[id_anue];
            // d_coeffs_kappa_s(j, id_nux)  = coeffs.kappa_s[id_nux];
            // d_coeffs_kappa_s(j, id_anux) = coeffs.kappa_s[id_anux];
        });

    Kokkos::fence();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double time_taken_seconds              = duration.count();

    // Copy results on host
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_eta_0(
        "h_coeffs_eta_0", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_eta(
        "h_coeffs_eta", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_0_a(
        "h_coeffs_kappa_0_a", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_a(
        "h_coeffs_kappa_a", npts, 4);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_s(
        "h_coeffs_kappa_s", npts, 4);

    Kokkos::deep_copy(h_coeffs_eta_0, d_coeffs_eta_0);
    Kokkos::deep_copy(h_coeffs_eta, d_coeffs_eta);
    Kokkos::deep_copy(h_coeffs_kappa_0_a, d_coeffs_kappa_0_a);
    Kokkos::deep_copy(h_coeffs_kappa_a, d_coeffs_kappa_a);
    Kokkos::deep_copy(h_coeffs_kappa_s, d_coeffs_kappa_s);

    // // Print results
    // for (int i = 0; i < npts; i++)
    // {
    //     printf("%d %d %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le
    //     "
    //            "%.15le %.15le %.15le %.15le %.15le %.15le %.15le "
    //            "%.15le %.15le %.15le %.15le %.15le\n",
    //            i, mb_h_idx(i), h_coeffs_eta_0(i, id_nue),
    //            h_coeffs_eta_0(i, id_anue), h_coeffs_eta_0(i, id_nux),
    //            h_coeffs_eta_0(i, id_anux), h_coeffs_eta(i, id_nue),
    //            h_coeffs_eta(i, id_anue), h_coeffs_eta(i, id_nux),
    //            h_coeffs_eta(i, id_anux), h_coeffs_kappa_0_a(i, id_nue),
    //            h_coeffs_kappa_0_a(i, id_anue),
    //            h_coeffs_kappa_0_a(i, id_nux),
    //            h_coeffs_kappa_0_a(i, id_anux), h_coeffs_kappa_a(i, id_nue),
    //            h_coeffs_kappa_a(i, id_anue), h_coeffs_kappa_a(i, id_nux),
    //            h_coeffs_kappa_a(i, id_anux), h_coeffs_kappa_s(i, id_nue),
    //            h_coeffs_kappa_s(i, id_anue), h_coeffs_kappa_s(i, id_nux),
    //            h_coeffs_kappa_s(i, id_anux));
    // }

    printf("Opacities computed.\n");
    printf("Number of quadrature points = %d\n", nx);
    printf("Meshblock Nx = %d\n", mb_nx);
    printf("[Nx x Nx x Nx] = %d\n", npts);
    printf("Total time taken = %lf\n", time_taken_seconds);
    printf("zone-cycles/second = %e\n", double(npts) / time_taken_seconds);

    return;
}

int main(int argc, char* argv[])
{

    Kokkos::initialize();

    {
        printf("=================================================== \n");
        printf("Benchmark tests for bns_nurates ... \n");
        printf("=================================================== \n");

        int nx, mb_nx;

        if (argc != 3)
        {
            printf("Usage: %s <n_quad> <mesh_size>\n", argv[0]);
            return 1;
        }

        nx    = atoi(argv[1]);
        mb_nx = atoi(argv[2]);

        if (nx == 0)
        {
            std::cerr << "Invalid input!" << std::endl;
            return 1;
        }

        if (nx * 2 > n_max)
        {
            std::cerr << "Number of quadrature points exceeds n_max!"
                      << std::endl;
            std::cerr << "2 * nx = " << 2 * nx << std::endl;
            std::cerr << "n_max = " << n_max << std::endl;
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


        OpacityParams opacity_pars = opacity_params_default_none;
        OpacityFlags opacity_flags = opacity_flags_default_none;

        printf("Only beta processes:\n");
        opacity_flags.use_abs_em = 1, opacity_flags.use_pair = 0,
        opacity_flags.use_brem = 0, opacity_flags.use_inelastic_scatt = 0,
        opacity_flags.use_iso = 0;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");


        printf("Only pair annihilations:\n");
        opacity_flags.use_abs_em = 0, opacity_flags.use_pair = 1,
        opacity_flags.use_brem = 0, opacity_flags.use_inelastic_scatt = 0,
        opacity_flags.use_iso = 0;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");

        printf("Only NN bremsstrahlung:\n");
        opacity_flags.use_abs_em = 0, opacity_flags.use_pair = 0,
        opacity_flags.use_brem = 1, opacity_flags.use_inelastic_scatt = 0,
        opacity_flags.use_iso = 0;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");

        printf("Only NEPS:\n");
        opacity_flags.use_abs_em = 0, opacity_flags.use_pair = 0,
        opacity_flags.use_brem = 0, opacity_flags.use_inelastic_scatt = 1,
        opacity_flags.use_iso = 0;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");

        printf("Only isoenergetic scattering:\n");
        opacity_flags.use_abs_em = 0, opacity_flags.use_pair = 0,
        opacity_flags.use_brem = 0, opacity_flags.use_inelastic_scatt = 0,
        opacity_flags.use_iso = 1;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");

        printf("All the reactions:\n");
        opacity_flags.use_abs_em = 1, opacity_flags.use_pair = 1,
        opacity_flags.use_brem = 1, opacity_flags.use_inelastic_scatt = 1,
        opacity_flags.use_iso = 1;
        TestM1OpacitiesBenchmarks(nx, mb_nx, &opacity_flags, &opacity_pars);
        printf("\n");
    }
    Kokkos::finalize();

    return 0;
}
