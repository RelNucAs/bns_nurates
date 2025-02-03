//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests.h
//  \brief header file for all test routines

#ifndef BNS_NURATES_TESTS_TESTS_HPP_
#define BNS_NURATES_TESTS_TESTS_HPP_

#include "bns_nurates.hpp"
#include "distribution.hpp"
#include "integration.hpp"
#include "m1_opacities.hpp"

using DevExeSpace  = Kokkos::DefaultExecutionSpace;
using DevMemSpace  = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace = Kokkos::HostSpace;
using LayoutWrapper =
    Kokkos::LayoutRight; // increments last index fastest views defined like

inline void TestM1Opacities(char filename[200], OpacityFlags* opacity_flags,
                            OpacityParams* opacity_pars)
{

    char filepath[300] = {'\0'};
    char filedir[300]  = SOURCE_DIR;
    char outname[200]  = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

    strcat(filepath, filedir);
    strcat(filepath, outname);

    printf("# Data_directory: %s\n", filepath);

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

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_zone("zone", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_r("r", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_rho("rho", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_T("T", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Ye("Ye", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_e("mu_e", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_hat("mu_hat",
                                                                num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yh("Yh", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Ya("Ya", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yp("Yp", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yn("Yn", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_em_nue("em_nue",
                                                                num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_l_nue_inv("l_nue_inv",
                                                                   num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_em_anue("em_anue",
                                                                 num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_l_anue_inv(
        "l_anue_inv", num_data);


    // read in the data file
    int i = 0;
    char line[1000];
    while (fgets(line, sizeof(line), fptr) != NULL)
    {
        if (line[1] == '#' && i == 0)
        {
	    if (std::is_same_v<BS_REAL, float>)
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

	if (std::is_same_v<BS_REAL, float>)
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

    printf("# Test for distribution function implementation:\n");

    Kokkos::deep_copy(d_zone, h_zone);
    Kokkos::deep_copy(d_r, h_r);
    Kokkos::deep_copy(d_rho, h_rho);
    Kokkos::deep_copy(d_T, h_T);
    Kokkos::deep_copy(d_Ye, h_Ye);
    Kokkos::deep_copy(d_mu_e, h_mu_e);
    Kokkos::deep_copy(d_mu_hat, h_mu_hat);
    Kokkos::deep_copy(d_Yh, h_Yh);
    Kokkos::deep_copy(d_Ya, h_Ya);
    Kokkos::deep_copy(d_Yp, h_Yp);
    Kokkos::deep_copy(d_Yn, h_Yn);
    Kokkos::deep_copy(d_em_nue, h_em_nue);
    Kokkos::deep_copy(d_em_anue, h_em_anue);
    Kokkos::deep_copy(d_l_nue_inv, h_l_nue_inv);
    Kokkos::deep_copy(d_l_anue_inv, h_l_anue_inv);

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

    printf("# Generating quadratures ...\n");
    MyQuadrature my_quad = {.type   = kGauleg,
                            .alpha  = -42.,
                            .dim    = 1,
                            .nx     = 10,
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
    printf("# Quadratures generated.\n");

    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_weights("h_weights",
                                                                  my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_points("h_points",
                                                                 my_quad.nx);

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_weights("weights",
                                                                 my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_points("points",
                                                                my_quad.nx);

    printf("# Copying quadrature data to Kokkos view on host\n");

    h_n_quad(0) = my_quad.nx;

    for (int i = 0; i < my_quad.nx; i++)
    {
        h_weights(i) = my_quad.w[i];
        h_points(i)  = my_quad.points[i];
    }

    printf("# Deep copying quadrature data from host to device\n");

    Kokkos::deep_copy(d_n_quad, h_n_quad);
    Kokkos::deep_copy(d_weights, h_weights);
    Kokkos::deep_copy(d_points, h_points);

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

    printf("# Generated tables:\n");
    printf("# r diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux "
           "j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux "
           "kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue "
           "kappa-s-anue kappa-s-nux kappa-s-anux\n");

    printf("# Entering Kokkos parallel_for loop\n");

    Kokkos::parallel_for(
        "loop_over_ccsn", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data),
        KOKKOS_LAMBDA(const int& i) {
            printf("# Entered in Kokkos parallel_for loop\n");

            GreyOpacityParams my_grey_opacity_params;

            my_grey_opacity_params.opacity_flags = {
                .use_abs_em          = d_opacity_flags(0),
                .use_pair            = d_opacity_flags(1),
                .use_brem            = d_opacity_flags(2),
                .use_inelastic_scatt = d_opacity_flags(3),
                .use_iso             = d_opacity_flags(4)};

            // Opacity parameters (corrections all switched off)
            my_grey_opacity_params.opacity_pars = {
                .use_dU             = d_opacity_pars(0),
                .use_dm_eff         = d_opacity_pars(1),
                .use_WM_ab          = d_opacity_pars(2),
                .use_WM_sc          = d_opacity_pars(3),
                .use_decay          = d_opacity_pars(4),
                .use_BRT_brem       = d_opacity_pars(5),
                .use_NN_medium_corr = d_opacity_pars(6),
                .neglect_blocking   = d_opacity_pars(7)};


            // populate EOS parameters from table
            my_grey_opacity_params.eos_pars.mu_e = d_mu_e(i);
            my_grey_opacity_params.eos_pars.mu_p = 0.;
            my_grey_opacity_params.eos_pars.mu_n = d_mu_hat(i) + kBS_Q;
            my_grey_opacity_params.eos_pars.temp = d_T(i);
            my_grey_opacity_params.eos_pars.yp   = d_Yp(i);
            my_grey_opacity_params.eos_pars.yn   = d_Yn(i);
            my_grey_opacity_params.eos_pars.nb   = d_rho(i) / kBS_Mu * 1e-21;

            // Distribution parameters
            my_grey_opacity_params.distr_pars = NuEquilibriumParams(
                &my_grey_opacity_params
                     .eos_pars); // consider neutrino distribution function at
                                 // equilibrium

            // M1 parameters
            ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                                 &my_grey_opacity_params.distr_pars,
                                 &my_grey_opacity_params.m1_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                my_grey_opacity_params.m1_pars.chi[idx] = 1. / 3.;
                my_grey_opacity_params.m1_pars.J[idx] =
                    my_grey_opacity_params.m1_pars.J[idx] * kBS_MeV;
            }

            printf("# Generating and populating quadrature on GPU\n");

            MyQuadrature gpu_quad;
            gpu_quad.nx = d_n_quad(0);
            for (int idx = 0; idx < gpu_quad.nx; idx++)
            {
                gpu_quad.w[idx]      = d_weights(idx);
                gpu_quad.points[idx] = d_points(idx);
            }

            BS_REAL distr_fermi =
                FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp,
                           my_grey_opacity_params.eos_pars.mu_e -
                               my_grey_opacity_params.eos_pars.mu_n +
                               my_grey_opacity_params.eos_pars.mu_p);
            BS_REAL distr_nuftot =
                TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);
            BS_REAL diff_distr = Kokkos::fabs(distr_fermi - distr_nuftot);

            printf("# Computing M1 coefficients\n");

            M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad,
                                                    &my_grey_opacity_params);

            printf("# Copying M1 coefficients to Kokkos view on device\n");
            d_diff_distribution(i) = diff_distr;

            d_coeffs_eta_0(i, id_nue)  = coeffs.eta_0[id_nue];
            d_coeffs_eta_0(i, id_anue) = coeffs.eta_0[id_anue];
            d_coeffs_eta_0(i, id_nux)  = coeffs.eta_0[id_nux];
            d_coeffs_eta_0(i, id_anux) = coeffs.eta_0[id_anux];

            d_coeffs_eta(i, id_nue)  = coeffs.eta[id_nue];
            d_coeffs_eta(i, id_anue) = coeffs.eta[id_anue];
            d_coeffs_eta(i, id_nux)  = coeffs.eta[id_nux];
            d_coeffs_eta(i, id_anux) = coeffs.eta[id_anux];

            d_coeffs_kappa_0_a(i, id_nue)  = coeffs.kappa_0_a[id_nue];
            d_coeffs_kappa_0_a(i, id_anue) = coeffs.kappa_0_a[id_anue];
            d_coeffs_kappa_0_a(i, id_nux)  = coeffs.kappa_0_a[id_nux];
            d_coeffs_kappa_0_a(i, id_anux) = coeffs.kappa_0_a[id_anux];

            d_coeffs_kappa_a(i, id_nue)  = coeffs.kappa_a[id_nue];
            d_coeffs_kappa_a(i, id_anue) = coeffs.kappa_a[id_anue];
            d_coeffs_kappa_a(i, id_nux)  = coeffs.kappa_a[id_nux];
            d_coeffs_kappa_a(i, id_anux) = coeffs.kappa_a[id_anux];

            d_coeffs_kappa_s(i, id_nue)  = coeffs.kappa_s[id_nue];
            d_coeffs_kappa_s(i, id_anue) = coeffs.kappa_s[id_anue];
            d_coeffs_kappa_s(i, id_nux)  = coeffs.kappa_s[id_nux];
            d_coeffs_kappa_s(i, id_anux) = coeffs.kappa_s[id_anux];
        });

    printf("# Back on CPU, waiting for Kokkos fence\n");

    Kokkos::fence();

    printf("# Deep copying output from device to host\n");

    Kokkos::deep_copy(h_diff_distribution, d_diff_distribution);
    Kokkos::deep_copy(h_coeffs_eta_0, d_coeffs_eta_0);
    Kokkos::deep_copy(h_coeffs_eta, d_coeffs_eta);
    Kokkos::deep_copy(h_coeffs_kappa_0_a, d_coeffs_kappa_0_a);
    Kokkos::deep_copy(h_coeffs_kappa_a, d_coeffs_kappa_a);
    Kokkos::deep_copy(h_coeffs_kappa_s, d_coeffs_kappa_s);

    printf("# Printing result\n");
    for (int i = 0; i < num_data; i++)
    {
	if (std::is_same_v<BS_REAL, float>)
	{
	    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e "
               "%e\n",
               h_r(i), h_diff_distribution(i), h_coeffs_eta_0(i, id_nue),
               h_coeffs_eta_0(i, id_anue), h_coeffs_eta_0(i, id_nux),
               h_coeffs_eta_0(i, id_anux), h_coeffs_eta(i, id_nue),
               h_coeffs_eta(i, id_anue), h_coeffs_eta(i, id_nux),
               h_coeffs_eta(i, id_anux), h_coeffs_kappa_0_a(i, id_nue),
               h_coeffs_kappa_0_a(i, id_anue), h_coeffs_kappa_0_a(i, id_nux),
               h_coeffs_kappa_0_a(i, id_anux), h_coeffs_kappa_a(i, id_nue),
               h_coeffs_kappa_a(i, id_anue), h_coeffs_kappa_a(i, id_nux),
               h_coeffs_kappa_a(i, id_anux), h_coeffs_kappa_s(i, id_nue),
               h_coeffs_kappa_s(i, id_anue), h_coeffs_kappa_s(i, id_nux),
               h_coeffs_kappa_s(i, id_anux));
       }
	else
	{	
        printf("%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le "
               "%le %le %le %le %le %le "
               "%le\n",
               h_r(i), h_diff_distribution(i), h_coeffs_eta_0(i, id_nue),
               h_coeffs_eta_0(i, id_anue), h_coeffs_eta_0(i, id_nux),
               h_coeffs_eta_0(i, id_anux), h_coeffs_eta(i, id_nue),
               h_coeffs_eta(i, id_anue), h_coeffs_eta(i, id_nux),
               h_coeffs_eta(i, id_anux), h_coeffs_kappa_0_a(i, id_nue),
               h_coeffs_kappa_0_a(i, id_anue), h_coeffs_kappa_0_a(i, id_nux),
               h_coeffs_kappa_0_a(i, id_anux), h_coeffs_kappa_a(i, id_nue),
               h_coeffs_kappa_a(i, id_anue), h_coeffs_kappa_a(i, id_nux),
               h_coeffs_kappa_a(i, id_anux), h_coeffs_kappa_s(i, id_nue),
               h_coeffs_kappa_s(i, id_anue), h_coeffs_kappa_s(i, id_nux),
               h_coeffs_kappa_s(i, id_anux));
       }
    }
}

inline void TestM1OpacitiesSelectedPoints(char filename[200],
                                          OpacityFlags* opacity_flags,
                                          OpacityParams* opacity_pars)
{

    char filepath[300] = {'\0'};
    char filedir[300]  = SOURCE_DIR;
    char outname[200] =
        "/inputs/DD2_M12980-12980_M1_LR/thermo_points_with_neutrinos.txt";

    strcat(filepath, filedir);
    strcat(filepath, outname);

    printf("# Data_directory: %s\n", filepath);

    FILE* fptr;
    fptr = fopen(filepath, "r");
    if (fptr == NULL)
    {
        printf("%s: The file %s does not exist!\n", __FILE__, filepath);
        exit(1);
    }

    // store columns here
    int num_data = 6;
    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_zone("zone", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_rho("rho", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_T("T", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Ye("Ye", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_e("mu_e",
                                                               num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_n("mu_n",
                                                               num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_p("mu_p",
                                                               num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yp("Yp", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yn("Yn", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_dU("dU", num_data);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_nnu(
        "nnu", num_data, total_num_species);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_jnu(
        "Jnu", num_data, total_num_species);
    Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_chinu(
        "chinu", num_data, total_num_species);

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_zone("zone", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_rho("rho", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_T("T", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Ye("Ye", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_e("mu_e", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_n("mu_n", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_p("mu_p", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yp("Yp", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yn("Yn", num_data);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_dU("dU", num_data);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_nnu(
        "nnu", num_data, total_num_species);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_jnu(
        "Jnu", num_data, total_num_species);
    Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_chinu(
        "chinu", num_data, total_num_species);

    // read in the data file
    int i = 0;
    char line[1000];
    while (fgets(line, sizeof(line), fptr) != NULL)
    {
        if (line[0] == '#')
        {
            continue;
        }

	if (std::is_same_v<BS_REAL, float>)
        {
		sscanf(line,
               "%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
               &h_zone(i), &h_rho(i), &h_T(i), &h_Ye(i), &h_Yn(i), &h_Yp(i),
               &h_mu_e(i), &h_mu_n(i), &h_mu_p(i), &h_dU(i), &h_nnu(i, 0),
               &h_nnu(i, 1), &h_nnu(i, 2), &h_jnu(i, 0), &h_jnu(i, 1),
               &h_jnu(i, 2), &h_chinu(i, 0), &h_chinu(i, 1), &h_chinu(i, 2));
	} else {
        sscanf(line,
               "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
               "%lf %lf %lf\n",
               &h_zone(i), &h_rho(i), &h_T(i), &h_Ye(i), &h_Yn(i), &h_Yp(i),
               &h_mu_e(i), &h_mu_n(i), &h_mu_p(i), &h_dU(i), &h_nnu(i, 0),
               &h_nnu(i, 1), &h_nnu(i, 2), &h_jnu(i, 0), &h_jnu(i, 1),
               &h_jnu(i, 2), &h_chinu(i, 0), &h_chinu(i, 1), &h_chinu(i, 2));
	}

        h_nnu(i, 3)   = h_nnu(i, 2);
        h_jnu(i, 3)   = h_jnu(i, 2);
        h_chinu(i, 3) = h_chinu(i, 2);

        i++;
    }

    fclose(fptr);


    printf("# Test for distribution function implementation:\n");

    Kokkos::deep_copy(d_zone, h_zone);
    Kokkos::deep_copy(d_rho, h_rho);
    Kokkos::deep_copy(d_T, h_T);
    Kokkos::deep_copy(d_Ye, h_Ye);
    Kokkos::deep_copy(d_mu_e, h_mu_e);
    Kokkos::deep_copy(d_mu_n, h_mu_n);
    Kokkos::deep_copy(d_mu_p, h_mu_p);
    Kokkos::deep_copy(d_Yp, h_Yp);
    Kokkos::deep_copy(d_Yn, h_Yn);
    Kokkos::deep_copy(d_dU, h_dU);
    Kokkos::deep_copy(d_nnu, h_nnu);
    Kokkos::deep_copy(d_jnu, h_jnu);
    Kokkos::deep_copy(d_chinu, h_chinu);


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

    printf("# Generating quadratures ...\n");
    MyQuadrature my_quad = {.type   = kGauleg,
                            .alpha  = -42.,
                            .dim    = 1,
                            .nx     = 50,
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
    printf("# Quadratures generated.\n");

    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_weights("h_weights",
                                                                  my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_points("h_points",
                                                                 my_quad.nx);

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_weights("weights",
                                                                 my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_points("points",
                                                                my_quad.nx);

    printf("# Copying quadrature data to Kokkos view on host\n");

    h_n_quad(0) = my_quad.nx;

    for (int i = 0; i < my_quad.nx; i++)
    {
        h_weights(i) = my_quad.w[i];
        h_points(i)  = my_quad.points[i];
    }

    printf("# Deep copying quadrature data from host to device\n");

    Kokkos::deep_copy(d_n_quad, h_n_quad);
    Kokkos::deep_copy(d_weights, h_weights);
    Kokkos::deep_copy(d_points, h_points);

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

    printf("# Generated tables:\n");
    printf("# r diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux "
           "j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux "
           "kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue "
           "kappa-s-anue kappa-s-nux kappa-s-anux\n");

    printf("# Entering Kokkos parallel_for loop\n");

    Kokkos::parallel_for(
        "loop_over_ccsn", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data),
        KOKKOS_LAMBDA(const int& i) {
            printf("# Entered in Kokkos parallel_for loop\n");

            GreyOpacityParams my_grey_opacity_params;

            my_grey_opacity_params.opacity_flags = {
                .use_abs_em          = d_opacity_flags(0),
                .use_pair            = d_opacity_flags(1),
                .use_brem            = d_opacity_flags(2),
                .use_inelastic_scatt = d_opacity_flags(3),
                .use_iso             = d_opacity_flags(4)};

            // Opacity parameters (corrections all switched off)
            my_grey_opacity_params.opacity_pars = {
                .use_dU             = d_opacity_pars(0),
                .use_dm_eff         = d_opacity_pars(1),
                .use_WM_ab          = d_opacity_pars(2),
                .use_WM_sc          = d_opacity_pars(3),
                .use_decay          = d_opacity_pars(4),
                .use_BRT_brem       = d_opacity_pars(5),
                .use_NN_medium_corr = d_opacity_pars(6),
                .neglect_blocking   = d_opacity_pars(7)};


            // populate EOS parameters from table
            my_grey_opacity_params.eos_pars.mu_e = d_mu_e(i);
            my_grey_opacity_params.eos_pars.mu_p = d_mu_p(i);
            my_grey_opacity_params.eos_pars.mu_n = d_mu_n(i);
            my_grey_opacity_params.eos_pars.temp = d_T(i);
            my_grey_opacity_params.eos_pars.yp   = d_Yp(i);
            my_grey_opacity_params.eos_pars.yn   = d_Yn(i);
            my_grey_opacity_params.eos_pars.nb   = d_rho(i) / kBS_Mu * 1e-21;
            my_grey_opacity_params.eos_pars.ye   = d_Ye(i);
            my_grey_opacity_params.eos_pars.dU   = d_dU(i);

            // M1 parameters
            for (int idx = 0; idx < total_num_species; ++idx)
            {
                my_grey_opacity_params.m1_pars.chi[idx] = d_chinu(i, idx);
                my_grey_opacity_params.m1_pars.n[idx]   = d_nnu(i, idx) * 1e-21;
                my_grey_opacity_params.m1_pars.J[idx] =
                    d_jnu(i, idx) * kBS_MeV * 1e-21;
            }

            // Distribution parameters
            my_grey_opacity_params.distr_pars =
                CalculateDistrParamsFromM1(&my_grey_opacity_params.m1_pars,
                                           &my_grey_opacity_params.eos_pars);

            printf("# Generating and populating quadrature on GPU\n");

            MyQuadrature gpu_quad;
            gpu_quad.nx = d_n_quad(0);
            for (int idx = 0; idx < gpu_quad.nx; idx++)
            {
                gpu_quad.w[idx]      = d_weights(idx);
                gpu_quad.points[idx] = d_points(idx);
            }

            BS_REAL distr_fermi =
                FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp,
                           my_grey_opacity_params.eos_pars.mu_e -
                               my_grey_opacity_params.eos_pars.mu_n +
                               my_grey_opacity_params.eos_pars.mu_p);
            BS_REAL distr_nuftot =
                TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);

            BS_REAL diff_distr = Kokkos::fabs(distr_fermi - distr_nuftot);

            printf("# Computing M1 coefficients\n");

            M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad,
                                                    &my_grey_opacity_params);

            printf("# Copying M1 coefficients to Kokkos view on device\n");
            d_diff_distribution(i) = diff_distr;

            d_coeffs_eta_0(i, id_nue)  = coeffs.eta_0[id_nue];
            d_coeffs_eta_0(i, id_anue) = coeffs.eta_0[id_anue];
            d_coeffs_eta_0(i, id_nux)  = coeffs.eta_0[id_nux];
            d_coeffs_eta_0(i, id_anux) = coeffs.eta_0[id_anux];

            d_coeffs_eta(i, id_nue)  = coeffs.eta[id_nue];
            d_coeffs_eta(i, id_anue) = coeffs.eta[id_anue];
            d_coeffs_eta(i, id_nux)  = coeffs.eta[id_nux];
            d_coeffs_eta(i, id_anux) = coeffs.eta[id_anux];

            d_coeffs_kappa_0_a(i, id_nue)  = coeffs.kappa_0_a[id_nue];
            d_coeffs_kappa_0_a(i, id_anue) = coeffs.kappa_0_a[id_anue];
            d_coeffs_kappa_0_a(i, id_nux)  = coeffs.kappa_0_a[id_nux];
            d_coeffs_kappa_0_a(i, id_anux) = coeffs.kappa_0_a[id_anux];

            d_coeffs_kappa_a(i, id_nue)  = coeffs.kappa_a[id_nue];
            d_coeffs_kappa_a(i, id_anue) = coeffs.kappa_a[id_anue];
            d_coeffs_kappa_a(i, id_nux)  = coeffs.kappa_a[id_nux];
            d_coeffs_kappa_a(i, id_anux) = coeffs.kappa_a[id_anux];

            d_coeffs_kappa_s(i, id_nue)  = coeffs.kappa_s[id_nue];
            d_coeffs_kappa_s(i, id_anue) = coeffs.kappa_s[id_anue];
            d_coeffs_kappa_s(i, id_nux)  = coeffs.kappa_s[id_nux];
            d_coeffs_kappa_s(i, id_anux) = coeffs.kappa_s[id_anux];
        });

    printf("# Back on CPU, waiting for Kokkos fence\n");

    Kokkos::fence();

    printf("# Deep copying output from device to host\n");

    Kokkos::deep_copy(h_diff_distribution, d_diff_distribution);
    Kokkos::deep_copy(h_coeffs_eta_0, d_coeffs_eta_0);
    Kokkos::deep_copy(h_coeffs_eta, d_coeffs_eta);
    Kokkos::deep_copy(h_coeffs_kappa_0_a, d_coeffs_kappa_0_a);
    Kokkos::deep_copy(h_coeffs_kappa_a, d_coeffs_kappa_a);
    Kokkos::deep_copy(h_coeffs_kappa_s, d_coeffs_kappa_s);

    printf("# Printing result\n");
    for (int i = 0; i < num_data; i++)
    {
	if (std::is_same_v<BS_REAL, float>)
	{
		printf("%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e "
               "%e\n",
               h_zone(i), h_diff_distribution(i), h_coeffs_eta_0(i, id_nue),
               h_coeffs_eta_0(i, id_anue), h_coeffs_eta_0(i, id_nux),
               h_coeffs_eta_0(i, id_anux), h_coeffs_eta(i, id_nue),
               h_coeffs_eta(i, id_anue), h_coeffs_eta(i, id_nux),
               h_coeffs_eta(i, id_anux), h_coeffs_kappa_0_a(i, id_nue),
               h_coeffs_kappa_0_a(i, id_anue), h_coeffs_kappa_0_a(i, id_nux),
               h_coeffs_kappa_0_a(i, id_anux), h_coeffs_kappa_a(i, id_nue),
               h_coeffs_kappa_a(i, id_anue), h_coeffs_kappa_a(i, id_nux),
               h_coeffs_kappa_a(i, id_anux), h_coeffs_kappa_s(i, id_nue),
               h_coeffs_kappa_s(i, id_anue), h_coeffs_kappa_s(i, id_nux),
               h_coeffs_kappa_s(i, id_anux));
	} else {
        printf("%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le "
               "%le %le %le %le %le "
               "%le\n",
               h_zone(i), h_diff_distribution(i), h_coeffs_eta_0(i, id_nue),
               h_coeffs_eta_0(i, id_anue), h_coeffs_eta_0(i, id_nux),
               h_coeffs_eta_0(i, id_anux), h_coeffs_eta(i, id_nue),
               h_coeffs_eta(i, id_anue), h_coeffs_eta(i, id_nux),
               h_coeffs_eta(i, id_anux), h_coeffs_kappa_0_a(i, id_nue),
               h_coeffs_kappa_0_a(i, id_anue), h_coeffs_kappa_0_a(i, id_nux),
               h_coeffs_kappa_0_a(i, id_anux), h_coeffs_kappa_a(i, id_nue),
               h_coeffs_kappa_a(i, id_anue), h_coeffs_kappa_a(i, id_nux),
               h_coeffs_kappa_a(i, id_anux), h_coeffs_kappa_s(i, id_nue),
               h_coeffs_kappa_s(i, id_anue), h_coeffs_kappa_s(i, id_nux),
               h_coeffs_kappa_s(i, id_anux));
	}
    }
}

inline void TestSpectralOpacities(OpacityFlags* opacity_flags,
                                  OpacityParams* opacity_pars)
{

    const int n_bins = 100;

    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_e_bins("e_bins",
                                                                 n_bins);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_nb("nb", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_T("T", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Ye("Ye", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_e("mu_e", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_n("mu_n", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mu_p("mu_p", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yp("Yp", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_Yn("Yn", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mn_eff("mn*", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_mp_eff("mp*", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_dU("dU", 6);

    // Energy bins
    const BS_REAL log_e_min = -1;
    const BS_REAL log_e_max = +3;

    const BS_REAL dE = (log_e_max - log_e_min) / ((BS_REAL)n_bins);

    for (int i = 0; i < n_bins; i++)
    {
        h_e_bins(i) = pow(10., log_e_min + i * dE);
    }

    // First thermodynamic point
    h_nb(0)     = 4.208366627847035e+38 * 1.0E-21; // nm-3
    h_T(0)      = 12.406403541564941;              // MeV
    h_Ye(0)     = 0.07158458232879639;
    h_mu_e(0)   = 187.1814489;   // MeV
    h_mu_n(0)   = 1221.59013681; // MeV
    h_mu_p(0)   = 1011.01797737; // MeV
    h_Yp(0)     = 0.07158458;
    h_Yn(0)     = 0.92841542;
    h_mn_eff(0) = 280.16495513; // MeV
    h_mp_eff(0) = 278.87162217; // MeV
    h_dU(0)     = 18.92714728;  // MeV

    // Second thermodynamic point
    h_nb(1)     = 6.051432235073535e+37 * 1.0E-21; // nm-3
    h_T(1)      = 16.887954711914062;              // MeV
    h_Ye(1)     = 0.05757330730557442;
    h_mu_e(1)   = 82.3367694;   // MeV
    h_mu_n(1)   = 943.76817391; // MeV
    h_mu_p(1)   = 843.50691102; // MeV
    h_Yp(1)     = 0.04352445;
    h_Yn(1)     = 0.91736565;
    h_mn_eff(1) = 733.97246814; // MeV
    h_mp_eff(1) = 732.67912955; // MeV
    h_dU(1)     = 35.09393792;  // MeV

    // Third thermodynamic point
    h_nb(2)     = 6.0146859596912e+36 * 1.0E-21; // nm-3
    h_T(2)      = 8.919591903686523;             // MeV
    h_Ye(2)     = 0.05690543353557587;
    h_mu_e(2)   = 36.56143955;  // MeV
    h_mu_n(2)   = 930.55123963; // MeV
    h_mu_p(2)   = 889.20482938; // MeV
    h_Yp(2)     = 0.02018168;
    h_Yn(2)     = 0.88856839;
    h_mn_eff(2) = 915.90721218; // MeV
    h_mp_eff(2) = 914.61387308; // MeV
    h_dU(2)     = 4.98247735;   // MeV

    // Fourth thermodynamic point
    h_nb(3)     = 6.128721244299242e+35 * 1.0E-21; // nm-3
    h_T(3)      = 6.698440074920654;               // MeV
    h_Ye(3)     = 0.11529815942049026;
    h_mu_e(3)   = 19.49956921;  // MeV
    h_mu_n(3)   = 920.74113958; // MeV
    h_mu_p(3)   = 903.09559414; // MeV
    h_Yp(3)     = 0.08023092;
    h_Yn(3)     = 0.8442462;
    h_mn_eff(3) = 937.06605368; // MeV
    h_mp_eff(3) = 935.77271381; // MeV
    h_dU(3)     = 0.46637764;   // MeV

    // Fifth thermodynamic point
    h_nb(4)     = 6.07728949804485e+34 * 1.0E-21; // nm-3
    h_T(4)      = 3.586622953414917;              // MeV
    h_Ye(4)     = 0.13497862219810486;
    h_mu_e(4)   = 8.98668385;   // MeV
    h_mu_n(4)   = 924.65692374; // MeV
    h_mu_p(4)   = 916.18996071; // MeV
    h_Yp(4)     = 0.11611014;
    h_Yn(4)     = 0.84457226;
    h_mn_eff(4) = 939.30610577; // MeV
    h_mp_eff(4) = 938.01276696; // MeV
    h_dU(4)     = 0.0454689;    // MeV

    // Sixth thermodynamic point
    h_nb(5)     = 6.192921723010272e+33 * 1.0E-21; // nm-3
    h_T(5)      = 2.144158124923706;               // MeV
    h_Ye(5)     = 0.18740089237689972;
    h_mu_e(5)   = 4.20812739;   // MeV
    h_mu_n(5)   = 927.30360027; // MeV
    h_mu_p(5)   = 922.72216948; // MeV
    h_Yp(5)     = 0.17166997;
    h_Yn(5)     = 0.79639333;
    h_mn_eff(5) = 939.56535; // MeV
    h_mp_eff(5) = 938.27201; // MeV
    h_dU(5)     = 0.0007919; // MeV

    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_e_bins("e_bins",
                                                                n_bins);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_nb("nb", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_T("T", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Ye("Ye", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_e("mu_e", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_n("mu_n", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mu_p("mu_p", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yp("Yp", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_Yn("Yn", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mn_eff("mn*", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_mp_eff("mp*", 6);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_dU("dU", 6);

    Kokkos::deep_copy(d_e_bins, h_e_bins);
    Kokkos::deep_copy(d_nb, h_nb);
    Kokkos::deep_copy(d_T, h_T);
    Kokkos::deep_copy(d_Ye, h_Ye);
    Kokkos::deep_copy(d_mu_e, h_mu_e);
    Kokkos::deep_copy(d_mu_n, h_mu_n);
    Kokkos::deep_copy(d_mu_p, h_mu_p);
    Kokkos::deep_copy(d_Yp, h_Yp);
    Kokkos::deep_copy(d_Yn, h_Yn);
    Kokkos::deep_copy(d_mn_eff, h_mn_eff);
    Kokkos::deep_copy(d_mp_eff, h_mp_eff);
    Kokkos::deep_copy(d_dU, h_dU);

    Kokkos::View<BS_REAL***, LayoutWrapper, DevMemSpace> d_j("j", 6, n_bins,
                                                             total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, DevMemSpace> d_j_s(
        "j_s", 6, n_bins, total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, DevMemSpace> d_kappa(
        "kappa", 6, n_bins, total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, DevMemSpace> d_kappa_s(
        "kappa_s", 6, n_bins, total_num_species);

    Kokkos::View<BS_REAL***, LayoutWrapper, HostMemSpace> h_j(
        "j", 6, n_bins, total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, HostMemSpace> h_j_s(
        "j_s", 6, n_bins, total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, HostMemSpace> h_kappa(
        "kappa", 6, n_bins, total_num_species);
    Kokkos::View<BS_REAL***, LayoutWrapper, HostMemSpace> h_kappa_s(
        "kappa_s", 6, n_bins, total_num_species);

    printf("# Generating quadratures ...\n");
    MyQuadrature my_quad = {.type   = kGauleg,
                            .alpha  = -42.,
                            .dim    = 1,
                            .nx     = 10,
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
    printf("# Quadratures generated.\n");

    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_weights("h_weights",
                                                                  my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, HostMemSpace> h_points("h_points",
                                                                 my_quad.nx);

    Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_n_quad("n_quad", 1);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_weights("weights",
                                                                 my_quad.nx);
    Kokkos::View<BS_REAL*, LayoutWrapper, DevMemSpace> d_points("points",
                                                                my_quad.nx);

    printf("# Copying quadrature data to Kokkos view on host\n");

    h_n_quad(0) = my_quad.nx;

    for (int i = 0; i < my_quad.nx; i++)
    {
        h_weights(i) = my_quad.w[i];
        h_points(i)  = my_quad.points[i];
    }

    printf("# Deep copying quadrature data from host to device\n");

    Kokkos::deep_copy(d_n_quad, h_n_quad);
    Kokkos::deep_copy(d_weights, h_weights);
    Kokkos::deep_copy(d_points, h_points);

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

    printf("# Generated tables:\n");
    printf("# Energy j-nue j-anue j0-nux j-anux j-s-nue j-s-anue j-s-nux "
           "j-s-anux kappa-nue kappa-anue kappa-nux kappa-anux "
           "kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");

    printf("# Entering Kokkos parallel_for loop\n");

    Kokkos::parallel_for(
        "loop_over_energy_bins",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, n_bins),
        KOKKOS_LAMBDA(const int& i) {
            printf("# Entered in Kokkos parallel_for loop\n");

            GreyOpacityParams my_grey_opacity_params;

            my_grey_opacity_params.opacity_flags = {
                .use_abs_em          = d_opacity_flags(0),
                .use_pair            = d_opacity_flags(1),
                .use_brem            = d_opacity_flags(2),
                .use_inelastic_scatt = d_opacity_flags(3),
                .use_iso             = d_opacity_flags(4)};


            // Opacity parameters (corrections all switched off)
            my_grey_opacity_params.opacity_pars = {
                .use_dU             = d_opacity_pars(0),
                .use_dm_eff         = d_opacity_pars(1),
                .use_WM_ab          = d_opacity_pars(2),
                .use_WM_sc          = d_opacity_pars(3),
                .use_decay          = d_opacity_pars(4),
                .use_BRT_brem       = d_opacity_pars(5),
                .use_NN_medium_corr = d_opacity_pars(6),
                .neglect_blocking   = d_opacity_pars(7)};

            printf("# Generating and populating quadrature on GPU\n");
            MyQuadrature gpu_quad;
            gpu_quad.nx = d_n_quad(0);
            for (int idx = 0; idx < gpu_quad.nx; idx++)
            {
                gpu_quad.w[idx]      = d_weights(idx);
                gpu_quad.points[idx] = d_points(idx);
            }

            for (int j = 0; j < 6; j++)
            {
                // populate EOS parameters from table
                my_grey_opacity_params.eos_pars.mu_e = d_mu_e(j);
                my_grey_opacity_params.eos_pars.mu_p = d_mu_p(j);
                my_grey_opacity_params.eos_pars.mu_n = d_mu_n(j);
                my_grey_opacity_params.eos_pars.temp = d_T(j);
                my_grey_opacity_params.eos_pars.yp   = d_Yp(j);
                my_grey_opacity_params.eos_pars.yn   = d_Yn(j);
                my_grey_opacity_params.eos_pars.nb   = d_nb(j);
                my_grey_opacity_params.eos_pars.dU   = d_dU(j);
                my_grey_opacity_params.eos_pars.dm_eff =
                    d_mn_eff(j) - d_mp_eff(j);

                // Distribution parameters
                my_grey_opacity_params.distr_pars = NuEquilibriumParams(
                    &my_grey_opacity_params
                         .eos_pars); // consider neutrino distribution function
                                     // at equilibrium

                // M1 parameters
                ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                                     &my_grey_opacity_params.distr_pars,
                                     &my_grey_opacity_params.m1_pars);

                for (int idx = 0; idx < total_num_species; ++idx)
                {
                    my_grey_opacity_params.m1_pars.chi[idx] = 1. / 3.;
                    my_grey_opacity_params.m1_pars.J[idx] =
                        my_grey_opacity_params.m1_pars.J[idx] * kBS_MeV;
                }

                printf("# Computing spectral rates\n");

                SpectralOpacities rates =
                    ComputeSpectralOpacitiesNotStimulatedAbs(
                        d_e_bins(i), &gpu_quad, &my_grey_opacity_params);

                printf("# Copying spectral rates to Kokkos view on device\n");

                d_j(j, i, id_nue)  = rates.j[id_nue];
                d_j(j, i, id_anue) = rates.j[id_anue];
                d_j(j, i, id_nux)  = rates.j[id_nux];
                d_j(j, i, id_anux) = rates.j[id_anux];

                d_j_s(j, i, id_nue)  = rates.j_s[id_nue];
                d_j_s(j, i, id_anue) = rates.j_s[id_anue];
                d_j_s(j, i, id_nux)  = rates.j_s[id_nux];
                d_j_s(j, i, id_anux) = rates.j_s[id_anux];

                d_kappa(j, i, id_nue)  = rates.kappa[id_nue];
                d_kappa(j, i, id_anue) = rates.kappa[id_anue];
                d_kappa(j, i, id_nux)  = rates.kappa[id_nux];
                d_kappa(j, i, id_anux) = rates.kappa[id_anux];

                d_kappa_s(j, i, id_nue)  = rates.kappa_s[id_nue];
                d_kappa_s(j, i, id_anue) = rates.kappa_s[id_anue];
                d_kappa_s(j, i, id_nux)  = rates.kappa_s[id_nux];
                d_kappa_s(j, i, id_anux) = rates.kappa_s[id_anux];
            }
        });

    printf("# Back on CPU, waiting for Kokkos fence\n");

    Kokkos::fence();

    printf("# Deep copying output from device to host\n");

    Kokkos::deep_copy(h_j, d_j);
    Kokkos::deep_copy(h_j_s, d_j_s);
    Kokkos::deep_copy(h_kappa, d_kappa);
    Kokkos::deep_copy(h_kappa_s, d_kappa_s);

    printf("# Printing result\n");
    for (int i = 0; i < n_bins; i++)
    {
	if (std::is_same_v<BS_REAL, float>)
	{
		printf("%.8e ", h_e_bins(i));
	} else {
        printf("%.8le ", h_e_bins(i));
	}

        for (int j = 0; j < 6; j++)
        {
	if (std::is_same_v<BS_REAL, float>)
	{
	    	printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",
                   h_j(j, i, id_nue), h_j(j, i, id_anue), h_j(j, i, id_nux),
                   h_j(j, i, id_anux), h_j_s(j, i, id_nue),
                   h_j_s(j, i, id_anue), h_j_s(j, i, id_nux),
                   h_j_s(j, i, id_anux), h_kappa(j, i, id_nue),
                   h_kappa(j, i, id_anue), h_kappa(j, i, id_nux),
                   h_kappa(j, i, id_anux), h_kappa_s(j, i, id_nue),
                   h_kappa_s(j, i, id_anue), h_kappa_s(j, i, id_nux),
                   h_kappa_s(j, i, id_anux));
	} else {
            printf("%le %le %le %le %le %le %le %le %le %le %le %le %le %le "
                   "%le %le ",
                   h_j(j, i, id_nue), h_j(j, i, id_anue), h_j(j, i, id_nux),
                   h_j(j, i, id_anux), h_j_s(j, i, id_nue),
                   h_j_s(j, i, id_anue), h_j_s(j, i, id_nux),
                   h_j_s(j, i, id_anux), h_kappa(j, i, id_nue),
                   h_kappa(j, i, id_anue), h_kappa(j, i, id_nux),
                   h_kappa(j, i, id_anux), h_kappa_s(j, i, id_nue),
                   h_kappa_s(j, i, id_anue), h_kappa_s(j, i, id_nux),
                   h_kappa_s(j, i, id_anux));
	   }
        }
        printf("\n");
    }
}

#endif // BNS_NURATES_TESTS_TESTS_HPP_
