//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_m1_opacities.c
//  \brief Generate a table for pair opacities

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#include <Kokkos_Core.hpp>

#include "../tests.hpp"
#include "opacities.hpp"
#include "m1_opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

void TestM1Opacities(char filename[200], OpacityFlags *opacity_flags, OpacityParams *opacity_pars) {

  using DevExeSpace = Kokkos::DefaultExecutionSpace;
  using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
  using HostMemSpace = Kokkos::HostSpace;
  using LayoutWrapper = Kokkos::LayoutRight;                // increments last index fastest views defined like

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data_directory: %s\n", filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  // store columns here
  int num_data = 102;
  double e_nu;
  Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_zone("zone", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_r("r", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_rho("rho", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_T("T", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_Ye("Ye", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_e("mu_e", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_hat("mu_hat", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_Yh("Yh", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_Ya("Ya", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_Yp("Yp", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_Yn("Yn", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_em_nue("em_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_l_nue_inv("l_nue_inv", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_em_anue("em_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_l_anue_inv("l_anue_inv", num_data);

  Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_zone("zone", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_r("r", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_rho("rho", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_T("T", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_Ye("Ye", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_mu_e("mu_e", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_mu_hat("mu_hat", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_Yh("Yh", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_Ya("Ya", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_Yp("Yp", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_Yn("Yn", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_em_nue("em_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_l_nue_inv("l_nue_inv", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_em_anue("em_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_l_anue_inv("l_anue_inv", num_data);

  // read in the data file
  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[1] == '#' && i == 0) {
      sscanf(line + 14, "%lf\n", &e_nu);
      continue;
    } else if (line[1] == '#' && i != 0) {
      continue;
    }

    sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &h_zone(i), &h_r(i), &h_rho(i), &h_T(i), &h_Ye(i), &h_mu_e(i), &h_mu_hat(i), &h_Yh[i], &h_Ya[i], &h_Yp(i), &h_Yn(i), &h_em_nue(i), &h_l_nue_inv(i), &h_em_anue(i), &h_l_anue_inv(i));
    i++;
  }

  fclose(fptr);
 
  printf("Test for distribution function implementation:\n");

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

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_diff_distribution("diff_distribution", num_data);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> d_coeffs_eta_0("coeffs_eta_0", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> d_coeffs_eta("coeffs_eta", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_0_a("coeffs_kappa_0_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_a("coeffs_kappa_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> d_coeffs_kappa_s("coeffs_kappa_s", num_data, 4);

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_diff_distribution("h_diff_distribution", num_data);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_eta_0("h_coeffs_eta_0", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_eta("h_coeffs_eta", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_0_a("h_coeffs_kappa_0_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_a("h_coeffs_kappa_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_s("h_coeffs_kappa_s", num_data, 4);
 
  printf("Generating quadratures ...\n");
  MyQuadrature my_quad = {.type = kGauleg, .alpha = -42., .dim = 1, .nx = 10, .ny = 1, .nz = 1, .x1 = 0., .x2 = 1., .y1 = -42., .y2 = -42., .z1 = -42., .z2 = -42., .points = {0}, .w = {0}};
  GaussLegendre(&my_quad);
  printf("Quadratures generated.\n");
   
  Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_n_quad("n_quad", 1);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_weights("h_weights", my_quad.nx);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_points("h_points", my_quad.nx);
 
  Kokkos::View<int*, LayoutWrapper, DevMemSpace> d_n_quad("n_quad", 1);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_weights("weights", my_quad.nx);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> d_points("points", my_quad.nx);
  
  printf("Copying quadrature data to Kokkos view on host\n");

  h_n_quad(0) = my_quad.nx;
 
  for (int i = 0; i < my_quad.nx; i ++) {
     h_weights(i) = my_quad.w[i];
     h_points(i) = my_quad.points[i];
  }
  
  printf("Deep copying quadrature data from host to device\n");
 
  Kokkos::deep_copy(d_n_quad, h_n_quad);
  Kokkos::deep_copy(d_weights, h_weights);
  Kokkos::deep_copy(d_points, h_points); 

  printf("Generated tables:\n");
  printf("r diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");

  printf("Entering Kokkos parallel_for loop\n");
  
  Kokkos::parallel_for("loop_over_ccsn", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data),
  KOKKOS_LAMBDA(const int &i) {
    printf("Entered in Kokkos parallel_for loop\n");

    GreyOpacityParams my_grey_opacity_params;
    
    my_grey_opacity_params.opacity_flags = {.use_abs_em          = 1,
                                            .use_pair            = 1,
                                            .use_brem            = 1,
                                            .use_inelastic_scatt = 1,
                                            .use_iso             = 1};
    

    // populate EOS parameters from table
    my_grey_opacity_params.eos_pars.mu_e = d_mu_e(i);
    my_grey_opacity_params.eos_pars.mu_p = 0.;
    my_grey_opacity_params.eos_pars.mu_n = d_mu_hat(i) + kQ;
    my_grey_opacity_params.eos_pars.temp = d_T(i);
    my_grey_opacity_params.eos_pars.yp = d_Yp(i);
    my_grey_opacity_params.eos_pars.yn = d_Yn(i);
    my_grey_opacity_params.eos_pars.nb = d_rho(i) / kMu;


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
    my_grey_opacity_params.distr_pars = NuEquilibriumParams(&my_grey_opacity_params.eos_pars); // consider neutrino distribution function at equilibrium

    // M1 parameters
    ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars, &my_grey_opacity_params.distr_pars, &my_grey_opacity_params.m1_pars);
  
    for (int idx = 0; idx < total_num_species; idx++) {
      my_grey_opacity_params.m1_pars.chi[idx] = 1. / 3.;
      my_grey_opacity_params.m1_pars.J[idx] = my_grey_opacity_params.m1_pars.J[idx] * kMeV;
    }

    printf("Generating and populating quadrature on GPU\n");
 
    MyQuadrature gpu_quad;
    gpu_quad.nx = d_n_quad(0);
    for (int idx = 0; idx < gpu_quad.nx; idx ++) {
       gpu_quad.w[idx] = d_weights(idx);
       gpu_quad.points[idx] = d_points(idx);
    }
 
    double distr_fermi =
        FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp, my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p);
    double distr_nuftot = TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);
    double diff_distr = Kokkos::fabs(distr_fermi - distr_nuftot);
 
    printf("Computing M1 coefficients\n");
  
    //const int n = gpu_quad.nx;
    //double t = 1.5 * d_T(i);
    //double nu_array[n_max];
    //for (int i = 0; i < n; i++)
    //{
        //nu_array[i]     = t * gpu_quad.points[i];
        //nu_array[n + i] = t / gpu_quad.points[i];
    //}
    //OpacityParams opacity_pars = grey_pars->opacity_pars;
    /* 
    M1MatrixKokkos2D pair = {0};
    if (my_grey_opacity_params.opacity_flags.use_pair == 1)
    {
        PairKernelsTable(2 * n, nu_array, &my_grey_opacity_params, &pair);
    }

    
    M1MatrixKokkos2D brem = {0};
    if (my_grey_opacity_params.opacity_flags.use_brem == 1)
    {
        if (my_grey_opacity_params.opacity_pars.use_BRT_brem == true)
        {
            BremKernelsTableBRT06(2 * n, nu_array, &my_grey_opacity_params, &brem);
        }
        else
        {
            BremKernelsTable(2 * n, nu_array, &my_grey_opacity_params, &brem);
	}
    }
   

    
    M1MatrixKokkos2D inel = {0};
    if (my_grey_opacity_params.opacity_flags.use_inelastic_scatt == 1)
    {
        InelasticKernelsTable(2 * n, nu_array, &my_grey_opacity_params, &inel);
    }
    

    printf("%.5e\n", pair.m1_mat_em[0][10][10]);
    printf("%.5e\n", brem.m1_mat_em[0][10][10]);
    printf("%.5e\n", inel.m1_mat_em[0][10][10]);
    */

    /*
    M1MatrixKokkos2D out_n = {0}, out_j = {0};

    ComputeM1DoubleIntegrandTest(&gpu_quad, &my_grey_opacity_params, 0.5 * 4.364 * d_T(i), &out_n, &out_j);

    printf("%.5e\n", out_n.m1_mat_em[0][10][10]);
    printf("%.5e\n", out_n.m1_mat_em[0][10][10]);
    printf("%.5e\n", out_j.m1_mat_em[0][10][10]);
    printf("%.5e\n", out_j.m1_mat_em[0][10][10]);
    */

    M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad, &my_grey_opacity_params);

    printf("Copying M1 coefficients to Kokkos view on device\n");
    d_diff_distribution(i) = diff_distr;

    d_coeffs_eta_0(i,id_nue) = coeffs.eta_0[id_nue];
    d_coeffs_eta_0(i,id_anue) = coeffs.eta_0[id_anue];
    d_coeffs_eta_0(i,id_nux) = coeffs.eta_0[id_nux];
    d_coeffs_eta_0(i,id_anux) = coeffs.eta_0[id_anux];

    d_coeffs_eta(i,id_nue) = coeffs.eta[id_nue];
    d_coeffs_eta(i,id_anue) = coeffs.eta[id_anue];
    d_coeffs_eta(i,id_nux) = coeffs.eta[id_nux];
    d_coeffs_eta(i,id_anux) = coeffs.eta[id_anux];

    d_coeffs_kappa_0_a(i,id_nue) = coeffs.kappa_0_a[id_nue];
    d_coeffs_kappa_0_a(i,id_anue) = coeffs.kappa_0_a[id_anue];
    d_coeffs_kappa_0_a(i,id_nux) = coeffs.kappa_0_a[id_nux];
    d_coeffs_kappa_0_a(i,id_anux) = coeffs.kappa_0_a[id_anux];

    d_coeffs_kappa_a(i,id_nue) = coeffs.kappa_a[id_nue];
    d_coeffs_kappa_a(i,id_anue) = coeffs.kappa_a[id_anue];
    d_coeffs_kappa_a(i,id_nux) = coeffs.kappa_a[id_nux];
    d_coeffs_kappa_a(i,id_anux) = coeffs.kappa_a[id_anux];

    d_coeffs_kappa_s(i,id_nue) = coeffs.kappa_s[id_nue];
    d_coeffs_kappa_s(i,id_anue) = coeffs.kappa_s[id_anue];
    d_coeffs_kappa_s(i,id_nux) = coeffs.kappa_s[id_nux];
    d_coeffs_kappa_s(i,id_anux) = coeffs.kappa_s[id_anux];
  });
    
  printf("Back on CPU, waiting for Kokkos fence\n");
 
  Kokkos::fence();
 
  printf("Deep copying output from device to host\n");
  
  Kokkos::deep_copy(h_diff_distribution, d_diff_distribution);
  Kokkos::deep_copy(h_coeffs_eta_0, d_coeffs_eta_0);
  Kokkos::deep_copy(h_coeffs_eta, d_coeffs_eta);
  Kokkos::deep_copy(h_coeffs_kappa_0_a, d_coeffs_kappa_0_a);
  Kokkos::deep_copy(h_coeffs_kappa_a, d_coeffs_kappa_a);
  Kokkos::deep_copy(h_coeffs_kappa_s, d_coeffs_kappa_s);

  printf("Printing result\n");
  for (int i = 0; i < num_data; i++) {
    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            h_r(i), h_diff_distribution(i),
            h_coeffs_eta_0(i,id_nue), h_coeffs_eta_0(i,id_anue), h_coeffs_eta_0(i,id_nux), h_coeffs_eta_0(i,id_anux),
            h_coeffs_eta(i,id_nue), h_coeffs_eta(i,id_anue), h_coeffs_eta(i,id_nux), h_coeffs_eta(i,id_anux),
            h_coeffs_kappa_0_a(i,id_nue), h_coeffs_kappa_0_a(i,id_anue), h_coeffs_kappa_0_a(i,id_nux), h_coeffs_kappa_0_a(i,id_anux),
            h_coeffs_kappa_a(i,id_nue), h_coeffs_kappa_a(i,id_anue), h_coeffs_kappa_a(i,id_nux), h_coeffs_kappa_a(i,id_anux),
            h_coeffs_kappa_s(i,id_nue), h_coeffs_kappa_s(i,id_anue), h_coeffs_kappa_s(i,id_nux), h_coeffs_kappa_s(i,id_anux));
  }

  return;
}

