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
#include "../../include/integration.hpp"
#include "../../include/distribution.hpp"
#include "../../include/constants.hpp"
#include "../../include/functions.hpp"

void TestM1Opacities(char filename[200], OpacityFlags *opacity_flags, OpacityParams *opacity_pars) {

  Kokkos::initialize();

  using DevExeSpace = Kokkos::DefaultExecutionSpace;
  using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
  using HostMemSpace = Kokkos::HostSpace;
  using LayoutWrapper = Kokkos::LayoutRight;                // increments last index fastest views defined like

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/inputs/nurates_CCSN/nurates_1.008E+01.txt";

  char data_fileout[200] = "/tests/tests_opacities_m1/output/";
  char data_filepath[300] = {'\0'};

  strcat(filepath, filedir);
  strcat(filepath, outname);
  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);
  strcat(data_filepath, filename);

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

  Kokkos::View<int*, LayoutWrapper, DevMemSpace> zone("zone", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> r("r", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> rho("rho", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> T("T", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> Ye("Ye", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_e("mu_e", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_hat("mu_hat", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> Yh("Yh", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> Ya("Ya", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> Yp("Yp", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> Yn("Yn", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> em_nue("em_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> l_nue_inv("l_nue_inv", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> em_anue("em_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> l_anue_inv("l_anue_inv", num_data);

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

  // print input data
  if (false) {
    printf("\n");
    printf("Printing out input table ... %s\n", outname);
    printf("\n");

    printf("Input data:\n");
    printf("E_nu: %lf\n", e_nu);
    for (int i = 0; i < num_data; i++) {
      printf("%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             h_zone(i), h_r(i), h_rho(i), h_T(i), h_Ye(i), h_mu_e(i), h_mu_hat(i), h_Yh[i], h_Ya[i], h_Yp(i), h_Yn(i), h_em_nue(i), h_l_nue_inv(i), h_em_anue(i), h_l_anue_inv(i));
    }

    printf("\n");
    printf("End of input table.\n");
    printf("\n");
  }

  printf("Test for distribution function implementation:\n");


  Kokkos::deep_copy(zone, h_zone);
  Kokkos::deep_copy(r, h_r);
  Kokkos::deep_copy(rho, h_rho);
  Kokkos::deep_copy(T, h_T);
  Kokkos::deep_copy(Ye, h_Ye);
  Kokkos::deep_copy(mu_e, h_mu_e);
  Kokkos::deep_copy(mu_hat, h_mu_hat);
  Kokkos::deep_copy(Yh, h_Yh);
  Kokkos::deep_copy(Ya, h_Ya);
  Kokkos::deep_copy(Yp, h_Yp);
  Kokkos::deep_copy(Yn, h_Yn);
  Kokkos::deep_copy(em_nue, h_em_nue);
  Kokkos::deep_copy(em_anue, h_em_anue);
  Kokkos::deep_copy(l_nue_inv, h_l_nue_inv);
  Kokkos::deep_copy(l_anue_inv, h_l_anue_inv);

  FILE *file;
  file = fopen(data_filepath, "w+");

  printf("Generating quadratures ...\n");
  MyQuadrature my_quadrature_1d = {.type = kGauleg, .alpha = -42., .dim = 1, .nx = 10, .ny = 1, .nz = 1, .x1 = 0., .x2 = 1., .y1 = -42., .y2 = -42., .z1 = -42., .z2 = -42., .points = NULL, .w = NULL};
  GaussLegendreMultiD(&my_quadrature_1d);
  MyQuadrature my_quadrature_2d = {.type = kGauleg, .alpha = -42., .dim = 1, .nx = 10, .ny = 1, .nz = 1, .x1 = 0., .x2 = 1., .y1 = -42., .y2 = -42., .z1 = -42., .z2 = -42., .points = NULL, .w = NULL};
  GaussLegendreMultiD(&my_quadrature_2d);
  printf("Quadratures generated.\n");

  printf("\n");

  /*  
  Kokkos::View<int> n_quad("n_quad");
  Kokkos::View<int>::HostMirror h_n_quad = Kokkos::create_mirror_view(n_quad);
  
  h_n_quad() = my_quadrature_1d.nx;

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_weights("h_weights", my_quadrature_1d.nx);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_points("h_points", my_quadrature_1d.nx);

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> weights("weights", my_quadrature_1d.nx);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> points("points", my_quadrature_1d.nx);

  for (int i = 0; i < my_quadrature_1d.nx; i ++) {
    h_weights(i) = my_quadrature_1d.w[i];
    h_points(i) = my_quadrature_1d.points[i];
  }

  Kokkos::deep_copy(n_quad, h_n_quad);
  Kokkos::deep_copy(weights, h_weights);
  Kokkos::deep_copy(points, h_points);
  */

  printf("Generated tables:\n");
  printf("r diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");

  fprintf(file, "#  diff_distr j0-nue j0-anue j0-nux j0-anux j-nue j-anue j-nux j-anux kappa0-a-nue kappa0-a-anue kappa0-a-nux kappa0-a-anux kappa-a-nue kappa-a-anue kappa-a-nux kappa-a-anux kappa-s-nue kappa-s-anue kappa-s-nux kappa-s-anux\n");
  // printf("r: [cm], diff_distr [should be zero, checks Fermi Dirac with NuFTotal], j-nue, ");

  clock_t start, end;
  double cpu_time_used;

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> diff_distribution("diff_distribution", num_data);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> coeffs_eta_0("coeffs_eta_0", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> coeffs_eta("coeffs_eta", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> coeffs_kappa_0_a("coeffs_kappa_0_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> coeffs_kappa_a("coeffs_kappa_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, DevMemSpace> coeffs_kappa_s("coeffs_kappa_s", num_data, 4);

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_diff_distribution("h_diff_distribution", num_data);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_eta_0("h_coeffs_eta_0", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_eta("h_coeffs_eta", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_0_a("h_coeffs_kappa_0_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_a("h_coeffs_kappa_a", num_data, 4);
  Kokkos::View<double**, LayoutWrapper, HostMemSpace> h_coeffs_kappa_s("h_coeffs_kappa_s", num_data, 4);

  start = clock();
 
  Kokkos::parallel_for("loop_over_ccsn", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data),
  KOKKOS_LAMBDA(const int &i) {
    GreyOpacityParams my_grey_opacity_params;
    //my_grey_opacity_params.opacity_flags = *opacity_flags;
    my_grey_opacity_params.opacity_flags = {.use_abs_em          = 1,
                                            .use_pair            = 1,
                                            .use_brem            = 1,
                                            .use_inelastic_scatt = 1,
                                            .use_iso             = 1};

    /*
    MyQuadrature my_quad;
    my_quad.nx = n_quad();
    for (int idx = 0; idx < my_quad.nx; idx ++) {
      my_quad.w[idx] = weights(idx);
      my_quad.points[idx] = points(idx);
    }
    */

    // populate EOS parameters from table
    my_grey_opacity_params.eos_pars.mu_e = mu_e(i);
    my_grey_opacity_params.eos_pars.mu_p = 0.;
    my_grey_opacity_params.eos_pars.mu_n = mu_hat(i) + kQ;
    my_grey_opacity_params.eos_pars.temp = T(i);
    my_grey_opacity_params.eos_pars.yp = Yp(i);
    my_grey_opacity_params.eos_pars.yn = Yn(i);
    my_grey_opacity_params.eos_pars.nb = rho(i) / kMu;

    // Opacity parameters (corrections all switched off)
    //my_grey_opacity_params.opacity_pars = *opacity_pars;
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

    double distr_fermi =
        FermiDistr(123.4, my_grey_opacity_params.eos_pars.temp, my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p);
    double distr_nuftot = TotalNuF(123.4, &my_grey_opacity_params.distr_pars, 0);
    double diff_distr = Kokkos::fabs(distr_fermi - distr_nuftot);

    M1Opacities coeffs = ComputeM1Opacities(&my_quadrature_1d, &my_quadrature_2d, &my_grey_opacity_params);
    //M1Opacities coeffs = ComputeM1Opacities(&my_quad, &my_quad, &my_grey_opacity_params);
   
    /*printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r(i), diff_distr,
            coeffs.eta_0[id_nue], coeffs.eta_0[id_anue], coeffs.eta_0[id_nux], coeffs.eta_0[id_anux],
            coeffs.eta[id_nue], coeffs.eta[id_anue], coeffs.eta[id_nux], coeffs.eta[id_anux],
            coeffs.kappa_0_a[id_nue], coeffs.kappa_0_a[id_anue], coeffs.kappa_0_a[id_nux], coeffs.kappa_0_a[id_anux],
            coeffs.kappa_a[id_nue], coeffs.kappa_a[id_anue], coeffs.kappa_a[id_nux], coeffs.kappa_a[id_anux],
            coeffs.kappa_s[id_nue], coeffs.kappa_s[id_anue], coeffs.kappa_s[id_nux], coeffs.kappa_s[id_anux]);
    */
    /*fprintf(file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            r(i), diff_distr,
            coeffs.eta_0[id_nue], coeffs.eta_0[id_anue], coeffs.eta_0[id_nux], coeffs.eta_0[id_anux],
            coeffs.eta[id_nue], coeffs.eta[id_anue], coeffs.eta[id_nux], coeffs.eta[id_anux],
            coeffs.kappa_0_a[id_nue], coeffs.kappa_0_a[id_anue], coeffs.kappa_0_a[id_nux], coeffs.kappa_0_a[id_anux],
            coeffs.kappa_a[id_nue], coeffs.kappa_a[id_anue], coeffs.kappa_a[id_nux], coeffs.kappa_a[id_anux],
            coeffs.kappa_s[id_nue], coeffs.kappa_s[id_anue], coeffs.kappa_s[id_nux], coeffs.kappa_s[id_anux]); */

    diff_distribution(i) = diff_distr;

    coeffs_eta_0(i,id_nue) = coeffs.eta_0[id_nue];
    coeffs_eta_0(i,id_anue) = coeffs.eta_0[id_anue];
    coeffs_eta_0(i,id_nux) = coeffs.eta_0[id_nux];
    coeffs_eta_0(i,id_anux) = coeffs.eta_0[id_anux];

    coeffs_eta(i,id_nue) = coeffs.eta[id_nue];
    coeffs_eta(i,id_anue) = coeffs.eta[id_anue];
    coeffs_eta(i,id_nux) = coeffs.eta[id_nux];
    coeffs_eta(i,id_anux) = coeffs.eta[id_anux];

    coeffs_kappa_0_a(i,id_nue) = coeffs.kappa_0_a[id_nue];
    coeffs_kappa_0_a(i,id_anue) = coeffs.kappa_0_a[id_anue];
    coeffs_kappa_0_a(i,id_nux) = coeffs.kappa_0_a[id_nux];
    coeffs_kappa_0_a(i,id_anux) = coeffs.kappa_0_a[id_anux];

    coeffs_kappa_a(i,id_nue) = coeffs.kappa_a[id_nue];
    coeffs_kappa_a(i,id_anue) = coeffs.kappa_a[id_anue];
    coeffs_kappa_a(i,id_nux) = coeffs.kappa_a[id_nux];
    coeffs_kappa_a(i,id_anux) = coeffs.kappa_a[id_anux];

    coeffs_kappa_s(i,id_nue) = coeffs.kappa_s[id_nue];
    coeffs_kappa_s(i,id_anue) = coeffs.kappa_s[id_anue];
    coeffs_kappa_s(i,id_nux) = coeffs.kappa_s[id_nux];
    coeffs_kappa_s(i,id_anux) = coeffs.kappa_s[id_anux];
  
  });
  


  Kokkos::deep_copy(h_diff_distribution, diff_distribution);
  Kokkos::deep_copy(h_coeffs_eta_0, coeffs_eta_0);
  Kokkos::deep_copy(h_coeffs_eta, coeffs_eta);
  Kokkos::deep_copy(h_coeffs_kappa_0_a, coeffs_kappa_0_a);
  Kokkos::deep_copy(h_coeffs_kappa_a, coeffs_kappa_a);
  Kokkos::deep_copy(h_coeffs_kappa_s, coeffs_kappa_s);

  for (int i = 0; i < num_data; i++) {
    printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            h_r(i), h_diff_distribution(i),
            h_coeffs_eta_0(i,id_nue), h_coeffs_eta_0(i,id_anue), h_coeffs_eta_0(i,id_nux), h_coeffs_eta_0(i,id_anux),
            h_coeffs_eta(i,id_nue), h_coeffs_eta(i,id_anue), h_coeffs_eta(i,id_nux), h_coeffs_eta(i,id_anux),
            h_coeffs_kappa_0_a(i,id_nue), h_coeffs_kappa_0_a(i,id_anue), h_coeffs_kappa_0_a(i,id_nux), h_coeffs_kappa_0_a(i,id_anux),
            h_coeffs_kappa_a(i,id_nue), h_coeffs_kappa_a(i,id_anue), h_coeffs_kappa_a(i,id_nux), h_coeffs_kappa_a(i,id_anux),
            h_coeffs_kappa_s(i,id_nue), h_coeffs_kappa_s(i,id_anue), h_coeffs_kappa_s(i,id_nux), h_coeffs_kappa_s(i,id_anux));
    fprintf(file, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            h_r(i), h_diff_distribution(i),
            h_coeffs_eta_0(i,id_nue), h_coeffs_eta_0(i,id_anue), h_coeffs_eta_0(i,id_nux), h_coeffs_eta_0(i,id_anux),
            h_coeffs_eta(i,id_nue), h_coeffs_eta(i,id_anue), h_coeffs_eta(i,id_nux), h_coeffs_eta(i,id_anux),
            h_coeffs_kappa_0_a(i,id_nue), h_coeffs_kappa_0_a(i,id_anue), h_coeffs_kappa_0_a(i,id_nux), h_coeffs_kappa_0_a(i,id_anux),
            h_coeffs_kappa_a(i,id_nue), h_coeffs_kappa_a(i,id_anue), h_coeffs_kappa_a(i,id_nux), h_coeffs_kappa_a(i,id_anux),
            h_coeffs_kappa_s(i,id_nue), h_coeffs_kappa_s(i,id_anue), h_coeffs_kappa_s(i,id_nux), h_coeffs_kappa_s(i,id_anux));
  }
  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  
  printf("Elapsed time: %.3lf sec\n", cpu_time_used);
  
  fclose(file);

  return;
}

