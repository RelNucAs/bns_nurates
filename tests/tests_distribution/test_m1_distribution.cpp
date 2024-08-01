//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_m1_distribution.cpp
//  \brief Tests the generation of the M1 distribution parameters

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Kokkos_Core.hpp>

#include "distribution.hpp"

using DevExeSpace = Kokkos::DefaultExecutionSpace;
using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;                // increments last index fastest views defined like

void test_distribution() {

  printf("=================================================== \n");
  printf("Testing distribution function reconstruction for M1 \n");
  printf("=================================================== \n");

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_distribution/data/distribution_function_m1_data.txt"; // file containing table generated from Boltzran

  char data_fileout[200] = "/tests/tests_distribution/data/distribution_function_params_from_m1.txt"; // file where data is written
  char data_filepath[300] = {'\0'};

  strcat(data_filepath, filedir);
  strcat(data_filepath, data_fileout);

  strcat(filepath, filedir);
  strcat(filepath, outname);

  printf("Data file path: %s\n", data_filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  const int num_data = 102;
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_T("h_T", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_e("h_mu_e", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_p("h_mu_p", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_n("h_mu_n", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_n_nue("h_n_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_j_nue("h_j_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_chi_nue("h_chi_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_n_anue("h_n_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_j_anue("h_j_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_chi_anue("h_chi_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_n_muon("h_n_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_j_muon("h_j_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_chi_muon("h_chi_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_n_amuon("h_n_amuon", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_j_amuon("h_j_amuon", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_chi_amuon("h_chi_amuon", num_data);

  // read in the data file
  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[0] == '#' && i == 0) {
      continue;
    }

    sscanf(line,
           "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &h_T(i), &h_mu_e(i), &h_mu_p(i), &h_mu_n(i),
           &h_n_nue(i), &h_j_nue(i), &h_chi_nue(i),
           &h_n_anue(i), &h_j_anue(i), &h_chi_anue(i),
           &h_n_muon(i), &h_j_muon(i), &h_chi_muon(i),
           &h_n_amuon(i), &h_j_amuon(i), &h_chi_amuon(i));
    printf("%lf %lf %lf %lf\n", h_T(i), h_mu_p(i), h_mu_e(i), h_mu_n(i));
    i++;
  }
  printf("\n");

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> T("T", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_e("mu_e", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_p("mu_p", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_n("mu_n", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> n_nue("n_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> j_nue("j_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> chi_nue("chi_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> n_anue("n_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> j_anue("j_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> chi_anue("chi_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> n_muon("n_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> j_muon("j_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> chi_muon("chi_muon", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> n_amuon("n_amuon", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> j_amuon("j_amuon", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> chi_amuon("chi_amuon", num_data);

  Kokkos::deep_copy(T, h_T);
  Kokkos::deep_copy(mu_e, h_mu_e);
  Kokkos::deep_copy(mu_p, h_mu_p);
  Kokkos::deep_copy(mu_n, h_mu_n);
  Kokkos::deep_copy(n_nue, h_n_nue);
  Kokkos::deep_copy(j_nue, h_j_nue);
  Kokkos::deep_copy(chi_nue, h_chi_nue);
  Kokkos::deep_copy(n_anue, h_n_anue);
  Kokkos::deep_copy(j_anue, h_j_anue);
  Kokkos::deep_copy(chi_anue, h_chi_anue);
  Kokkos::deep_copy(n_muon, h_n_muon);
  Kokkos::deep_copy(j_muon, h_j_muon);
  Kokkos::deep_copy(chi_muon, h_chi_muon);
  Kokkos::deep_copy(n_amuon, h_n_amuon);
  Kokkos::deep_copy(j_amuon, h_j_amuon);
  Kokkos::deep_copy(chi_amuon, h_chi_amuon);

  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_t_id_nue("h_my_nudistributionparams_w_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_t_id_nue("h_my_nudistributionparams_temp_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_eta_t_id_nue("h_my_nudistributionparams_eta_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_f_id_nue("h_my_nudistributionparams_w_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_f_id_nue("h_my_nudistributionparams_temp_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_c_f_id_nue("h_my_nudistributionparams_c_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_beta_f_id_nue("h_my_nudistributionparams_beta_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_t_id_anue("h_my_nudistributionparams_w_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_t_id_anue("h_my_nudistributionparams_temp_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_eta_t_id_anue("h_my_nudistributionparams_eta_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_f_id_anue("h_my_nudistributionparams_w_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_f_id_anue("h_my_nudistributionparams_temp_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_c_f_id_anue("h_my_nudistributionparams_c_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_beta_f_id_anue("h_my_nudistributionparams_beta_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_t_id_nux("h_my_nudistributionparams_w_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_t_id_nux("h_my_nudistributionparams_temp_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_eta_t_id_nux("h_my_nudistributionparams_eta_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_w_f_id_nux("h_my_nudistributionparams_w_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_temp_f_id_nux("h_my_nudistributionparams_temp_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_c_f_id_nux("h_my_nudistributionparams_c_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_my_nudistributionparams_beta_f_id_nux("h_my_nudistributionparams_beta_f_id_nux", num_data);

  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_t_id_nue("my_nudistributionparams_w_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_t_id_nue("my_nudistributionparams_temp_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_eta_t_id_nue("my_nudistributionparams_eta_t_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_f_id_nue("my_nudistributionparams_w_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_f_id_nue("my_nudistributionparams_temp_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_c_f_id_nue("my_nudistributionparams_c_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_beta_f_id_nue("my_nudistributionparams_beta_f_id_nue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_t_id_anue("my_nudistributionparams_w_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_t_id_anue("my_nudistributionparams_temp_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_eta_t_id_anue("my_nudistributionparams_eta_t_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_f_id_anue("my_nudistributionparams_w_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_f_id_anue("my_nudistributionparams_temp_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_c_f_id_anue("my_nudistributionparams_c_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_beta_f_id_anue("my_nudistributionparams_beta_f_id_anue", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_t_id_nux("my_nudistributionparams_w_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_t_id_nux("my_nudistributionparams_temp_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_eta_t_id_nux("my_nudistributionparams_eta_t_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_w_f_id_nux("my_nudistributionparams_w_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_temp_f_id_nux("my_nudistributionparams_temp_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_c_f_id_nux("my_nudistributionparams_c_f_id_nux", num_data);
  Kokkos::View<double*, LayoutWrapper, DevMemSpace> my_nudistributionparams_beta_f_id_nux("my_nudistributionparams_beta_f_id_nux", num_data);

  char fileline[1000];
  fptr = fopen(data_filepath, "w");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }
  fprintf(fptr,
      "#w_t_nue temp_t_nue eta_t_nue w_f_nue temp_f_nue c_f_nue beta_f_nue w_t_anue temp_t_anue eta_t_anue w_f_anue temp_f_anue c_f_anue beta_f_anue w_t_nux temp_t_nux eta_t_nux w_f_nux temp_f_nux c_f_nux beta_f_nux\n");


  printf("got here! just before start of loop\n");
  Kokkos::parallel_for("par_for_test_m1_distribution", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data),
  KOKKOS_LAMBDA(const int &i) {
    printf("got here! just beginning of loop\n");
    GreyOpacityParams my_grey_opacity_params;
    my_grey_opacity_params.opacity_flags = {.use_abs_em = 0,
                                            .use_pair = 0,
                                            .use_brem = 0,
                                            .use_inelastic_scatt = 0,
                                            .use_iso = 0};

    my_grey_opacity_params.eos_pars.mu_e = mu_e(i);
    my_grey_opacity_params.eos_pars.mu_p = mu_p(i);
    my_grey_opacity_params.eos_pars.mu_n = mu_n(i);
    my_grey_opacity_params.eos_pars.temp = T(i);
    my_grey_opacity_params.m1_pars.n[id_nue] = n_nue(i);
    my_grey_opacity_params.m1_pars.J[id_nue] = j_nue(i);
    my_grey_opacity_params.m1_pars.chi[id_nue] = chi_nue(i);
    my_grey_opacity_params.m1_pars.n[id_anue] = n_anue(i);
    my_grey_opacity_params.m1_pars.J[id_anue] = j_anue(i);
    my_grey_opacity_params.m1_pars.chi[id_anue] = chi_anue(i);
    my_grey_opacity_params.m1_pars.n[id_nux] = n_muon(i);
    my_grey_opacity_params.m1_pars.J[id_nux] = j_muon(i);
    my_grey_opacity_params.m1_pars.chi[id_nux] = chi_muon(i);
    my_grey_opacity_params.m1_pars.n[id_anux] = n_muon(i);
    my_grey_opacity_params.m1_pars.J[id_anux] = j_muon(i);
    my_grey_opacity_params.m1_pars.chi[id_anux] = chi_muon(i);

    printf("got here! before calculate distribution params\n");
    NuDistributionParams my_nudistributionparams = CalculateDistrParamsFromM1(&my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

    my_nudistributionparams_w_t_id_nue(i) = my_nudistributionparams.w_t[id_nue];
    my_nudistributionparams_temp_t_id_nue(i) = my_nudistributionparams.temp_t[id_nue];
    my_nudistributionparams_eta_t_id_nue(i) = my_nudistributionparams.eta_t[id_nue];
    my_nudistributionparams_w_f_id_nue(i) = my_nudistributionparams.w_f[id_nue];
    my_nudistributionparams_temp_f_id_nue(i) = my_nudistributionparams.temp_f[id_nue];
    my_nudistributionparams_c_f_id_nue(i) = my_nudistributionparams.c_f[id_nue];
    my_nudistributionparams_beta_f_id_nue(i) = my_nudistributionparams.beta_f[id_nue];

    my_nudistributionparams_w_t_id_anue(i) = my_nudistributionparams.w_t[id_anue];
    my_nudistributionparams_temp_t_id_anue(i) = my_nudistributionparams.temp_t[id_anue];
    my_nudistributionparams_eta_t_id_anue(i) = my_nudistributionparams.eta_t[id_anue];
    my_nudistributionparams_w_f_id_anue(i) = my_nudistributionparams.w_f[id_anue];
    my_nudistributionparams_temp_f_id_anue(i) = my_nudistributionparams.temp_f[id_anue];
    my_nudistributionparams_c_f_id_anue(i) = my_nudistributionparams.c_f[id_anue];
    my_nudistributionparams_beta_f_id_anue(i) = my_nudistributionparams.beta_f[id_anue];

    my_nudistributionparams_w_t_id_nux(i) = my_nudistributionparams.w_t[id_nux];
    my_nudistributionparams_temp_t_id_nux(i) = my_nudistributionparams.temp_t[id_nux];
    my_nudistributionparams_eta_t_id_nux(i) = my_nudistributionparams.eta_t[id_nux];
    my_nudistributionparams_w_f_id_nux(i) = my_nudistributionparams.w_f[id_nux];
    my_nudistributionparams_temp_f_id_nux(i) = my_nudistributionparams.temp_f[id_nux];
    my_nudistributionparams_c_f_id_nux(i) = my_nudistributionparams.c_f[id_nux];
    my_nudistributionparams_beta_f_id_nux(i) = my_nudistributionparams.beta_f[id_nux];

    printf("got here! 5\n");
  });

  Kokkos::deep_copy(h_my_nudistributionparams_w_t_id_nue, my_nudistributionparams_w_t_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_t_id_nue, my_nudistributionparams_temp_t_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_eta_t_id_nue, my_nudistributionparams_eta_t_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_w_f_id_nue, my_nudistributionparams_w_f_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_f_id_nue, my_nudistributionparams_temp_f_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_c_f_id_nue, my_nudistributionparams_c_f_id_nue);
  Kokkos::deep_copy(h_my_nudistributionparams_beta_f_id_nue, my_nudistributionparams_beta_f_id_nue);

  Kokkos::deep_copy(h_my_nudistributionparams_w_t_id_anue, my_nudistributionparams_w_t_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_t_id_anue, my_nudistributionparams_temp_t_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_eta_t_id_anue, my_nudistributionparams_eta_t_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_w_f_id_anue, my_nudistributionparams_w_f_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_f_id_anue, my_nudistributionparams_temp_f_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_c_f_id_anue, my_nudistributionparams_c_f_id_anue);
  Kokkos::deep_copy(h_my_nudistributionparams_beta_f_id_anue, my_nudistributionparams_beta_f_id_anue);

  Kokkos::deep_copy(h_my_nudistributionparams_w_t_id_nux, my_nudistributionparams_w_t_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_t_id_nux, my_nudistributionparams_temp_t_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_eta_t_id_nux, my_nudistributionparams_eta_t_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_w_f_id_nux, my_nudistributionparams_w_f_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_temp_f_id_nux, my_nudistributionparams_temp_f_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_c_f_id_nux, my_nudistributionparams_c_f_id_nux);
  Kokkos::deep_copy(h_my_nudistributionparams_beta_f_id_nux, my_nudistributionparams_beta_f_id_nux);

  for (int i = 0; i < num_data; i++) {
    fprintf(fptr, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            h_my_nudistributionparams_w_t_id_nue(i), h_my_nudistributionparams_temp_t_id_nue(i),
            h_my_nudistributionparams_eta_t_id_nue(i), h_my_nudistributionparams_w_f_id_nue(i),
            h_my_nudistributionparams_temp_f_id_nue(i), h_my_nudistributionparams_c_f_id_nue(i),
            h_my_nudistributionparams_beta_f_id_nue(i),
            h_my_nudistributionparams_w_t_id_anue(i), h_my_nudistributionparams_temp_t_id_anue(i),
            h_my_nudistributionparams_eta_t_id_anue(i), h_my_nudistributionparams_w_f_id_anue(i),
            h_my_nudistributionparams_temp_f_id_anue(i), h_my_nudistributionparams_c_f_id_anue(i),
            h_my_nudistributionparams_beta_f_id_anue(i),
            h_my_nudistributionparams_w_t_id_nux(i), h_my_nudistributionparams_temp_t_id_nux(i),
            h_my_nudistributionparams_eta_t_id_nux(i), h_my_nudistributionparams_w_f_id_nux(i),
            h_my_nudistributionparams_temp_f_id_nux(i), h_my_nudistributionparams_c_f_id_nux(i),
            h_my_nudistributionparams_beta_f_id_nux(i));
  }

  fclose(fptr);
}
int main() {

  Kokkos::initialize();

  test_distribution();

  Kokkos::finalize();

  return 0;
}
