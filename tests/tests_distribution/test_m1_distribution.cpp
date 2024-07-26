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

#include <Kokkos_Core.hpp>

#include "opacities.hpp"
#include "m1_opacities.hpp"
#include "../../include/integration.hpp"
#include "../../include/distribution.hpp"
#include "../../include/constants.hpp"
#include "../../include/functions.hpp"

int main() {

  Kokkos::initialize();

  using DevExeSpace = Kokkos::DefaultExecutionSpace;
  using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
  using HostMemSpace = Kokkos::HostSpace;
  using LayoutWrapper = Kokkos::LayoutRight;                // increments last index fastest views defined like
  
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

  __attribute__((unused))
  double energy_vals[] = {3.36092278, 4.28273982, 5.45738823, 6.9542133, 8.86158006, 11.2920898, 14.38922757, 18.33583276,
                          23.36489302, 29.77329872, 37.9393698, 48.34518991, 61.60506618, 78.50179483, 100.03287349, 127.46938843,
                          162.43105312, 206.98182789, 263.75176578, 336.09227757};
  int num_data = 102;
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

  Kokkos::parallel_for("par_for_test_m1_distribution", Kokkos::RangePolicy<>(DevExeSpace(), 0, num_data - 1),
  KOKKOS_LAMBDA(const int &i) {
    GreyOpacityParams my_grey_opacity_params;
    my_grey_opacity_params.opacity_flags = opacity_flags_default_none;

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

    NuDistributionParams my_nudistributionparams = CalculateDistrParamsFromM1(&my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

  });



  char fileline[1000];
  fptr = fopen(data_filepath, "w");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }
  fputs(
      "#w_t_nue temp_t_nue eta_t_nue w_f_nue temp_f_nue c_f_nue beta_f_nue w_t_anue temp_t_anue eta_t_anue w_f_anue temp_f_anue c_f_anue beta_f_anue w_t_nux temp_t_nux eta_t_nux w_f_nux temp_f_nux c_f_nux beta_f_nux\n",
      fptr);
  /*
  for (i = 0; i < num_data; i++) {


    sprintf(fileline, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
            my_nudistributionparams.w_t[id_nue], my_nudistributionparams.temp_t[id_nue], my_nudistributionparams.eta_t[id_nue], my_nudistributionparams.w_f[id_nue],
            my_nudistributionparams.temp_f[id_nue], my_nudistributionparams.c_f[id_nue], my_nudistributionparams.beta_f[id_nue],
            my_nudistributionparams.w_t[id_anue], my_nudistributionparams.temp_t[id_anue], my_nudistributionparams.eta_t[id_anue], my_nudistributionparams.w_f[id_anue],
            my_nudistributionparams.temp_f[id_anue], my_nudistributionparams.c_f[id_anue], my_nudistributionparams.beta_f[id_anue],
            my_nudistributionparams.w_t[id_nux], my_nudistributionparams.temp_t[id_nux], my_nudistributionparams.eta_t[id_nux], my_nudistributionparams.w_f[id_nux],
            my_nudistributionparams.temp_f[id_nux], my_nudistributionparams.c_f[id_nux], my_nudistributionparams.beta_f[id_nux]);
    fputs(fileline, fptr);

  }
  */
  fclose(fptr);

  //Kokkos::finalize();

  return 0;
}
