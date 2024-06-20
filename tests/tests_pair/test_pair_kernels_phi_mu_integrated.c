//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_kernels_phi_mu_integrated.c
//  \brief Test pair kernels which have been integrated over phi' and mu'

#include <stdio.h>
#include <math.h>
#include "../../include/bns_nurates.h"
#include "../../include/constants.h"
#include "../../include/integration.h"
#include "opacities.h"
#include "../../include/kernels.h"

#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))

#define M_PI_VAL 3.14159265358979323846264338328 /* pi */

int main() {
  double max_percentage_error = -42.;

  double temp = 10.;
  double mu_e = 84.181299999999993;
  double pair_kernels_python[18][7] = { // nutype eta omega(MeV)  omega_prime(MeV) cos_theta R_p(cm^3s^-1) R_a(cm^3s^-1)
      {0, 8.4181299999999997, 2, 2, 0, 3.7344281862494617e-37, 5.5711121998146404e-37},
      {0, 8.4181299999999997, 2, 5, 0, 7.5952594997506983e-37, 1.5294974381563826e-36},
      {0, 8.4181299999999997, 2, 50, 0, 7.6816308521376666e-37, 1.3924664458243231e-34},
      {0, 8.4181299999999997, 5, 2, 0, 8.5842096421064313e-37, 1.7286475408285996e-36},
      {0, 8.4181299999999997, 5, 5, 0, 1.7426501600901373e-36, 4.7370142635342659e-36},
      {0, 8.4181299999999997, 5, 50, 0, 1.6703173370166646e-36, 4.0871317668903464e-34},
      {0, 8.4181299999999997, 50, 2, 0, 3.6039543137528758e-36, 6.5329687806960597e-34},
      {0, 8.4181299999999997, 50, 5, 0, 7.222632008133751e-36, 1.7673197821036537e-33},
      {0, 8.4181299999999997, 50, 50, 0, 3.7405176147900196e-36, 8.2390383297044382e-32},
      {1, 8.4181299999999997, 2, 2, 0, 7.9822789333374034e-38, 1.1908160856214353e-37},
      {1, 8.4181299999999997, 2, 5, 0, 1.7098488384966924e-37, 3.443212727887964e-37},
      {1, 8.4181299999999997, 2, 50, 0, 4.1187070226847197e-37, 7.4660725562898844e-35},
      {1, 8.4181299999999997, 5, 2, 0, 1.7484861258538989e-37, 3.5210186699128533e-37},
      {1, 8.4181299999999997, 5, 5, 0, 3.7248861049956034e-37, 1.0125290212289139e-36},
      {1, 8.4181299999999997, 5, 50, 0, 8.4196498590303211e-37, 2.0602203929942999e-34},
      {1, 8.4181299999999997, 50, 2, 0, 5.2266221327540312e-37, 9.4744151143860744e-35},
      {1, 8.4181299999999997, 50, 5, 0, 1.0588883354811727e-36, 2.5910143286093217e-34},
      {1, 8.4181299999999997, 50, 50, 0, 7.995294986861833e-37, 1.7610809154750178e-32}
  };
  double prefactor_correction = (kClight * POW3(kH * kClight));

  printf("=======================================================================================\n");
  printf("Pair process: Comparison of pair kernels, integrated pair kernels and old Python result\n");
  printf("=======================================================================================\n");

  printf(
      "R_em_e R_em_e (python) R_em_e (integrated) R_abs_e R_abs_e (python) R_abs_e (integrated) R_em_x R_em_x (python) R_em_x (integrated) R_abs_x R_abs_x (python) R_abs_x (integrated)\n");

  for (int i = 0; i < 9; i++) {
    MyEOSParams eos_params = {.temp = temp, .mu_e= mu_e};
    PairKernelParams pair_kernel_params =
        {.omega_prime = pair_kernels_python[i][3], .omega = pair_kernels_python[i][2], .mu_prime = -42., .mu = -42., .cos_theta = pair_kernels_python[i][4], .lmax = 0, .filter = 0};
    MyKernelQuantity result = PairKernels(&eos_params, &pair_kernel_params);
    MyKernelQuantity result_phi_mu_integrated = PairKernelsPhiMuIntegrated(&eos_params, &pair_kernel_params);

    result_phi_mu_integrated.em_e = result_phi_mu_integrated.em_e * prefactor_correction;
    result_phi_mu_integrated.abs_e = result_phi_mu_integrated.abs_e * prefactor_correction;
    result_phi_mu_integrated.em_x = result_phi_mu_integrated.em_x * prefactor_correction;
    result_phi_mu_integrated.abs_x = result_phi_mu_integrated.abs_x * prefactor_correction;

    printf("%0.6e %0.6e %0.6e | %0.6e %0.6e %0.6e | %0.6e %0.6e %0.6e | %0.6e %0.6e %0.6e\n",
           result.em_e, pair_kernels_python[i][5], result_phi_mu_integrated.em_e / (4. * M_PI_VAL),
           result.abs_e, pair_kernels_python[i][6], result_phi_mu_integrated.abs_e / (4. * M_PI_VAL),
           result.em_x, pair_kernels_python[i + 9][5], result_phi_mu_integrated.em_x / (4. * M_PI_VAL),
           result.abs_x, pair_kernels_python[i + 9][6], result_phi_mu_integrated.abs_x / (4. * M_PI_VAL));

    if (fabs((result_phi_mu_integrated.em_e / (4. * M_PI_VAL) - result.em_e) / result.em_e) * 100. > max_percentage_error) {
      max_percentage_error = fabs((result_phi_mu_integrated.em_e / (4. * M_PI_VAL) - result.em_e) / result.em_e) * 100.;
    }
    if (fabs((result_phi_mu_integrated.abs_e / (4. * M_PI_VAL) - result.abs_e) / result.abs_e) * 100 > max_percentage_error) {
      max_percentage_error = fabs((result_phi_mu_integrated.abs_e / (4. * M_PI_VAL) - result.abs_e) / result.abs_e) * 100;
    }
    if (fabs((result_phi_mu_integrated.em_x / (4. * M_PI_VAL) - result.em_x) / result.em_x) * 100. > max_percentage_error) {
      max_percentage_error = fabs((result_phi_mu_integrated.em_x / (4. * M_PI_VAL) - result.em_x) / result.em_x) * 100.;
    }
    if (fabs((result_phi_mu_integrated.abs_x / (4. * M_PI_VAL) - result.abs_x) / result.abs_x) * 100. > max_percentage_error) {
      max_percentage_error = fabs((result_phi_mu_integrated.abs_x / (4. * M_PI_VAL) - result.abs_x) / result.abs_x) * 100.;
    }
  }

  printf("Maximum percentage error: %.16e\n", max_percentage_error);

  if (max_percentage_error < 0.01) {
    return 0;
  } else {
    return 1;
  }
}
