//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair.c
//  \brief test routines for the pair process

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests.h"
#include "../src/bns_nurates.h"
#include "../src/opacities/kernels/kernels.h"

// Compare the PairT functions with the old python code
void TestPairT() {
  double pair_t[5];
  pair_t[0] = PairT(1, 0., 1e-6);
  pair_t[1] = PairT(1, 1., 1e-6);
  pair_t[2] = PairT(2, 0., 1e-6);
  pair_t[3] = PairT(11, 0., 1e-6);
  pair_t[4] = PairT(11, 17., 1e-6);

  double pair_t_python[5] = {
      0.6931471805599453,
      0.3132615578738979,
      0.8224670334241132,
      0.9995171434980608,
      4.139937718785167e-08
  };

  printf("PairT comparison\n");
  printf("================\n");
  printf("T(ours) \t\t T(python) \t\t\t |diff|\n");
  for (int i = 0; i < 5; i++) {
    printf("%0.16e \t %0.16e \t %0.16e\n", pair_t[i], pair_t_python[i], fabs(pair_t[i] - pair_t_python[i]));
  }

}

// Compare the pairF functions with the old Python code
void TestPairF() {
  double pair_f[6];

  pair_f[0] = PairF(0, 2.2, 5.0);
  pair_f[1] = PairF(1, 2.2, 5.0);
  pair_f[2] = PairF(0, -2.2, 5.0);
  pair_f[3] = PairF(1, -2.2, 5.0);
  pair_f[4] = PairF(0, 6.2, 5.0);
  pair_f[5] = PairF(5, 6.2, 5.0);

  double pair_f_python[6] = {
      2.246050458762737,
      3.601983945397414,
      0.10433698594797482,
      0.1033982406148323,
      4.738745038671483,
      2233.422169785745
  };

  printf("PairF comparison\n");
  printf("================\n");
  printf("F(ours) \t\t F(python) \t\t\t |diff|\n");
  for (int i = 0; i < 6; i++) {
    printf("%0.16e \t %0.16e \t %0.16e\n", pair_f[i], pair_f_python[i], fabs(pair_f[i] - pair_f_python[i]));
  }

}

// Compare the PairG functions with the old Python results
void TestPairG() {
  double pair_g[5];

  pair_g[0] = PairG(0, 3, 7, 8.2, 3.0, 5.67);
  pair_g[1] = PairG(1, 3, 7, 8.2, 3.0, 5.67);
  pair_g[2] = PairG(6, 3, 7, 8.2, 3.0, 5.67);
  pair_g[3] = PairG(6, 3, 7, -8.2, 3.0, 5.67);
  pair_g[4] = PairG(3, 3, 15, -8.2, 3.0, 5.67);

  double pair_g_python[5] = {
      -0.2577301768928941,
      -1.5505542691166285,
      -15358.570715562628,
      -586.9670564136454,
      -6.071262421520971
  };

  printf("PairG comparison\n");
  printf("================\n");
  printf("G(ours) \t\t G(python) \t\t\t |diff|\n");
  for (int i = 0; i < 6; i++) {
    printf("%0.16e \t %0.16e \t %0.16e\n", pair_g[i], pair_g_python[i], fabs(pair_g[i] - pair_g_python[i]));
  }

}

// compare PairPsi with old Python result
void TestPairPsi() {
  double pair_psi[8];

  pair_psi[0] = PairPsi(0, 0.3, 3.46, 44);
  pair_psi[1] = PairPsi(1, 0.3, 3.46, 44);
  pair_psi[2] = PairPsi(1, 8.3, 3.6, 4);
  pair_psi[3] = PairPsi(1, 5.3, 5.3, 7);
  pair_psi[4] = PairPsi(2, 3.5 / 13., 1.9 / 13., 12);
  pair_psi[5] = PairPsi(2, 1.9 / 13., 3.5 / 13., 12);
  pair_psi[6] = PairPsi(3, 3.5 / 13., 1.9 / 13., 12);
  pair_psi[7] = PairPsi(3, 1.9 / 13., 3.5 / 13., 12);

  double pair_psi_python[8] = {
      9.414416002242998e-17,
      -4.9404125853812454e-17,
      12.152335242563368,
      2.3811601968304785,
      -9.756549503572235e-09,
      -9.210731175879528e-09,
      -2.257559103965137e-06,
      -2.367788655295864e-06
  };

  printf("PairPsi comparison\n");
  printf("==================\n");
  printf("Psi(ours) \t\t Psi(python) \t\t |diff|\n");
  for (int i = 0; i < 8; i++) {
    printf("%0.16e  %0.16e  %0.16e\n", pair_psi[i], pair_psi_python[i], fabs(pair_psi[i] - pair_psi_python[i]));
  }

}

// Compare PairPhi with old python result
void TestPairPhi() {
  double pair_phi_e[3];

  pair_phi_e[0] = PairPhi(0, 0.5, 0.9, 8, 13, 0);
  pair_phi_e[1] = PairPhi(1, 1.5, 0.9, 8, 13, 0);
  pair_phi_e[2] = PairPhi(2, 3.5, 1.9, 12, 13, 0);

  double pair_phi_e_python[3] = {
      1.2020437773872976e-37,
      -1.743683568656641e-37,
      2.999896383246949e-39
  };

  printf("PairPhi 'e' comparison\n");
  printf("==================\n");
  printf("Phi(ours) \t\t Phi(python) \t\t |diff|\n");
  for (int i = 0; i < 3; i++) {
    printf("%0.16e  %0.16e  %0.16e\n", pair_phi_e[i], pair_phi_e_python[i], fabs(pair_phi_e[i] - pair_phi_e_python[i]));
  }

  double pair_phi_x[3];

  pair_phi_x[0] = PairPhi(0, 0.5, 0.9, 8, 13, 1);
  pair_phi_x[1] = PairPhi(1, 1.5, 0.9, 8, 13, 1);
  pair_phi_x[2] = PairPhi(2, 3.5, 1.9, 12, 13, 1);

  double pair_phi_x_python[3] = {
      -4.233132673385275e-39,
      7.8204497976434e-39,
      -1.5966596697987114e-40
  };

  printf("PairPhi 'x' comparison\n");
  printf("======================\n");
  printf("Phi(ours) \t\t Phi(python) \t\t |diff|\n");
  for (int i = 0; i < 3; i++) {
    printf("%0.16e  %0.16e  %0.16e\n", pair_phi_x[i], pair_phi_x_python[i], fabs(pair_phi_x[i] - pair_phi_x_python[i]));
  }

}

void TestPairKernels() {
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

  double pair_kernels[18][7];

  printf("Pair Kernel comparison\n");
  printf("======================\n");
  printf(
      "omega omega' cos_theta R_em_e (new) R_em_e (python) %% err  R_abs_e (new)    R_abs_e (python) %% err  R_em_x (new)     R_em_x (python)  %% err  R_abs_x (new)    R_abs_x (python) %% err\n");

  for (int i = 0; i < 9; i++) {
    MyEOSParams eos_params = {.temp = temp, .mu_e= mu_e};
    PairKernelParams pair_kernel_params =
        {.omega_prime = pair_kernels_python[i][3], .omega = pair_kernels_python[i][2], .mu_prime = -42., .mu = -42., .cos_theta = pair_kernels_python[i][4], .lmax = 0, .filter = 0};
    MyOpacityQuantity result = PairKernels(&eos_params, &pair_kernel_params);

    printf("%0.0e %0.0e %0.0e %0.10e %0.10e |%0.2f| %0.10e %0.10e |%0.2f| %0.10e %0.10e |%0.2f| %0.10e %0.10e |%0.2f|\n",
           pair_kernels_python[i][2],
           pair_kernels_python[i][3],
           pair_kernels_python[i][4],
           result.em_e,
           pair_kernels_python[i][5],
           fabs((result.em_e - pair_kernels_python[i][5]) / pair_kernels_python[i][5]) * 100.,
           result.abs_e,
           pair_kernels_python[i][6], fabs((result.abs_e - pair_kernels_python[i][6]) / result.abs_e - pair_kernels_python[i][6]) * 100,
           result.em_x,
           pair_kernels_python[i + 9][5],
           fabs((result.em_x - pair_kernels_python[i + 9][5]) / pair_kernels_python[i + 9][5]) * 100.,
           result.abs_x,
           pair_kernels_python[i + 9][6],
           fabs((result.abs_x - pair_kernels_python[i + 9][6]) / pair_kernels_python[i + 9][6]) * 100.
    );
  }
}