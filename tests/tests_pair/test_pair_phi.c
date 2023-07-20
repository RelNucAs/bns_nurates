//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_phi.c
//  \brief Compare PairPhi with old python results

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../src/opacities/kernels/kernels.h"

int main() {
  double pair_phi_e[3];

  pair_phi_e[0] = PairPhi(0, 0.5, 0.9, 8, 13, 0);
  pair_phi_e[1] = PairPhi(1, 1.5, 0.9, 8, 13, 0);
  pair_phi_e[2] = PairPhi(2, 3.5, 1.9, 12, 13, 0);

  double max_error = -42.;

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
    if (fabs(pair_phi_e[i] - pair_phi_e_python[i]) > max_error) {
      max_error = fabs(pair_phi_e[i] - pair_phi_e_python[i]);
    }
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
    if (fabs(pair_phi_x[i] - pair_phi_x_python[i]) > max_error) {
      max_error = fabs(pair_phi_x[i] - pair_phi_x_python[i]);
    }
  }

  if (max_error < 1e-14) {
    return 0;
  } else {
    return 1;
  }

}
