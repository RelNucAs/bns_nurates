//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_psi.c
//  \brief Compare the PairG functions with the old Python results

#include <stdio.h>
#include <math.h>
#include "../../src/opacities/kernels/kernels.h"

int main() {
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

  double max_error = -42.;

  printf("PairPsi comparison\n");
  printf("==================\n");
  printf("Psi(ours) \t\t Psi(python) \t\t |diff|\n");
  for (int i = 0; i < 8; i++) {
    printf("%0.16e  %0.16e  %0.16e\n", pair_psi[i], pair_psi_python[i], fabs(pair_psi[i] - pair_psi_python[i]));
    if (fabs(pair_psi[i] - pair_psi_python[i]) > max_error) {
      max_error = fabs(pair_psi[i] - pair_psi_python[i]);
    }
  }

  if (max_error < 1e-6) {
    return 0;
  } else {
    return 1;
  }

}