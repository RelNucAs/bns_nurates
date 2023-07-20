//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_f.c
//  \brief Compare the pairF functions with the old Python code

#include <stdio.h>
#include <math.h>
#include "../../src/opacities/kernels/kernels.h"

int main() {
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

  double max_error = -42.;

  printf("PairF comparison\n");
  printf("================\n");
  printf("F(ours) \t\t F(python) \t\t\t |diff|\n");
  for (int i = 0; i < 6; i++) {
    printf("%0.16e \t %0.16e \t %0.16e\n", pair_f[i], pair_f_python[i], fabs(pair_f[i] - pair_f_python[i]));
    if (fabs(pair_f[i] - pair_f_python[i]) > max_error) {
      max_error = fabs(pair_f[i] - pair_f_python[i]);
    }
  }

  if (max_error < 1e-10) {
    return 0;
  } else {
    return 1;
  }

}