//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_f.c
//  \brief Compare the pairF functions with the old Python code

#include <stdio.h>
#include <math.h>
#include "../../include/kernels.hpp"

int main() {
  double pair_f[6];

  pair_f[0] = PairFOptimized(0, 2.2, 5.0);
  pair_f[1] = PairFOptimized(1, 2.2, 5.0);
  pair_f[2] = PairFOptimized(0, -2.2, 5.0);
  pair_f[3] = PairFOptimized(1, -2.2, 5.0);
  pair_f[4] = PairFOptimized(0, 6.2, 5.0);
  pair_f[5] = PairFOptimized(5, 6.2, 5.0);

  double pair_f_python[6] = {
      2.246050458762737,
      3.601983945397414,
      0.10433698594797482,
      0.1033982406148323,
      4.738745038671483,
      2233.422169785745
  };

  double max_error = -42.;

  printf("==================================================\n");
  printf("Pair procrss: Comparison of F with old Python code\n");
  printf("==================================================\n");

  printf("F(ours) \t\t F(python) \t |diff|\n");
  for (int i = 0; i < 6; i++) {
    printf("%0.6e \t %0.6e \t %0.6e\n", pair_f[i], pair_f_python[i], fabs(pair_f[i] - pair_f_python[i]));
    if (fabs(pair_f[i] - pair_f_python[i]) > max_error) {
      max_error = fabs(pair_f[i] - pair_f_python[i]);
    }
  }

  printf("Maximum error: %.16e\n", max_error);

  if (max_error < 1e-10) {
    return 0;
  } else {
    return 1;
  }

}