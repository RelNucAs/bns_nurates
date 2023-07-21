//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_t.c
//  \brief Compare the PairT functions with the old python code

#include <stdio.h>
#include <math.h>
#include "../../src/opacities/kernels/kernels.h"

int main() {
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

  double max_error = -42.;

  printf("==================================================\n");
  printf("Pair process: Comparison of T with old Python code\n");
  printf("==================================================\n");

  printf("T(ours) \t\t T(python) \t |diff|\n");
  for (int i = 0; i < 5; i++) {
    printf("%0.6e \t %0.6e \t %0.6e\n", pair_t[i], pair_t_python[i], fabs(pair_t[i] - pair_t_python[i]));
    if (fabs(pair_t[i] - pair_t_python[i]) > max_error) {
      max_error = fabs(pair_t[i] - pair_t_python[i]);
    }
  }

  printf("Maximum error: %0.16e\n", max_error);

  if (max_error < 1e-10) {
    return 0;
  } else {
    return 1;
  }

}