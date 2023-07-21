//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  tests_pair_g.c
//  \brief Compare the PairG functions with the old Python results

#include <stdio.h>
#include <math.h>
#include "../../src/opacities/kernels/kernels.h"

int main() {
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

  double max_error = -42.;

  printf("==================================================\n");
  printf("Pair process: Comparison of G with old Python code\n");
  printf("==================================================\n");
  printf("G(ours) \t\t G(python) \t |diff|\n");
  for (int i = 0; i < 6; i++) {
    printf("%0.6e \t %0.6e \t %0.6e\n", pair_g[i], pair_g_python[i], fabs(pair_g[i] - pair_g_python[i]));
    if (fabs(pair_g[i] - pair_g_python[i]) > max_error) {
      max_error = fabs(pair_g[i] - pair_g_python[i]);
    }
  }

  printf("Maximum error: %.16e\n", max_error);

  if (max_error < 1e-10) {
    return 0;
  } else {
    return 1;
  }

}