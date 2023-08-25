#include <math.h>
#include <stdio.h>

#include "functions.h"

/*
 * Computation of Lambert function using logarithmic recursive formula
 * in "Loczi, Applied Mathematics and Computation, 433 (2022)"
 */

// @TODO: the following constant has already been defined in src/distribution/opt_thin.c
#define kExp exp(1.)

// Initial guess for recursion
double BetaZero(double x) {
  if (x > kExp) {
    double logx = log(x);
    return logx - log(logx);
  } else if ((x > 0.) && (x <= kExp)) {
    return x / kExp;
  } else if ((x > 1./kExp) && (x <= 0.)) {
    double tmp_1 = kExp * x;
    double tmp_2 = sqrt(1 + tmp_1);
    double tmp_3 = 1. + tmp_2;
    return tmp_1 * log(tmp_3) / (tmp_2 * tmp_3);
  }

  return -1.;
}

// Recursive formula for W0 real branch of Lambert function
double W0(double x) {
  static const int kMaxIterW0 = 5;
  int n = 0;

  // @TODO: debug-only + kExp defined above
  if (x < 1./exp(1.)) {   
    printf("WARNING -> lambert.c: Lambert function is not defined for x < 1/e, returning -1");
    return -1;
  }
  
  double beta = BetaZero(x);

  while (n < kMaxIterW0) {
    beta = beta * (1. + log(x / beta)) / (1. + beta);
    n++;
  }

  return beta;
}
