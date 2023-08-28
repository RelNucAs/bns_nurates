#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "functions.h"

/*
 * Computation of Lambert function using logarithmic recursive formula
 * in "Loczi, Applied Mathematics and Computation, 433 (2022)"
 */

// Recursive formula for W0 real branch of Lambert function
double W0(const double x) {
  static const double e  = 2.718281828459045235360287471352;
  static const double ie = 1. / e;
  double beta;

  assert(x > -ie);

  // Initial guess for recursion
  if (x > e) {
    double log_x = log(x);
    beta = log_x - log(log_x);
  } else if (x > 0.) {
    beta = x / e;
  } else { // if(x > -ie)
    double tmp_1 = e * x;
    double tmp_2 = sqrt(1 + tmp_1);
    double tmp_3 = 1. + tmp_2;
    beta = tmp_1 * log(tmp_3) / (tmp_2 * tmp_3);
  }

  beta = beta / (1. + beta) * (1. + log(x / beta));
  beta = beta / (1. + beta) * (1. + log(x / beta));
  beta = beta / (1. + beta) * (1. + log(x / beta));

  return beta;
}

