#include <math.h>
#include <float.h>

#include "functions.h"

// Implementations of the Heaviside Step function
//
// Continuous approximation is taken from Bierng-Chearl Ahn (2013)

#define a_heaviside 20.

// Implementation with if statement
double HeavisidePiecewise(const double x) {
  if (x < 0.) {
    return 0.;
  } else {
    return 1.;
  }
}

// Continuous approximation with tanh - (Eq.5)
double HeavisideTanhApprox(const double x) {
  return 0.5 * (1. + tanh(2. * a_heaviside * x / (fabs(x) + DBL_MIN))); //DBL_TRUE_MIN
}