#include "functions.h"

// Computation of theta step function 

double StepFunction(const double x) {
  if (x < 0.) {
    return 0.;
  } else {
    return 1.;
  }
}