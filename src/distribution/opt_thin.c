#include <math.h>

#include "distribution.h"
#include "../constants.h"
#include "../bns_nurates.h"
#include "../functions/functions.h"

#define kExp exp(1.)
#define kSqrt32PiThree sqrt(32. * kPi * kPi * kPi)

// Neutrino distribution function in optically thin regime
double NuFThin(double omega, NuDistributionParams *distr_pars) {
  double temp_f = distr_pars->temp_f;
  double c_f = distr_pars->c_f;
  return pow(omega, c_f) * exp( -omega / temp_f);
}

// Recover parameters of thin distribution function from M1 quantities
// Input  -> n, J from M1Quantities
// Output -> temp_f, c_f in NuDistributionParams
void ThinFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars) {
  double n = M1_pars->n;
  double j = M1_pars->J;

  double j_over_n = j / n;
  double A = j_over_n / kExp;
  double logA = log(A);
  double B = n / kSqrt32PiThree;

  double lambert = W0(- 2 * logA / (B * B) );

  // @TODO: change the following to work only in debug mode
  if (lambert < 0.) {
    printf("WARNING -> opt_thin.c: negative value for Lambert function");
  }
  
  lambert = (lambert > 0.) ? lambert : 0.;
  double x = - 0.5 * lambert / logA;
  
  out_distr_pars->w_f = 1.5 * (1. - M1_pars->chi);
  out_distr_pars->c_f = x - 3.;
  out_distr_pars->temp_f  = j_over_n / x;
                            
  return;
}