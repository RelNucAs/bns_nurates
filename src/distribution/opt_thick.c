#include <math.h>

#include "distribution.h"
#include "../bns_nurates.h"
#include "../functions/functions.h"

// Definition of constants
const double kFourThirds = 4./3.;
const double kMinERatio = 3.1;

// Neutrino distribution function in optically thick regime
double NuFThick(double omega, NuDistributionParams *distr_pars) {
  double temp = distr_pars->temp_t;
  double mu = temp * distr_pars->eta_t;
  return FermiDistr(omega, temp, mu);
}

// Recover parameters of thick distribution function from M1 quantities
// Input  -> n, J from M1Quantities, T from MyEOSParams
// Output -> temp_t, eta_t in NuDistributionParams
void ThickFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars,
                 NuDistributionParams *out_distr_pars) {
  static double kExpThree = exp(3.);

  double arg = M1_pars->J / (M1_pars->n * eos_pars->temp); // average neutrino energy over thermal energy

  // @TODO: change the following to work only in debug mode
  if (arg <= kMinERatio) {
    printf("WARNING -> opt_thick.c: Neutrino average enerrgy - thermal energy ratio is too small");
  }
  arg = (arg > kMinERatio) ? arg : kMinERatio;
  
  out_distr_pars->w_t = 0.5 * (3. * M1_pars->chi - 1.);
  out_distr_pars->temp_t = eos_pars->temp;
  out_distr_pars->eta_t  = kFourThirds * log( SafeExp(arg)- kExpThree );
                            
  return;
}