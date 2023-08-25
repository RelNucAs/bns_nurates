#include <math.h>
#include "functions.h"

// Computation of Fermi-Dirac distribution function

/* Inputs:
 * 	e     [MeV] : energy
 * 	temp  [MeV] : temperature
 * 	mu    [MeV] : chemical potential
*/

double FermiDistr(const double e, const double temp, const double mu) {
  const double arg = (e - mu) / temp;
  double tmp;

  // Handle differently arg>0 and arg<0 cases
  if (arg > 0.) {
    tmp = SafeExp(-arg);
    return tmp / (tmp + 1.);
  } else {
    tmp = SafeExp(arg);
    return 1. / (tmp + 1.);
  }

}

double pair_fermi_dirac(double omega_prime, double eta_prime_electron, double Temp) {
  return 1. / (exp(omega_prime / Temp - eta_prime_electron) + 1.);
}