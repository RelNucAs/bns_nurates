#ifndef FERMI_H
#define FERMI_H

// \file  fermi.hpp
// \brief Compute Fermi-Dirac distribution function

/* Inputs:
 * 	E  [MeV] : energy
 * 	T  [MeV] : temperature
 * 	mu [MeV] : chemical potential
*/

/* Fermi-Dirac distribution function */
double Fermi(const double mu, const double E, const double T) {
  double tmp;
  const double arg = (E-mu) / T;

  if (arg > 0.) {
	  tmp = exp(-arg);
	  return tmp / (tmp+1.);
  } else {
	  tmp = exp(arg);
	  return 1.  / (tmp+1.);
  }
}

#endif
