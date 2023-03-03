#ifndef ELASTIC_SC_H
#define ELASTIC_SC_H

// \file nu_elastic_scatt.hpp
// \brief Computation of first two Legendre coefficients for neutrino elastic
//        scattering on nucleons
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

#include "constants.hpp" // Header file for physical constants

using namespace constants;

namespace elastic_scatt {
  const double hbar = 0.5*h/pi; //[MeV*s]	
	
  const double c0 = (4.*pi*GF*GF)/h;
  const double c1 = (4./3.*pi*GF*GF)/h;	
}	

/* Computation of degeneracy parameter eta_NN */
double eta_NN(const double nb, const double T, const double yN);

/* Legendre expansion (up to l=1) of scattering kernel */
double nu_N_scatt(const double omega, const double nb, const double T, const double yN, const int reacflag);

/* Sum proton and neutron scattering kernels */
double nu_N_scatt_tot(const double omega, const double nb, const double T, const double yn, const double yp);

#endif
