//Compute scattering kernel for neutrino elastic scattering on nucleons

#pragma once //compile only once

#include "constants.hpp"

using namespace constants;

namespace elastic_scatt
{
	const double hbar = 0.5*h/pi; //[MeV*s]	
	
	const double c0 = (4.*pi*GF*GF)/h;
	const double c1 = (4./3.*pi*GF*GF)/h;
	
}	

// degeneracy parameter eta_NN
double eta_NN(const double nb, const double T, const double yN);

double nu_N_scatt(const double Enu, const double nb, const double T, const double yN, const int reacflag);

double nu_N_scatt_tot(double Enu, double nb, double T, double yn, double yp);
