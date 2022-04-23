// This contains all the main functions which implement the various opactities due to neutral current interactions (elastic scattering)
// Based on Bruenn(1985) and Burrows et al(2004). List of interactions considered: 
// 1: electron neutrino absorption on neutrons
// 2: electron antineutrino absorption on protons

//Libraries to be loaded

#include <iostream> //Input&Output stream model based library
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates
#include "weak_rates.h" //Header file containing all rates

using namespace constants;
using namespace corrections;
using namespace weakrates;

int main (){
	double nb = 10e2; //fm^3
	double temp = 1.2040;
	double ye = 0.134;
	double yn = 6.8660*(10e-1);
	double yp = 3.1340*(10e-1);
	double nn = nb*yn;
	double np = nb*yp;
	double e_nu = 10.079454858851342;
	double e_nu_bar = 10.079454858851342;
	double mu_e = 249.51;
	double mu_np = 6.015;
    double nun;
    double nup;

	nun = nu_n_abs(nn, np, temp, ye, e_nu, mu_e, mu_np);
	nup = nu_p_abs(nn, np, temp, ye, e_nu_bar, mu_e, mu_np);
	printf ("nu_n_abs %f\n", nun);
	printf ("nu_p_abs %f\n", nup);
}
