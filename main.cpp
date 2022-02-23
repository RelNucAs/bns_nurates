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

using namespace constants;
using namespace corrections;


int main (){
	double nb = ;
	double temp = ;
	double ye = ;
	double e_nu = ;
	double e_nu_bar = ;
	double mu_e = ;
	double mu_np = ;
	nun = nu_n_abs(nb, temp, ye, e_nu, mu_e, mu_np);
	nup = nu_p_abs(nb, temp, ye, e_nu_bar, mu_e, mu_np);
	cout << "Neutrino absorption on neutrons for" << "rho =" << rho << "," << "temp =" << temp << "is" << nup <<endl;
	cout << "Neutrino absorption on protonss for" << "rho =" << rho << "," << "temp =" << temp << "is" << nun <<endl;
}
