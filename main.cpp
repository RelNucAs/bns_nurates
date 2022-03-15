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
	double nb = 3.7300E+12;
	double temp = 1.2040;
	double ye = 0.134;
  double ee = 10;
  double pe = 10;
  double ne = 10;
	double e_nu = 10.079454858851342;
	double e_nu_bar = 10.079454858851342;
	double mu_e = 249.51;
	double mu_np = 6.015;
  double nun;
  double nup;

	nun = nu_n_abs(nb, temp, ye, e_nu, e_nu_bar, mu_e, mu_np, ne, pe, ee);
	nup = nu_p_abs(nb, temp, ye, e_nu, e_nu_bar, mu_e, mu_np, ne, pe, ee);
	std::cout << "Neutrino absorption on neutrons for" << "nb =" << nb << "," << "temp =" << temp << "is" << nup <<std::endl;
	std::cout << "Neutrino absorption on protons for" << "nb =" << nb << "," << "temp =" << temp << "is" << nun <<std::endl;
}
