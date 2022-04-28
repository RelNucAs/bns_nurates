// This contains all the main functions which implement the various opactities due to neutral current interactions (elastic scattering)
// Based on Bruenn(1985) and Burrows et al(2004). List of interactions considered: 
// 1: electron neutrino absorption on neutrons
// 2: electron antineutrino absorption on protons

//Libraries to be loaded

#include <iostream> //Input&Output stream model based library
#include <fstream> // both read and write from/to files
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates
#include "weak_rates.h" //Header file containing all rates

using namespace constants;
using namespace corrections;
using namespace weakrates;
using namespace std;

int main (){
	ifstream input_file("./Input/nurates_1.008E+01.txt");
	ofstream output_file("./Output/nu_rates.txt")
	if((input_file.is_open()) && (output_file.is_open()))
	{
		int zone;
	  	double r, rho, temp, ye, mu_e, mu_np, yh, ya, yp, yn;

	  	while (!input_file.eof())
		{
			input_file >> zone >> r >> rho >> temp >> ye >> mu_e >> mu_np >> yh >> ya >> yp >> yn
			double mb = 1.674e-24; // g
 			double nn = rho*yn/mb;
 			double np = rho*yp/mb;
			double e_nu = 10.079454858851342;
			double e_nu_bar = 10.079454858851342;
			double mu_e = 249.51; // MeV
			double mu_np = 60.156; // MeV
    		double nun;
    		double nup;

			nun = nu_n_abs(nn, np, temp, ye, e_nu, mu_e, mu_np);
			nup = nu_p_abs(nn, np, temp, ye, e_nu_bar, mu_e, mu_np);
			output_file << r << rho << temp << ye << mu_e << mu_np << yp << yn << nun << nup	
     	}
    input_file.close();
    output_file.close();
	}

	//double rho = 3.73e14; // g/cm3
 	//double temp = 12.040; // MeV
 	//double ye = 0.31340;
 	//double yn = 6.8660e-1;
 	//double yp = 3.1340e-1;
 	//double mb = 1.674e-24; // g
 	//double nn = rho*yn/mb;
 	//double np = rho*yp/mb;
	//double e_nu = 10.079454858851342;
	//double e_nu_bar = 10.079454858851342;
	//double mu_e = 249.51; // MeV
	//double mu_np = 60.156; // MeV
    //double nun;
    //double nup;

	//nun = nu_n_abs(nn, np, temp, ye, e_nu, mu_e, mu_np);
	//nup = nu_p_abs(nn, np, temp, ye, e_nu_bar, mu_e, mu_np);
	//printf ("nu_n_abs %.3e cm^-1\n", nun);
	//printf ("anu_p_abs %.3e cm^-1\n", nup);
}
