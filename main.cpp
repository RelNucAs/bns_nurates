// This contains all the main functions which implement the various opactities due to neutral current interactions (elastic scattering)
// Based on Bruenn(1985) and Burrows et al(2004). List of interactions considered: 
// 1: electron neutrino absorption on neutrons
// 2: electron antineutrino absorption on protons

//Libraries to be loaded

#include <iostream> //Input&Output stream model based library
#include <string>
#include <sstream>
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

int main ()
{
	ifstream inputfile;
    inputfile.open("./Input/nurates_1.008E+01.txt");
	ofstream outputfile("./Output/nu_rates.txt");
		int zone;
		double d1;
	  	double r, rho, temp, ye, mu_e, mu_np, yh, ya, yp, yn, em_nue, l_nue, em_anue, l_anue;
        string header1;
        string header2;

        if (!inputfile)
        {
        	cout << "File not found" << '\n';
        }

        //two header lines
        getline(inputfile,header1); getline(inputfile,header2);
        outputfile << header1 << '\n'; 
        outputfile << "#" << " " << "r[cm]" << " " << "rho[g/cm^3]" << " " << "T[MeV]" << " " \
        << "Ye" << " " << "mu_e[MeV]" << " " << "mu_hat[MeV]" << " " << "Yp" << " " << \
        "Yn" << " " << "1/l_nue[1/cm]" << " " << "1/l_anue[1/cm]" << '\n';

	  	while (true)
		{
			inputfile >> zone >> r >> rho >> temp >> ye >> mu_e >> mu_np \
			>> yh >> ya >> yp >> yn >> em_nue >> l_nue >> em_anue >> l_anue;
 			double nn = rho*yn/mb;
 			double np = rho*yp/mb;
			double e_nu = 10.079454858851342;
			double e_nu_bar = 10.079454858851342;
    		double nun;
    		double nup;

			if(inputfile.eof()) break;
			nun = nu_n_abs(nn, np, temp, ye, e_nu, mu_e, mu_np);
			nup = nu_p_abs(nn, np, temp, ye, e_nu_bar, mu_e, mu_np);
			outputfile << r  << " " << rho << " " << temp << " " << \
			ye << " " << mu_e << " " << mu_np << " " << yp << " " << yn << " " << nun << " " << nup <<'\n';
			cout << zone << " " << r  << " " << rho << " " << temp << " " << \
			ye << " " << mu_e << " " << mu_np << " " << yp << " " << yn << '\n';
			//printf ("nu_n_abs %.3e cm^-1\n", nun);
			//printf ("anu_p_abs %.3e cm^-1\n", nup);
     	}
    inputfile.close();
    outputfile.close();

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
