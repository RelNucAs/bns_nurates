// Main file to compute opacities and emissivities of the various reactions
//Based on Bruenn(1985) and Burrows+(2006)
//List of interactions considered: 
// 1: electron neutrino absorption on neutrons
// 2: electron antineutrino absorption on protons
// 3: neutrino scattering on nucleons

//Libraries to be loaded

#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <string> 
#include <vector>
#include <sstream>
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include <cstdlib> //libabry for abs()

#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates
#include "parameters.h" //Header file containing all parameters
#include "weak_rates.h" //Header file containing all rates
#include "nu_elastic_scatt.h" //Header file containing elastic scattering rates

using namespace constants;
using namespace corrections;
using namespace parameters;
using namespace weakrates;
using namespace elastic_scatt;
using namespace std;

int main (){
	int zone;
	double r, d, T;
	double yh, ya, ye, yn, yp;
	double e_nu = 10.079454858851342;
	double e_nu_bar = 10.079454858851342;
	double mu_e, mu_hat;
	double np, nn;
	double em_nue, ab_nue, em_anue, ab_anue, B_IS_nue, B_IS_anue;
	string dummyLine;
	
	// open input profile
	//string filename = "nurates_1.008E+01.txt";
	std::fstream fin("nurates_1.008E+01.txt"); //read file
        if (!fin)
        {
                cout << "File not found" << '\n';
        }
	for (int i = 0; i < 3; i++) { //skip first two lines and read first one with data
		getline(fin, dummyLine);
	}
	
	// open output file
	string outname;
	if ((use_WM_ab == 1) && (use_WM_sc == 1)) {
		outname = "newrates_WM.txt";
 	} else {
		outname = "newrates.txt";
	}
	ofstream fout(outname);
	//myfile.open("check_rates.txt");
	// file header
	string fhead;
	fhead = "# id, d [g/cm^3], T [MeV], Ye, mu_e [MeV], mu_hat [MeV], Yp, Yn, em_nue, ab_nue[1/cm], em_anue[1/cm], ab_anue[1/cm], B_IS_nue[1/cm], B_IS_anue[1/cm]\n";
	fout << fhead;


	while (!fin.fail()){
		// read data from input table
		std::stringstream ss(dummyLine);
		ss >> zone >> r >> d >> T >> ye >> mu_e >> mu_hat >> yh >> ya >> yp >> yn;

		// electron (anti)neutrino absorption 
		std::tie(em_nue,ab_nue)   = nu_n_abs(d, T, ye, yp, yn, e_nu, mu_e, mu_hat);
		std::tie(em_anue,ab_anue) = nu_p_abs(d, T, ye, yp, yn, e_nu_bar, mu_e, mu_hat);

		// neutrino-nucleon scattering
		B_IS_nue  = nu_N_scatt_tot(e_nu, d, T, yn, yp);
                B_IS_anue = anu_N_scatt_tot(e_nu_bar, d, T, yn, yp);

		// write data in myfile
		fout << std::scientific << zone << "\t"  << d << "\t" << T << "\t" << ye << "\t" << mu_e << "\t" << mu_hat << "\t";
		fout << yp << "\t" << yn << "\t" << em_nue << "\t" << ab_nue << "\t" << em_anue << "\t" << ab_anue << "\t";
		fout << B_IS_nue << "\t" << B_IS_anue << "\n";
		
		std::getline(fin, dummyLine);
		// comparison with Fortran code
		//nuscatt_mod = 2.*pi*nuscatt;
	}

	fin.close();
	fout.close();

}
