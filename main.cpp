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

#include "./source/constants.h" //Header file containing all relevant constants
#include "./source/parameters.h" //Header file containing all parameters
#include "./source/weak_magnetism.h" //Header file containing functions to compute WM
#include "./source/nu_abs_em.h" //Header file containing all rates
#include "./source/nu_elastic_scatt.h" //Header file containing elastic scattering rates

using namespace constants;
using namespace parameters;
using namespace nuabsem;
using namespace elastic_scatt;
using namespace weakmag;
using namespace std;

int main (){
	int zone;
	double r, d, T;
	double nb;
	double yh, ya, ye, yn, yp;
	double dU;
	double R, Rbarc, Rn, Rbarn, Rp, Rbarp;
	double E[ne];
	double e_nu = 10.079454858851342;
	double e_nu_bar = e_nu;
	double mu_e, mu_hat;
	double np, nn;
	double em_nue, ab_nue, em_anue, ab_anue, B_IS_nue, B_IS_anue;
	string dummyLine;
	
	// define energy array
	//double f1 = pow(emax/emin)
      
	//	llocate(e(ne),de(ne))

      //f1 = (emax/emin)**(1./float(ne-1))
      //f2 = (f1-1.)/sqrt( (1.+f1*(1.+f1))/3. )
      //e(1) = emin
      //de(1) = f2*e(1)
      //do ie=2,ne
        //e(ie) = f1*e(ie-1)
        //de(ie) = f2*e(ie)
      //enddo
      //egroup_empty = .false.

	// open input profile
	//string filename = "nurates_1.008E+01.txt";
	std::fstream fin("./input/nurates_1.008E+01.txt"); //read file
        if (!fin)
        {
                cout << "File not found" << '\n';
        }
	for (int i = 0; i < 3; i++) { //skip first two lines and read first one with data
		getline(fin, dummyLine);
	}
	
	// open output file
	//string outname;
	char outname[35];
	//if ((use_WM_ab == 1) && (use_WM_sc == 1)) {
	if (use_dU == 1) {
		sprintf(outname, "./input/newrates_%.3e_dU.txt", e_nu);
 	} else {
		sprintf(outname, "./input/newrates_%.3e.txt", e_nu);
	}
	ofstream fout(outname);
	//myfile.open("check_rates.txt");
	// file header
	string fhead;
	fhead = "# id, d [g/cm^3], T [MeV], Ye, mu_e [MeV], mu_hat [MeV], Yp, Yn, dU, em_nue, ab_nue[1/cm], em_anue[1/cm], ab_anue[1/cm], R, Rbar, B_IS_nue[1/cm], B_IS_anue[1/cm], Rn, Rbarn, Rp, Rbarp\n";
	fout << fhead;

        R    = WM_nue_abs(e_nu);
	Rbarc = WM_anue_abs(e_nu_bar);
	std::tie(Rp,Rbarp) = WM_scatt(e_nu, 1);
        std::tie(Rn,Rbarn) = WM_scatt(e_nu, 2);


	while (!fin.fail()){
		// read data from input table
		std::stringstream ss(dummyLine);
		ss >> zone >> r >> d >> T >> ye >> mu_e >> mu_hat >> yh >> ya >> yp >> yn >> dU;

		nb = d/mb;

		//if (use_dU == 0) dU = 0.;

		// electron (anti)neutrino absorption 
		std::tie(em_nue,ab_nue)   = nu_n_abs(e_nu, nb, T, me, yp, yn, mu_e, mu_hat, dU);
		std::tie(em_anue,ab_anue) = nu_p_abs(e_nu_bar, nb, T, me, yp, yn, mu_e, mu_hat, dU);

		// neutrino-nucleon scattering
		B_IS_nue  = nu_N_scatt_tot(e_nu, nb, T, yn, yp);
                B_IS_anue = anu_N_scatt_tot(e_nu_bar, nb, T, yn, yp);

		// write data in myfile
		fout << std::scientific << zone << "\t"  << d << "\t" << T << "\t" << ye << "\t" << mu_e << "\t" << mu_hat << "\t";
		fout << yp << "\t" << yn << "\t" << dU << "\t" << em_nue << "\t" << ab_nue << "\t" << em_anue << "\t" << ab_anue << "\t";
		fout << R << "\t" << Rbarc << "\t";
		fout << B_IS_nue << "\t" << B_IS_anue << "\t";
		fout << Rn << "\t" << Rbarn << "\t" << Rp << "\t" << Rbarp << "\n";
		
		std::getline(fin, dummyLine);
		// comparison with Fortran code
		//nuscatt_mod = 2.*pi*nuscatt;
	}

	fin.close();
	fout.close();

}
