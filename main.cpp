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

#include "constants.hpp" //Header file containing all relevant constants
#include "parameters.hpp" //Header file containing all parameters
#include "weak_magnetism.hpp" //Header file containing functions to compute WM
#include "nu_abs_em.hpp" //Header file containing all rates
#include "nu_elastic_scatt.hpp" //Header file containing elastic scattering rates
#include "tools/integration.hpp" //Header file contatining integration routines
#include "spectral_function.hpp" //Header file contatining spectral functions
#include "tools/nr3.h"
//#include "tools/NewtonRaphson.h"
#include "tools/fermi_integrals.h"

using namespace constants;
using namespace parameters;
using namespace nuabsem;
using namespace elastic_scatt;
using namespace std;

//struct my_f_params { double nb; double T; double me; double yp; double yn; double mu_e; double mu_hat; double dU; };

int main (){
	int n = 32;
	int zone;
	double r, d, T;
	double nb;
	double yh, ya, ye, yn, yp;
	double dU;
	double R, Rbar, R0_n, R1_n, R0_p, R1_p;
	double E[ne];
	double e_nu = 3.; //10.079454858851342;
	double e_nu_bar = e_nu;
	double mu_e, mu_hat;
	double np, nn;
	double em_nue, ab_nue, em_anue, ab_anue, B_IS_nu;
	string dummyLine;

	double a = 0., b = 1.;
	double alpha = 0., beta = 0.;
	//double GSL_lag, NR_lag;
	my_function F1, F2, F3, F4;
	double eta_nue, eta_anue;
	double eta0_nue, eta0_anue;
	double t1, t2, t3, t4;

	save_gauleg(n, a, b);
	save_gaulag(n, alpha);

	F1.function = &j_nue;
	F2.function = &j_anue;
	F3.function = &j0_nue;
	F4.function = &j0_anue;
	
	// open input profile
	std::fstream fin("./input/nurates_1.008E+01.txt"); //read file
        if (!fin)
        {
                cout << "File not found" << '\n';
        }
	for (int i = 0; i < 3; i++) { //skip first two lines and read first one with data
		getline(fin, dummyLine);
	}
	
	// open output file
	char outname[35];
	//if ((use_WM_ab == 1) && (use_WM_sc == 1)) {
	if (use_dU == 1) {
		sprintf(outname, "./input/newrates_%.3e_dU.txt", e_nu);
 	} else {
		sprintf(outname, "./input/newrates_%.3e.txt", e_nu);
	}
	ofstream fout(outname);

	string Ioutname = "./output/j_integral_2";
	if (use_WM_ab == 1) Ioutname = Ioutname + "_WM";
	if (use_dU == 1)    Ioutname = Ioutname + "_dU";
	Ioutname = Ioutname + ".txt";
	ofstream Iout(Ioutname);
	Iout << "# r [cm], T [MeV], mu_e [MeV], [Laguerre] n=8, n=16, n=24, n=32, n=100\n";


	// file header
	string fhead;
	fhead = "# id, d [g/cm^3], T [MeV], Ye, mu_e [MeV], mu_hat [MeV], Yp, Yn, dU, em_nue [1/s], ab_nue [1/cm], em_anue [1/s], ab_anue [1/cm], R, Rbar, B_IS_nu [1/cm], R0_n, R1_n, R0_p, R1_p\n";
	fout << fhead;

        R    = WM_nue_abs(e_nu);
	Rbar = WM_anue_abs(e_nu_bar);
	std::tie(R0_p,R1_p) = WM_scatt(e_nu, 1);
        std::tie(R0_n,R1_n) = WM_scatt(e_nu, 2);

	while (!fin.fail()){
		// read data from input table
		std::stringstream ss(dummyLine);
		ss >> zone >> r >> d >> T >> ye >> mu_e >> mu_hat >> yh >> ya >> yp >> yn >> dU;

		nb = d/mb;

		// electron (anti)neutrino absorption 
		std::tie(em_nue,ab_nue)   = nu_n_abs(e_nu, nb, T, me, yp, yn, mu_e, mu_hat, dU);
		std::tie(em_anue,ab_anue) = nu_p_abs(e_nu_bar, nb, T, me, yp, yn, mu_e, mu_hat, dU);

		// neutrino-nucleon scattering
		B_IS_nu  = nu_N_scatt_tot(e_nu, nb, T, yn, yp);
		
		struct my_f_params params = {nb, T, me, yp, yn, mu_e, mu_hat, dU};

		//double eta_e = mu_e/T;

		
		//gsl_function F_lag;		
		//F_lag.function = &j_lag;
		//F_lag.params = &params;
		//const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;
		//GSL_lag = GSL_integration(n, a, b, alpha, beta, F_lag, lag_type);
		
	        //NR_lag = eta_lag(n, F_lag);

        	//gsl_function F_leg;
		//F_leg.function = &j_leg;
		//F_leg.params = &params1

                //x = x_NRaphson(max_iter, x0, FDF);
		
		F1.params = &params;
		F2.params = &params;
		F3.params = &params;
		F4.params = &params;
		
		t1 = max(mu_e-Q,5.);
		//t = find_max(1000, 0.1, 400., &j_nue, params);
		t2 = max(-(mu_e-Q),5.); 
		t3 = max(mu_e-Q,4.); 
		t4 = max(-(mu_e-Q),4.); 
		
		Iout << std::scientific << setprecision(10) << r << "\t" << T << "\t" << mu_e << "\t"; 

		for (int nn=8;nn<=32;nn+=8) {
                        eta_nue   = gleg_integ(nn, F1, t1);
                        eta_anue  = gleg_integ(nn, F2, t2);
                        eta0_nue  = gleg_integ(nn, F3, t3);
                        eta0_anue = gleg_integ(nn, F4, t4);

                        Iout << eta_nue << "\t"; 

			printf("t1 = %.3e, eta_nue = %.6e;	t2 = %.3e, eta_anue  = %.6e\n", t1, eta_nue, t2, eta_anue);
			printf("t3 = %.3e, eta0_nue = %.6e;	t4 = %.3e, eta0_anue = %.6e\n\n", t3, eta0_nue, t4, eta0_anue);
		}

		eta_nue   = gleg_integ(100, F1, t1);
		eta_anue  = gleg_integ(100, F2, t2);
		eta0_nue  = gleg_integ(100, F3, t3);
		eta0_anue = gleg_integ(100, F4, t4);

		Iout << eta_nue << "\n";


		//printf("%.3e, %.3e\n", GSL_lag, NR_lag, NR_leg);
		//Iout << std::scientific << r << "\t" << T << "\t" << mu_e << "\t" << GSL_lag << "\t" << NR_lag << "\t" << NR_leg << "\n"; 
			
		// write data in myfile
		fout << std::scientific << zone << "\t"  << d << "\t" << T << "\t" << ye << "\t" << mu_e << "\t" << mu_hat << "\t";
		fout << yp << "\t" << yn << "\t" << dU << "\t" << em_nue << "\t" << ab_nue << "\t" << em_anue << "\t" << ab_anue << "\t";
		fout << R << "\t" << Rbar << "\t";
		fout << B_IS_nu << "\t" <<  R0_n << "\t" << R1_n << "\t" << R0_p << "\t" << R1_p << "\n";
		
		std::getline(fin, dummyLine);
		// comparison with Fortran code
		//nuscatt_mod = 2.*pi*nuscatt;
	}

	
	double k = c * (GF*GF/pi) / pow(h*c/(2*pi),4) * (gV*gV+3*gA*gA);
	
	fin.close();
	fout.close();
	Iout.close();

	return 0;
}

