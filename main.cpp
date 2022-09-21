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
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "./source/constants.h" //Header file containing all relevant constants
#include "./source/parameters.h" //Header file containing all parameters
#include "./source/weak_magnetism.h" //Header file containing functions to compute WM
#include "./source/nu_abs_em.h" //Header file containing all rates
#include "./source/nu_elastic_scatt.h" //Header file containing elastic scattering rates
#include "./source/integration.h" //Header file contatining integration routines
#include "./source/tools/nr3.h"
#include "./source/tools/NewtonRaphson.h"
#include "./source/tools/fermi_integrals.h"

using namespace constants;
using namespace parameters;
using namespace nuabsem;
using namespace elastic_scatt;
using namespace weakmag;
using namespace std;

struct my_f_params { double nb; double T; double me; double yp; double yn; double mu_e; double mu_hat; double dU; };

double j_lag(double x, void * p) {
	struct my_f_params * params = (struct my_f_params *)p;
	double nb = (params->nb);
	double T  = (params->T);
	double me = (params->me);
	double yp = (params->yp);
	double yn = (params->yn);
	double mu_e = (params->mu_e);
	double mu_hat = (params->mu_hat);
	double dU = (params->dU);
                          
	double fEm, fAb;
	
	std::tie(fEm,fAb) = nu_n_abs(x, nb, T, me, yp, yn, mu_e, mu_hat, dU);
        
    	return 4. * pi * pow(x,3.) * fEm * exp(x) / pow(x,5.);
}

double j_leg(double x, struct my_f_params params) {
	double nb = (params.nb);
	double T  = (params.T);
	double me = (params.me);
	double yp = (params.yp);
	double yn = (params.yn);
	double mu_e = (params.mu_e);
	double mu_hat = (params.mu_hat);
	double dU = (params.dU);
                          
	double fEm, fAb;
	
	std::tie(fEm,fAb) = nu_n_abs(x, nb, T, me, yp, yn, mu_e, mu_hat, dU);
        
    	return 4. * pi * pow(x,3.) * fEm;
}


double eta_lag(int n, gsl_function F) {
        double x[n], w[n];
	double f[n];

        read_wgts(n, x, w, 1);

        void *p = F.params;

        for (int i=0;i<n;i++) f[i] = j_lag(x[i], p);

        return do_integration(n, w, f);
}


void eta_leg(int n, double *f1, double *f2, struct my_f_params p, double t) {
        double x[n], w[n];

        read_wgts(n, x, w, 0);

        for (int i=0;i<n;i++) {
                f1[i] = j_leg(t*x[i], p);
                f2[i] = j_leg(t/x[i], p);
        }

        return;
}


double find_max(int nslice, struct my_f_params p) {
        double val, max = 0;
        double guess;
	double x0 = log10(0.1), x1 = log10(400.);
        double dx = (x1-x0) / (double) nslice;
        for (double x=x0;x<x1;x+=dx) {
                val = j_leg(pow(10.,x), p);
		//printf("%.3e  ",val);
                if (val > max) {
                        max = val;
                        guess = pow(10.,x);
                }
        }
	printf("\n");
        return guess;
}



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
	double GSL_lag, NR_lag, NR_leg;

	//double f1[n] = {0.}, f2[n] = {0.};
	double t;

	save_gauleg(n, a, b);
	save_gaulag(n, alpha);

	
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

	string Ioutname = "./output/j_integral_1_WM_dU.txt";
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
		//F_leg.params = &params;

  		//x0 = find_guess(1000, F_leg);
                //printf("id = %d, scanned_value = %.3lf\n", i, x);     

                //x = x_NRaphson(max_iter, x0, FDF);

		t = max(mu_e-Q,5.);
		//t = find_max(1000, params);

		Iout << std::scientific << setprecision(10) << r << "\t" << T << "\t" << mu_e << "\t"; 

		for (int nn=8;nn<=32;nn+=8) {
		        double f1[nn] = {0.};
		        double f2[nn] = {0.};

                        eta_leg(nn, f1, f2, params, t);
                        NR_leg = split_gauleg(nn, f1, f2, t);

                        Iout << NR_leg << "\t"; 

			//printf("t = %.3e, max = %.3e, Result = %.6e, Exact = %.6e\n", t, Mscan, NR_leg, ciao);
		}

		double f1[100] = {0.};
		double f2[100] = {0.};

		eta_leg(100, f1, f2, params, t);
		NR_leg = split_gauleg(100, f1, f2, t);

		Iout << NR_leg << "\n";


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

