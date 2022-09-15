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

#include "./source/constants.h" //Header file containing all relevant constants
#include "./source/parameters.h" //Header file containing all parameters
#include "./source/weak_magnetism.h" //Header file containing functions to compute WM
#include "./source/nu_abs_em.h" //Header file containing all rates
#include "./source/nu_elastic_scatt.h" //Header file containing elastic scattering rates
#include "./source/integration.h" //Header file contatining integration routines
#include "./source/tools/nr3.h"

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

double j_leg(double x, void * p) {
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


double eta_leg(int n, gsl_function F, double t) {
        double x[n], w[n];
        double f1[n], f2[n];

        read_wgts(n, x, w, 0);

        void *p = F.params;

        for (int i=0;i<n;i++) {
                f1[i] = j_leg(t*x[i], p);
                f2[i] = j_leg(t/x[i], p) / (x[i]*x[i]);
        }

        return t * (do_integration(n, w, f1) + do_integration(n, w, f2));
}


double GSL_integration(int n, double a, double b, double alpha, double beta, gsl_function F, const gsl_integration_fixed_type * gauss_type) {
        gsl_integration_fixed_workspace *w;
        double result;
        w = gsl_integration_fixed_alloc(gauss_type, n, a, b, alpha, beta);
        gsl_integration_fixed(&F, &result, w);
        gsl_integration_fixed_free (w);
        return result;
}



int main (){
	int n = 24;
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
	double alpha = 5., beta = 0.;
	double GSL_lag, NR_lag, NR_leg;

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

	string Ioutname = "./output/j_integral.txt";
	ofstream Iout(Ioutname);
	Iout << "# r [cm], eta_e, Laguerre [GSL], Laguerre [NR], Legendre [NR]\n";


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

		double eta_e = mu_e/T;
		double t = eta_e;
		
		gsl_function F_lag;		
		F_lag.function = &j_lag;
		F_lag.params = &params;
		const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;
		GSL_lag = GSL_integration(n, a, b, alpha, beta, F_lag, lag_type);
		
	        NR_lag = eta_lag(n, F_lag);

        	gsl_function F_leg;
		F_leg.function = &j_leg;
		F_leg.params = &params;

		NR_leg = eta_leg(n, F_leg, t);

		printf("%.3e, %.3e, %.3e\n", GSL_lag, NR_lag, NR_leg);

		Iout << std::scientific << r << "\t" << eta_e << "\t" << GSL_lag << "\t" << NR_lag << "\t" << NR_leg << "\n"; 
			
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

