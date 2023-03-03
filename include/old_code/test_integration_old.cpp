// Test integration routines

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

#include "tools/integration.hpp" //Header file contatining integration routines
#include "integration_GSL.h" //Header file contatining integration routines
#include "tools/nr3.h"
#include "tools/fermi_integrals.h"
#include "tools/NewtonRaphson_GSL.h"
//struct FD_params { int k; double eta; double alf; };

int main (){
	int n = 32;
	double a = 0., b = 1.;
	double beta = 0., alf = 5.;
	double result;
	double GSL_lag[2], NR_lag[2], NR_leg[3];

	int k = 5.;
	gsl_function F_lag, F_leg;
	my_function FD;
	gsl_function_fdf FDF;
        const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;
	double t;

	int max_iter = 100;
	double x, x0;

	int nstep = 500;
	double emin = 0.1, emax = 120.;
	double de = (log10(emax)-log10(emin)) / (double) nstep;
	double eta[2*nstep];

	for (int i=0;i<nstep;i++) {
		eta[nstep+i] = pow(10., log10(emin)+i*de);
		eta[(nstep-1)-i] = -eta[nstep+i];
	}

	save_gauleg(n, a, b);
	save_gauleg(100, a, b);
	save_gaulag(n, 0.);

	if (n<32) {
		save_gaulag(n, 5.);
	}

        char Ioutname[50];
        sprintf(Ioutname, "../output/FD_integral_n_%d.txt", n);
	ofstream Iout(Ioutname);
	Iout << "# eta, t_split, Exact result, Lag (alf=0) [GSL], Lag (alpha=5) [GSL], Lag (alpha=0) [NR], Lag (alpha=5) [NR], Legendre [NR]\n";
	
	string IIoutname = "../output/FD_integral_n_100.txt";
	ofstream IIout(IIoutname);
	IIout << "# eta, t_split, Exact result, [Legendre] t=eta/2, t=eta, t=2eta, t=x\n";

	for (int i=0;i<2*nstep;i++) {
		result = Fermi_integral_p5(eta[i]);

		// Gauss-Laguerre
		//alf = 0
		alf = 0.;
		struct FD_params params = {k,eta[i],alf};

		F_lag.function = &FD_laguerre;
		F_lag.params = &params;

		FD.function = &FD_legendre;
		FD.params = &params;

		GSL_lag[0] = GSL_integration(n, a, b, alf, beta, F_lag, lag_type);
		
		NR_lag[0] = glag_integ(n, FD, alf);

		//alf = 5
		alf = 5.;
		params.alpha = alf;

		//F_lag.function = &FD_laguerre;
		F_lag.params = &params;
		FD.params = &params;

		GSL_lag[1] = GSL_integration(n, a, b, alf, beta, F_lag, lag_type);
		
		if (n<32) {
			NR_lag[1] = glag_integ(n, FD, alf);
		} else {
			NR_lag[1] = GSL_lag[1];
		}


		// Gauss-Legendre
		F_leg.function = &FD_legendre;
		F_leg.params = &params;

		FDF.f = &f_NRaphson; 
		FDF.df = &df_NRaphson;
		FDF.fdf = &fdf_NRaphson;
		FDF.params = &params;
		
		x0 = find_guess(1000, FD);

		//printf("id = %d, scanned_value = %.3lf\n", i, x);	

		x = x_NRaphson(max_iter, x0, FDF);
		
		x = max(5., eta[i]);

		t = eta[i];

		NR_leg[0] = gleg_integ(n, FD, 0.5*t);
		NR_leg[1] = gleg_integ(n, FD, t);
		NR_leg[2] = gleg_integ(n, FD, 2.0*t);
		NR_leg[3] = gleg_integ(n, FD, x);

		//double ciao = gleg_integ(100, F_leg, x);
		//printf("%.6e, %.6e, %.10e, %.10e, %.10e, %.10e\n", eta[i], x, gleg_integ(100, F_leg, 0.5*t), f_leg(100, F_leg, t), f_leg(100, F_leg, 2.0*t), f_leg(100, F_leg, x));
		printf("Eta = %.3e, x = %.3e, Exact result: %.6e, Lag_0 = %.6e, Lag_5 = %.6e, Leg = %.6e\n", eta[i], x, result, GSL_lag[0], NR_lag[1], NR_leg[3]);
		//printf("GSL_Lag_0 = %.6e, NR_Lag_0 = %.6e\n", GSL_lag[0], NR_lag[0]);
		//printf("GSL_Lag_5 = %.6e, NR_Lag_5 = %.6e\n", GSL_lag[1], NR_lag[1]);
		Iout << std::scientific << setprecision(10) << eta[i] << "\t" << x << "\t" << result << "\t";
	        Iout << GSL_lag[0] << "\t" << GSL_lag[1] << "\t";
		Iout << NR_lag[0]  << "\t" << NR_lag[1]  << "\t"; 
                Iout << NR_leg[0] << "\t" << NR_leg[1] << "\t" << NR_leg[2] << "\t" << NR_leg[3] << "\n"; 

		IIout << std::scientific << setprecision(10) << eta[i] << "\t" << x << "\t" << result << "\t";
		IIout << gleg_integ(100, FD, 0.5*t) << "\t";
		IIout << gleg_integ(100, FD, t) << "\t";
		IIout << gleg_integ(100, FD, 2.0*t) << "\t";
		IIout << gleg_integ(100, FD, x) << "\n";

	}

        struct FD_params pars = {k,100.,5.};
	FD.params = &pars;
        x0 = find_guess(1000, FD);
        printf("scanned_value = %.3lf\n", x_NRaphson(max_iter, x0, FDF));


	Iout.close();
	IIout.close();


	return 0;
}

