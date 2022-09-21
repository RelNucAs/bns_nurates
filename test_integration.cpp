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

#include "./source/integration.h" //Header file contatining integration routines
#include "./source/tools/nr3.h"
#include "./source/tools/fermi_integrals.h"
#include "./source/tools/NewtonRaphson.h"

struct my_f_params { int k; double eta; double alpha; };

double FD_laguerre(double x, void *p) {
	struct my_f_params * params = (struct my_f_params *)p;
	int k = (params->k);
	double eta = (params->eta);
	double alpha = (params->alpha);
	return exp(x) * pow(x,-alpha) * pow(x,k) / (exp(x-eta)+1.);
}

double FD_legendre(double x, void *p) {
	struct my_f_params * params = (struct my_f_params *)p;
	int k = (params->k);
	double eta = (params->eta);
	return pow(x,k) / (exp(x-eta)+1.);
}

double f_NRaphson(double x, void *p) {
	struct my_f_params * params = (struct my_f_params *)p;
	double eta = (params->eta);
	return 5.*exp(eta) - exp(x)*(x-5.);
}

double df_NRaphson(double x, void *p) {
	struct my_f_params * params = (struct my_f_params *)p;
	double eta = (params->eta);
	return exp(x)*(4.-x);
}

void fdf_NRaphson(double x, void *p, double * f, double * df) {
	struct my_f_params * params = (struct my_f_params *)p;
	double eta = (params->eta);
   	double t  = 5.*exp(eta) - exp(x)*(x-5.);
	double dt = exp(x)*(4.-x); 
	*f = t;
	*df = dt;   /* uses existing value */
}

double find_guess(int nslice, gsl_function F) {
	double val, max = 0;
	double guess;
	void *p = F.params;
        struct my_f_params * params = (struct my_f_params *)p;
        double eta = (params->eta);
	double x0, x1;

	if (abs(eta) < 10.) {
		x0 = 0.1*abs(eta),
		x1 = 100*abs(eta);
	} else {
		x0 = 0.5*abs(eta);
		x1 = 2.0*abs(eta);
	}
       	
	double dx = (x1-x0) / (double) nslice;
	for (double x=x0;x<x1;x+=dx) {
		val = GSL_FN_EVAL(&F,x);
		//printf("x = %.3e, Val = %.3e\n", x, val);
		if (val > max) {
			max = val;
			guess = x;
		}
	}
	return guess;
}
		

double f_lag(int n, gsl_function F, int g_type) {
	double x[n], w[n];
	double f[n];

	read_wgts(n, x, w, g_type);

	void *p = F.params;

	for (int i=0;i<n;i++) f[i] = FD_laguerre(x[i], p);

        return do_integration(n, w, f);
}

double f_leg(int n, gsl_function F, double t) {
	double x[n], w[n];
	double f1[n], f2[n];

	read_wgts(n, x, w, 0);

	void *p = F.params;

	for (int i=0;i<n;i++) {
		f1[i] = FD_legendre(t*x[i], p);
       		f2[i] = FD_legendre(t/x[i], p) / (x[i]*x[i]);
	}
	//struct my_f_params * params = (struct my_f_params *)p;
        //double eta = (params->eta);
	//printf("%.6e, %.6e, %.10e, %.10e\n", eta, t, t*do_integration(n, w, f1), t*do_integration(n, w, f2));
	return t * (do_integration(n, w, f1) + do_integration(n, w, f2));
}

void f_leg_old(int n, double *f1, double *f2, gsl_function F, double t) {
        double x[n], w[n];

        read_wgts(n, x, w, 0);

        void *p = F.params;

        for (int i=0;i<n;i++) {
                f1[i] = FD_legendre(t*x[i], p);
                f2[i] = FD_legendre(t/x[i], p) / (x[i]*x[i]);
        }

	return;
}


int main (){
	int n = 32;
	double a = 0., b = 1.;
	double beta = 0.; //alpha = 5.;
	double result;
	double GSL_lag[2], NR_lag[2], NR_leg[3];

	int k = 5.;
	gsl_function F_lag, F_leg;
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
        sprintf(Ioutname, "./output/FD_integral_n_%d.txt", n);
	ofstream Iout(Ioutname);
	Iout << "# eta, t_split, Exact result, Lag (alpha=0) [GSL], Lag (alpha=5) [GSL], Lag (alpha=0) [NR], Lag (alpha=5) [NR], Legendre [NR]\n";
	
	string IIoutname = "./output/FD_integral_n_100.txt";
	ofstream IIout(IIoutname);
	IIout << "# eta, t_split, Exact result, [Legendre] t=eta/2, t=eta, t=2eta, t=x\n";

	for (int i=0;i<2*nstep;i++) {
		result = Fermi_integral_p5(eta[i]);

		// Gauss-Laguerre
		//alpha = 0
		struct my_f_params params_0 = {k,eta[i],0.};

		F_lag.function = &FD_laguerre;
		F_lag.params = &params_0;

		GSL_lag[0] = GSL_integration(n, a, b, 0., beta, F_lag, lag_type);
		
		NR_lag[0] = f_lag(n, F_lag, 1);

		//alpha = 5
		struct my_f_params params_1 = {k,eta[i],5.};

		//F_lag.function = &FD_laguerre;
		F_lag.params = &params_1;

		GSL_lag[1] = GSL_integration(n, a, b, 5., beta, F_lag, lag_type);
		
		if (n<32) {
			NR_lag[1] = f_lag(n, F_lag, 2);
		} else {
			NR_lag[1] = GSL_lag[1];
		}

		// Gauss-Legendre
		F_leg.function = &FD_legendre;
		F_leg.params = &params_0;

		FDF.f = &f_NRaphson; 
		FDF.df = &df_NRaphson;
		FDF.fdf = &fdf_NRaphson;
		FDF.params = &params_0;
		
		x0 = find_guess(1000, F_leg);
		//printf("id = %d, scanned_value = %.3lf\n", i, x);	

		x = x_NRaphson(max_iter, x0, FDF);
		
		x = max(5., eta[i]);

		t = eta[i];

		NR_leg[0] = f_leg(n, F_leg, 0.5*t);
		NR_leg[1] = f_leg(n, F_leg, t);
		NR_leg[2] = f_leg(n, F_leg, 2.0*t);
		NR_leg[3] = f_leg(n, F_leg, x);

		//double ciao = f_leg(100, F_leg, x);
		//printf("%.6e, %.6e, %.10e, %.10e, %.10e, %.10e\n", eta[i], x, f_leg(100, F_leg, 0.5*t), f_leg(100, F_leg, t), f_leg(100, F_leg, 2.0*t), f_leg(100, F_leg, x));
		//printf("Eta = %.3e, x = %.3e, Exact result: %.6e, Lag_0 = %.6e, Lag_5 = %.6e, Leg = %.6e\n", eta[i], x, result, NR_lag[0], NR_lag[1], NR_leg[1]);
		Iout << std::scientific << setprecision(10) << eta[i] << "\t" << x << "\t" << result << "\t";
	        Iout << GSL_lag[0] << "\t" << GSL_lag[1] << "\t";
		Iout << NR_lag[0]  << "\t" << NR_lag[1]  << "\t"; 
                Iout << NR_leg[0] << "\t" << NR_leg[1] << "\t" << NR_leg[2] << "\t" << NR_leg[3] << "\n"; 

		IIout << std::scientific << setprecision(10) << eta[i] << "\t" << x << "\t" << result << "\t";
		IIout << f_leg(100, F_leg, 0.5*t) << "\t";
		IIout << f_leg(100, F_leg, t) << "\t";
		IIout << f_leg(100, F_leg, 2.0*t) << "\t";
		IIout << f_leg(100, F_leg, x) << "\n";

	}

        struct my_f_params pars = {k,100.,5.};
        FDF.params = &pars;
	F_leg.params = &pars;
        x0 = find_guess(1000, F_leg);
        //printf("scanned_value = %.3lf\n", x_NRaphson(max_iter, x0, FDF));


	Iout.close();
	IIout.close();

	return 0;
}

