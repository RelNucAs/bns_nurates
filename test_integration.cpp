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

#include "./source/integration.h" //Header file contatining integration routines
#include "./source/tools/nr3.h"
#include "./source/tools/fermi_integrals.h"


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

double f_lag(int n, gsl_function F) {
	double x[n], w[n];
	double f[n];

	read_wgts(n, x, w, 1);

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


double GSL_integration(int n, double a, double b, double alpha, double beta, gsl_function F, const gsl_integration_fixed_type * gauss_type) {
	gsl_integration_fixed_workspace *w;
        double result;
        w = gsl_integration_fixed_alloc(gauss_type, n, a, b, alpha, beta);
        gsl_integration_fixed(&F, &result, w);
        gsl_integration_fixed_free (w);
	return result;
}

int main (){
	int n = 32;
	double a = 0., b = 1.;
	double alpha = 5., beta = 0.;
	double result, GSL_lag, NR_lag, NR_leg_old;
	double NR_leg[3];

	int k = 5.;
	gsl_function F_lag, F_leg;
        const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;
	double t;

	int nstep = 500;
	double emin = 0.1, emax = 300.;
	double de = (log10(emax)-log10(emin)) / (double) nstep;
	double eta[nstep];

	for (int i=0;i<nstep;i++) {
		eta[i] = pow(10., log10(emin)+i*de);
	}

	save_gauleg(n, a, b);
	if ((n==32) && (alpha==5.)) {
		save_gaulag(24, alpha);
	} else {
		save_gaulag(n, alpha);
	}

        char Ioutname[50];
        sprintf(Ioutname, "./output/FD_integral_alpha_%.0lf_n_%d.txt", alpha, n);
	ofstream Iout(Ioutname);
	Iout << "# eta, Exact result, Laguerre [GSL], Laguerre [NR], Legendre [NR]\n";


	for (int i=0;i<nstep;i++) {
		struct my_f_params params = {k,eta[i],alpha};
		t = eta[i];

		result = Fermi_integral_p5(eta[i]);

		//double F_lag, F_leg;

		F_lag.function = &FD_laguerre;
		F_lag.params = &params;

		GSL_lag = GSL_integration(n, a, b, alpha, beta, F_lag, lag_type);
		
		if ((n==32) && (alpha==5.)) {
			NR_lag = GSL_lag;
		} else {
			NR_lag = f_lag(n, F_lag);
		}

		F_leg.function = &FD_legendre;
		F_leg.params = &params;

		NR_leg[0] = f_leg(n, F_leg, 0.5*t);
		NR_leg[1] = f_leg(n, F_leg, t);
		NR_leg[2] = f_leg(n, F_leg, 2.0*t);

        	//f_leg_old(n, leg_arr1, leg_arr2, F_leg, t);
        	//NR_leg_old = compute_split_gauleg(n, leg_arr1, leg_arr2, t);

		printf("Exact result: %.6e, GSL_lag = %.6e, NR_lag = %.6e, NR_leg = %.6e\n", result, GSL_lag, NR_lag, NR_leg[1]);
		Iout << std::scientific << eta[i] << "\t" << result << "\t" << GSL_lag << "\t" << NR_lag << "\t"; 
                Iout << NR_leg[0] << "\t" << NR_leg[1] << "\t" << NR_leg[2] << "\n"; 
	}

	Iout.close();

	return 0;
}

