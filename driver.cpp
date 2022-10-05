#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "./source/eos/complete_FG.h"
#include "./source/constants.h"
#include "./source/parameters.h"
#include "./source/eos/find_eta.h"
#include "./source/tools/nr3.h"

int main (){
        int ngle = 20;
        int ngla = 20;
        double a = 0., b = 1.;
        double beta = 0.;
	int max_iter = 150;
	gsl_function_fdf FDF;

        const gsl_integration_fixed_type * leg_type = gsl_integration_fixed_legendre;
        const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;

        gsl_integration_fixed_workspace *w_lag, *w_leg;

        w_lag = gsl_integration_fixed_alloc(lag_type, ngla, a, b, alpha, beta);
        w_leg = gsl_integration_fixed_alloc(leg_type, ngle, a, b, alpha, beta);

	double eta = 20.;
	float k = 1.5;
	double theta = 20.;
	double guess;
	double T = 10., ne = 1.e10/mb;
	double eta1 = -1.3e4; //minimum non relat degeneracy parameter
        double eta2 =  1.3e4; //maximum non relat degeneracy parameter

	double h1 = compute_res(eta, theta, k, w_leg, w_lag);
	double h2 = compute_res_ed(eta, theta, k, w_leg, w_lag);
	//printf("%.15e, %.15e\n", h1, h2);

	struct eta_NR_params p = {ne, theta, me, w_leg, w_lag};

	FDF.f = &n_GSL_f;
        FDF.df = &n_GSL_df;
        FDF.fdf = &n_GSL_fdf;
        FDF.params = &p;

	guess = find_guess_e(ne,T);
	printf("guess = %.3e\n", guess);
	double result = rtsafe(ne, T, me, guess, eta1, eta2, w_leg, w_lag);
	//double result = NRaphson_1D_GSL(max_iter, guess,  FDF);
	printf("eta = %.6e\n", result);
	gsl_integration_fixed_free (w_leg);
	gsl_integration_fixed_free (w_lag);
	return 0;
}

