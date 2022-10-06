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
	double T = 10., ne = 1.e35;
	double eta1 = -5.e3; //-1.3e4; //minimum non relat degeneracy parameter
        double eta2 =  5.e3; //1.3e4; //maximum non relat degeneracy parameter

	double h1 = compute_res(eta, theta, k, w_leg, w_lag);
	double h2 = compute_res_ed(eta, theta, k, w_leg, w_lag);
	//printf("%.15e, %.15e\n", h1, h2);

	struct eta_NR_params p = {ne, theta, me, w_leg, w_lag};

	FDF.f = &n_GSL_f;
        FDF.df = &n_GSL_df;
        FDF.fdf = &n_GSL_fdf;
        FDF.params = &p;

	double yemin = 5.e-3;
        double yemax = 5.5e-1;
        double denmin = 5.e+3;
        double denmax = 2.e+16;

        double ne_min = yemin*denmin/mb;
        double ne_max = yemax*denmax/mb;
        double t_min = 5.e-02;
        double t_max = 1.05e+2;

	int nne=1500;
        int nt=340;

	double result;
	double lne_edges[nne+1];
	double lt_edges[nt+1];
	double ne_array[nne];
	double t_array[nt];
	//.......compute the grid
        for (int i=0;i<nne+1;i++) lne_edges[i] = log10(ne_min)+ ((double) (i-1)) *(log10(ne_max)-log10(ne_min))/( (double) nne);
        for (int i=0;i<nne;i++) ne_array[i] = pow(10.,(0.5*(lne_edges[i+1]+lne_edges[i])));

        for (int j=0;j<nt+1;j++) lt_edges[j] = log10(t_min)+ ((double) (j-1)) *(log10(t_max)-log10(t_min))/((double) nt);
        for (int j=0;j<nt;j++) t_array[j] = pow(10.,(0.5*(lt_edges[j+1]+lt_edges[j])));

        for (int i=0;i<nne;i++) {
		for (int j=0;j<nt;j++) {
			ne = ne_array[i];
			T = t_array[j];
			guess = find_guess_e(ne,T);
			printf("guess = %.3e\n", guess);
			double result = rtsafe(ne, T, me, guess, eta1, eta2, w_leg, w_lag);
			printf("eta = %.6e\n\n", result);
		}
	}
	//double result = NRaphson_1D_GSL(max_iter, guess,  FDF);

	gsl_integration_fixed_free (w_leg);
	gsl_integration_fixed_free (w_lag);
	return 0;
}

