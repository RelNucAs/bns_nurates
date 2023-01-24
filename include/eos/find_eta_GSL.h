#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "../constants_bis.h"
#include "../tools/nr3.h"
#include "complete_FG.h"

using namespace constants;
struct eta_NR_params { double nLep; double T; double mLep; gsl_integration_fixed_workspace *w_leg; gsl_integration_fixed_workspace *w_lag; };

//.......Given the density rho, the non relativistic degeneracy eta, the temperature T
//.......and the particle species, the subroutine computes the net fraction ynet and the first eta derivative.
//.......UNITS: rho in g/cm**3; T in MeV
double n_net_f(double eta, double T, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        double f1, f2, f3, f4;
        double theta = T/mLep;

        f1 = compute_res( eta,          theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 0.5);
        f3 = compute_res( eta,          theta, 1.5);
        f4 = compute_res(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}


double n_net_df(double eta, double T, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        double f1, f2, f3, f4;
        double theta = T/mLep;

        f1 = compute_res_ed( eta,          theta, 0.5);
        f2 = compute_res_ed(-eta-2./theta, theta, 0.5);
        f3 = compute_res_ed( eta,          theta, 1.5);
        f4 = compute_res_ed(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}


double find_guess_e(double ne, double T) {
	double nq, guess;
	double t = T/kB; //temp in K
	double nfm = ne*1.e-39; //number density in fm^{-3}
	double pF = h*c*pow(3.*ne/(8.*pi),1./3.); //[MeV]

//.........Find the guess
        if ((t<=1.e10) && (nfm <= 1.e-17)) { //classical NR
          	nq = pow(2.*pi*me*T/(h*h*c*c),1.5); //cm^{-3}
		guess = -log(fabs(2.*nq/ne));
	} else if ((t<=1.e10) && (nfm>1.e-17) && (nfm <= 1.e-10)) { //degenerate NR
        	guess = (pF*pF/(2.*me*T)) - (me/T);       
        } else if ((t<=1.e10) && (nfm>1.e-10)) { //degenerate UR
        	guess = (pF/T) - (me/T);
        } else if (t>1.e10) {
        	nq = 8.*pi*pow(T/(h*c),3.); //classical UR
        	guess = -log(fabs(2.*nq/ne)) - (me/T);
        }
	return guess;
}

double find_guess_mu(double nmu, double T) {
        double nq, guess;
//.........Find the guess
	nq = pow(2.*pi*mmu*T/(h*h*c*c),1.5); //cm^{-3}
	guess = -log(fabs(2.*nq/nmu));
        return guess;
}


double NRaphson_1D_GSL(int max_iter, double guess, gsl_function_fdf FDF) {
        int status;
        int iter = 0;
        double x0, x = guess;
        const gsl_root_fdfsolver_type *T;
        gsl_root_fdfsolver *s;
        T = gsl_root_fdfsolver_newton;
        s = gsl_root_fdfsolver_alloc (T);
        gsl_root_fdfsolver_set (s, &FDF, x);

        do {
                iter++;
                status = gsl_root_fdfsolver_iterate (s);
                x0 = x;
                x = gsl_root_fdfsolver_root (s);
                status = gsl_root_test_delta (x, x0, 0, 1e-3);

                //if (status == GSL_SUCCESS)
                        //printf ("Converged:\n");

        }
        while (status == GSL_CONTINUE && iter < max_iter);
        gsl_root_fdfsolver_free(s);
        return x;
}
