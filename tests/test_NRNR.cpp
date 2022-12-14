#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <string> 
#include <vector>
#include <sstream>
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include <cstdlib> //libabry for abs()
#include <functional> 

#include "../source/constants.h" //Header file containing all relevant constants
#include "../source/parameters.h" //Header file containing all parameters
#include "../source/weak_magnetism.h" //Header file containing functions to compute WM
#include "../source/nu_abs_em.h" //Header file containing all rates
#include "../source/nu_elastic_scatt.h" //Header file containing elastic scattering rates
#include "../source/spectral_function.h"
#include "../source/tools/fermi_integrals.h"
#include "../source/tools/NewtonRaphson.h"
#include "../source/mnewt.h"
#include "../source/tools/nr3.h"

using namespace constants;
using namespace parameters;
using namespace nuabsem;
using namespace elastic_scatt;
using namespace weakmag;
using namespace std;

int main (){
	double mu_e = 2.495089e+02;
	double mu_hat = 6.015553e+01;
        double mu_nu = mu_e - (Q+mu_hat);
	double T = 1.204000E+01;
	double A = 1./T;
	double B = mu_nu/T;
	double n_nu, j_nu;
	int opt = 1; //0: optically thick, 1: optically thin

        if (opt == 0) { //optically thick
                 n_nu = 4.*pi*Fermi_integral_p2(B)/pow(A,3.);
                 j_nu = 4.*pi*Fermi_integral_p3(B)/pow(A,4.);
        } else if (opt == 1) { //optically thin
                 n_nu = 4.*pi*gsl_sf_gamma(A+3.)/pow(B,A+3.);
                 j_nu = 4.*pi*gsl_sf_gamma(A+4.)/pow(B,A+4.);
        }


	gsl_multiroot_function_fdf FDF;
        gsl_multiroot_function F;
	
	double x[2], C[2] = {n_nu/pow(h*c,3.),j_nu/pow(h*c,3.)};
	x[0] = A * 0.95;
	x[1] = B * 1.05;
	printf("%.3lf, %.3lf\n", A, B);
	
	printf("NR method:\n");
        if (opt == 0) { //optically thick
                 mnewt2d(x,C,&thick_fdf);
        } else if (opt == 1) { //optically thin
                 mnewt2d(x,C,&thin_fdf);
        }

	printf("%.3lf, %.3lf\n", x[0], x[1]);

  
	printf("GSL method:\n");
	struct M1_values p = {n_nu/pow(h*c,3.),j_nu/pow(h*c,3.)};

	if (opt == 0) {
		FDF = {&thick_f, &thick_df, &thick_fdf_GSL, 2, &p};
		F = {&thick_f, 2, &p};
	} else if (opt == 1) {
		FDF = {&thin_f, &thin_df, &thin_fdf_GSL, 2, &p};
		F = {&thin_f, 2, &p};
	}

	double x0[2];
	double xf[2];
	
	x0[0] = A*0.95;
	x0[1] = B*1.05;

	NRaphson_2D(x0, xf, FDF);
	printf("%.3lf, %.3lf\n", xf[0], xf[1]);

	return 0;
}

