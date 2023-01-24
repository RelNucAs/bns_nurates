// Test integration routines

#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <string> 
#include <vector>
#include <sstream>
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include <cstdlib> //libabry for abs()

#include "tools/integration.hpp" //Header file contatining integration routines
#include "tools/nr3.h"
#include "tools/fermi_integrals.hpp"
//#include "tools/NewtonRaphson.hpp"
//struct FD_params { int k; double eta; double alpha; };


int main (){
	int n = 32;
	double a = 0., b = 1.;
	double beta = 0.; //alpha = 5.;
	double result;
	double GSL_lag[2], NR_lag[2], NR_leg[3];

	int k = 5.;
	double t;

	int max_iter = 100;
	double x, x0;

	//Test integration for f1 reconstruction in optically thick limit
	printf("Testing integration for f1 in optically thick limit:\n");
        my_function F_leg;
	F_leg.function = &FD_laguerre;
        struct FD_params pp1 = {k,100.,5.};
	F_leg.params = &pp1;
	double test_1 = glag_integ(24,F_leg,5.);

	F_leg.function = &FD_legendre;
	F_leg.params = &pp1;
	double test_2 = use_gaulag(24,F_leg,5.);

	printf("1: %.5e, 2: %.5e\n", test_1, test_2);

	return 0;
}

