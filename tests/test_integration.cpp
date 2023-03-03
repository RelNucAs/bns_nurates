// \file  test_integration.cpp
// \brief Test integration using Gaussian Quadrature rule (Gauss-Laguerre
//        or Gauss-Legendre split into two pieces)

#include <iostream>
#include <fstream>
#include <string> 
#include <vector>
#include <sstream>
#include <cmath> 
#include <cstdio>  
#include <cstdlib>

#include "tools/integration.hpp" //Header file for integration routines

int main (){
	const int n = 32;
	
	const double a = 0., b = 1.;
	const double alpha = 0.;

	const int k = 5;
        const double eta = 5.;	
	
	my_function FD;

	save_gauleg(n, a, b);
	save_gaulag(n, 0.);

	
	// Gauss-Laguerre
	struct FD_params params = {k,eta,alpha};

	FD.function = &FD_laguerre;
	FD.params = &params;

	cout << "Gauss-Laguerre integration: " << glag_integ(n, FD, alpha) << endl;

	// Gauss-Legendre
	FD.function = &FD_legendre;
	FD.params = &params;

	//const double x0 = find_guess(1000, FD);
	
	cout << "Gauss-Legendre integration: " << gleg_integ(n, FD, eta)   << endl;

	return 0;
}

