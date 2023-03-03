// Calculation of single nucleon form factors

#include <cmath>
#include <cstdlib>
#include <tuple>
#include <iostream>
#include <ostream>

#include "constants.hpp"          //Header file for physical constants
#include "physics/nucfrmfac.hpp"  //Header file for nuclear form factors

//! \file nucfrmfac.hpp
//  \brief Calculation of single nucleon form factors as in C.J. Horowitz, 2002
//         (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001).
//         These are needed to compute the recoil and weak magnetism correction
//         for (anti)neutrino absorption on nucleons and elastic scattering on
//         nucleons.

// Reactions are distinguished using the following indices:
//   reacflag = 0: (anti)neutrino scattering on proton  (nu p -> nu p)
//   reacflag = 1: (anti)neutrino scattering on neutron (nu n -> nu n)
//   reacflag = 2: (anti)neutrino absorption on nucleon (nue n -> e- p, anue p -> e+ n)

using namespace constants;
using namespace formfactors;

// Computation of single nucleon form factors for reaction reacflag,
// given the (anti)neutrino energy
std::tuple<double,double,double> nucfrmfac(const double E, const int reacflag) {
	/* (Anti)neutrino energy rescaled by the nucleon mass */
	const double ehor = E * MeV/(mb*c*c); //Eq.(4), dimensionless

	const double tau = 0.5*ehor*ehor/(1.+ehor);       //Eq.(B10) 
	const double eta = 1./(1.+5.6*tau);               //Eq.(B16)
	const double G   = 1./pow(1.+4.97*tau,2.);        //Eq.(B17)
	const double Fp1 = (1.+tau*(1.+lamp))*G/(1.+tau); //Eq.(B11) 
	const double Fp2 = lamp*G/(1.+tau);               //Eq.(B12)
	const double Fn1 = tau*lamn*(1.-eta)*G/(1.+tau);  //Eq.(B13)
	const double Fn2 = lamn*(1.+tau*eta)*G/(1.+tau);  //Eq.(B14)

	double cv, ca, F2;

	/* Different parametrization depending on the reaction */
	if (reacflag == 1) {
		cv  = (0.5-2.*sinsqthetaw)*Fp1 - 0.5*Fn1; //Eq.(B1)
		ca  = 0.5*(gA-gS)/pow(1.+3.53*tau,2.);    //Eq.(B2)
		F2  = (0.5-2.*sinsqthetaw)*Fp2 - 0.5*Fn2; //Eq.(B3)
	} else if (reacflag == 2) {
		cv = (0.5-2.*sinsqthetaw)*Fn1 - 0.5*Fp1;  //Eq.(B4)
		ca = -0.5*(gA+gS)/pow(1.+3.53*tau,2.);    //Eq.(B5)
		F2 = (0.5-2.*sinsqthetaw)*Fn2 - 0.5*Fp2;  //Eq.(B6)
	} else if (reacflag == 3) {
		cv = Fp1 - Fn1;                           //Eq.(B7)
		ca = gA/pow(1.+3.53*tau,2.);              //Eq.(B8)
		F2 = Fp2 - Fn2;                           //Eq.(B9)
	} else {
		std::cout << "Error: reacflag out of range in nucfrmfac" << std::end;
		std::exit(EXIT_FAILURE);
	}

	return std::make_tuple(cv,ca,F2);
}	
