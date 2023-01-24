//Compute scattering kernel for neutrino elastic scattering on nucleons

#include <cmath> //Needed for sqrt computation
#include <tuple>

#include "constants.hpp" //Header file containing all relevant constants
#include "parameters.hpp" //Header file containing all parameters
#include "weak_magnetism.hpp" //Header file for implementation of weak magnetism corrections
#include "nu_elastic_scatt.hpp"

using namespace constants;
using namespace parameters;
using namespace elastic_scatt;

//double eta_nn(double rho, double T, double Ye, double yn){
	//dmudrho = ....; //[MeV*cm^3/g]
	//etann = yn / (T*mu*dmudrho); //[1/cm^3]
	//return etann;

// degeneracy parameter eta_NN
double eta_NN(const double nb, const double T, const double yN) {
	double nN = yN*nb;
	double eFN, tmp, etaNN;

	if (nN <= 0.) {
		etaNN = 0.;
	} else {
		eFN = 0.5*hbar*hbar * pow(3.*pi*pi*nN,2./3.) / mb * MeV;
		tmp = 1.5*T/eFN;
		etaNN = nN*tmp/sqrt(1.+tmp*tmp);
	}

	return etaNN;
}

double nu_N_scatt(const double Enu, const double nb, const double T, const double yN, const int reacflag){
	double R0, R1;
	double h0, h1;
	double phi0, phi1;
	double B_IS_n;
	
	double etaNN = eta_NN(nb, T, yN); //degeneracy parameter eta_NN

	if (use_WM_sc != 0) {
		std::tie(R0,R1) = WM_scatt(Enu, reacflag);
	} else {
		R0 = 1.;
		R1 = 1.;
	}

	if (reacflag == 1) {
		h0 = hpv*hpv + 3.*hpa*hpa;
		h1 = hpv*hpv -    hpa*hpa;			
	} else if (reacflag == 2) {
		h0 = hnv*hnv + 3.*hna*hna;
		h1 = hnv*hnv -    hna*hna;			
	}

	phi0 = R0*c0*etaNN*h0;
	phi1 = R1*c1*etaNN*h1;
	
	B_IS_n = (2.*pi/c)/pow(2.*pi*hbar*c,3.) * Enu*Enu * (phi1-phi0);
	
	return B_IS_n;			
}


double nu_N_scatt_tot(const double Enu, const double nb, const double T, const double yn, const double yp){
	return nu_N_scatt(Enu,nb,T,yn,2) + nu_N_scatt(Enu,nb,T,yp,1); //neutron + proton contribution
}
