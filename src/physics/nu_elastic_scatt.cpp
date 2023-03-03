// \file nu_elastic_scatt.cpp
// \brief Computation of first two Legendre coefficients for neutrino elastic
//        scattering on nucleons
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)

#include <cmath> 
#include <tuple>

#include "constants.hpp"                 // Header file for physical constants
#include "parameters.hpp"                // Header file for code parameters
#include "physics/weak_magnetism.hpp"    // Header file for weak magnetism correction
#include "physics/nu_elastic_scatt.hpp"  // Header file for elastic scattering kernel

using namespace constants;
using namespace parameters;
using namespace elastic_scatt;

/* Inputs:
 *      omega    [MeV] : (anti)neutrino energy
 *      nb      [cm-3] : baryon number density
 *      T        [MeV] : temperature
 *      yp             : proton fraction
 *      yn             : neutron fraction
*/


/* Computation of degeneracy parameter eta_NN */
double eta_NN(const double nb, const double T, const double yN) {
	const double nN = yN*nb;

	if (nN <= 0.) {
		/* Enforce zero rate if no nucleons are present */
		return 0.;
	} else {
		/* Linear interpolation between degenerate and nondegnerate limit in Eq.(C37) */
		const double eFN = 0.5*hbar*hbar * pow(3.*pi*pi*nN,2./3.) / mb * MeV; // [MeV]
		const double tmp = 1.5*T/eFN;
		return nN*tmp/sqrt(1.+tmp*tmp); // [cm-3]
	}
	
	/* Alternative computation: evaluate numerical derivative (even if some quick tests showed that could be unstable) */
	//dmudrho = ....; //[MeV*cm^3/g]
	//etann = yn / (T*mu*dmudrho); //[1/cm^3]
}

/* Legendre expansion (up to l=1) of scattering kernel */
double nu_N_scatt(const double omega, const double nb, const double T, const double yN, const int reacflag){
	double R0, R1;
	double h0, h1;
	
	const double etaNN = eta_NN(nb, T, yN); //degeneracy parameter eta_NN, Eq.(C37)

	/* Phase space, recoil and weak magnetism corrections */
	if (use_WM_sc != 0) {
		// R0 (R1) is the correction to the zeroth (first) Legendre coefficient
		std::tie(R0,R1) = WM_scatt(omega, reacflag);
	} else {
		R0 = 1.;
		R1 = 1.;
	}

	if (reacflag == 1) {
		// Scattering on proton
		h0 = hpv*hpv + 3.*hpa*hpa;
		h1 = hpv*hpv -    hpa*hpa;			
	} else if (reacflag == 2) {
		// Scattering on neutron
		h0 = hnv*hnv + 3.*hna*hna;
		h1 = hnv*hnv -    hna*hna;			
	}

	const double phi0 = R0*c0*etaNN*h0; // zeroth Legendre coefficient, Eq.(C38)
	const double phi1 = R1*c1*etaNN*h1; // first  Legendre coefficient, Eq.(C39)
	
	return  (2.*pi/c)/pow(2.*pi*hbar*c,3.) * omega*omega * (phi1-phi0); // scattering kernel, Eq.(A40)
	
}

/* Sum proton and neutron scattering kernels */
double nu_N_scatt_tot(const double omega, const double nb, const double T, const double yn, const double yp){
	return nu_N_scatt(omega,nb,T,yn,2) + nu_N_scatt(omega,nb,T,yp,1); //neutron + proton contribution
}
