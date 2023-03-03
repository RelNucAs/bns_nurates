// \file nu_abs_em.cpp
// \brief Computation of emissivity and absorptivity for neutrino absorption on neutron
//        (and inverse) and for antineutrino absorption on proton (and inverse)
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
//        
// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001) 
//
// Possible inclusion of nucleon interation correction as in 
// Hempel, 2015 (https://journals.aps.org/prc/abstract/10.1103/PhysRevC.91.055807)

#include <tuple>

#include "constants.hpp"              //Header file for physical constants
#include "parameters.hpp"             //Header file for code parameters
#include "tools/fermi.hpp"            //Header file for Fermi-Dirac distribution function
#include "physics/nucfrmfac.hpp"      //Header file for nucleon form factors
#include "physics/weak_magnetism.hpp" //Header file for weak magnetism corrections
#include "physics/nu_abs_em.hpp"      //Header file for nu_abs_em rates

using namespace constants;
using namespace parameters;
using namespace nuabsem;

/* Inputs:
 * 	omega    [MeV] : (anti)neutrino energy
 * 	nb      [cm-3] : baryon number density
 * 	T        [MeV] : temperature
 * 	lep_mass [MeV] : lepton mass (electron/muon) 
 * 	yp             : proton fraction
 * 	yn             : neutron fraction
 * 	mu_l     [MeV] : lepton chemical potential, rest mass included
 * 	mu_hat   [MeV] : neutron-proton chemical potential difference, rest mass NOT included
 * 	dU       [MeV] : nuclear interaction correction on mu_hat
*/


/* Factors resulting from nucleon phase space integration */
double eta_np(const double np, const double nn, const double mu_hat, const double T) {
	if (nn == 0.) return 0; //enforce zero rates if no neutrons available

	if (fabs(mu_hat) < mu_thres) return nn; //backup if mu_hat too small

	const double etanp = (np-nn) / (exp(-mu_hat/T)-1.); // Eq.(C14), [cm-3]

	if (etanp < 0.) return nn; //backup if etanp is negative
	
	return etanp;
}

double eta_pn(const double np, const double nn, const double mu_hat, const double T) {
	if (np == 0.) return 0; //enforce zero rates if no protons available

	if (fabs(mu_hat) < mu_thres) return np; //backup if mu_hat too small

	const double etapn = (nn-np) / (exp(mu_hat/T)-1.); // Eq.(C14), [cm-3]

	if (etapn < 0.) return np; //backup if etapn is negative

	return etapn;
}

/* Theta step function */
double theta(const double x) {
	if (x < 0.) return 0.;
	else        return 1.;
}


/* Neutrino absorption on neutron (nul + n -> l- + p) */
std::tuple<double,double> nu_n_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU) {
	const double nn = nb*yn; //neutron number density [cm-3]
	const double np = nb*yp; //proton number density  [cm-3]
	double em, ab;

	/* Nucleon interaction correction to chemical potential */
	if (use_dU == 1) {
		const double dU = deltaU; // [MeV]
	} else {
		const double dU = 0.;
	}
	const double Qprime = Q + dU;     // [MeV], Eq.(79) in Hempel
	const double E = omega + Qprime;  // [MeV]
	const double mu_np = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

	/* Phase space, recoil and weak magnetism correction */
	if (use_WM_ab == 1) {
		const double R = WM_nue_abs(omega);
	} else {
		const double R = 1.;
	}

	
	/* Check kinematics constraint for the reaction */
	if (E-lep_mass > 0.) {
		// Absoprtivity [s-1]
		ab = R * c * g1 * g2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
		     pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(mu_l,E,T)); //*theta(E-lep_mass); // Eq.(C13)
		
		// Emissivity [s-1]
		em = R * c * g1 * g2 * eta_pn(np,nn,mu_np,T) * pow(E,2) *\
		     pow((1-pow(lep_mass/E,2)),0.5) * Fermi(mu_l,E,T); //*theta(E-lep_mass);      // Eq.(C15)
	} else {
		ab = 0;
		em = 0;
	}

	/* Emissivity from detailed balance (NOT TESTED) */
	//double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/T); //to be checked

	return std::make_tuple(em,ab);
}

/* Antineutrino absorption on neutron (anul + p -> l+ + n) */
std::tuple<double,double> nu_p_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU) {
	const double nn = nb*yn; //neutron number density [cm-3]
	const double np = nb*yp; //proton number density  [cm-3]
	double ab, em;

	/* Nucleon interaction correction to chemical potential */
	if (use_dU == 1) {
		dU = deltaU; // [MeV]
	} else {
		dU = 0;    
	}
	const double Qprime = Q + dU;     // [MeV], Eq.(79) in Hempel
	const double E = omega - Qprime;  // [MeV]
	const double mu_np = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

	/* Phase space, recoil and weak magnetism correction */
	if (use_WM_ab == 1) {
		const double Rbar = WM_anue_abs(omega);
	} else {
		const double Rbar = 1.;
	}
	
	/* Check kinematics constraint for the reaction */
	if (E-lep_mass > 0.) {
		// Absorptivity [s-1]
		ab = Rbar * c * g1 * g2 * eta_pn(np,nn,mu_np,T) * pow(E,2) *\
		     pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(-mu_l,E,T)); // * theta(E-lep_mass); // Eq.(C19)
		
		// Emissivity [s-1]
		em = Rbar * c * g1 * g2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
		    pow((1-pow(lep_mass/E,2)),0.5) * Fermi(-mu_l,E,T); // * theta(E-lep_mass);       // Eq.(C20)
	} else {
		ab = 0.;
		em = 0.;
	}

	/* Emissivity from detailed balance (NOT TESTED) */
	//double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/T);

	return std::make_tuple(em,ab);
}


/* Stimulated absoption versions */
std::tuple<double,double> nu_n_abs_stim(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU) {
        double em, ab;
        std::tie(em,ab) = nu_n_abs(omega, nb, T, lep_mass, yp, yn, mu_l, mu_hat, deltaU);
	return std::make_tuple(em,em+ab)	
}


std::tuple<double,double> nu_p_abs_stim(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU) {
        double em, ab;
        std::tie(em,ab) = nu_p_abs(omega, nb, T, lep_mass, yp, yn, mu_l, mu_hat, deltaU);
	return std::make_tuple(em,em+ab)	
}
