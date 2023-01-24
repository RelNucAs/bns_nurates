//Opacities

#include <tuple>

#include "constants.hpp" //Header file for physical constants
#include "fermi.hpp" //Header file for Fermi distribution
#include "parameters.hpp" // Header file for code parameters
#include "nucfrmfac.hpp" //Header file for nucleon form factors
#include "weak_magnetism.hpp" //Header file for weak magnetism corrections
#include "nu_abs_em.hpp"

using namespace constants;
using namespace parameters;
using namespace nuabsem;

double eta_np(const double np, const double nn, const double mu_hat, const double T)
{

	if (nn == 0.){
		return 0; //enforce zero rates if no neutrons available
	}

	if (fabs(mu_hat) < 1.e-2) {
		return nn;
	}

	double etanp = (np-nn)/(exp(-mu_hat/T)-1.);

	if (etanp < 0.){
		etanp = nn;
	}

		return etanp;
  }

double eta_pn(const double np, const double nn, const double mu_hat, const double T)
{

	if (np == 0.){
		return 0; //enforce zero rates if no protons available
	}

	if (fabs(mu_hat) < 1.e-2) {
		return np;
	}

	double etapn = (nn-np)/(exp(mu_hat/T)-1.);

	if (etapn < 0.) {
		etapn = np;
	}

	return etapn;
}

double theta(const double x)
{
	if (x < 0.)
		return 0.;
	else
		return 1.;
}


//electron neutrino absorption on neutrons
std::tuple<double,double> nu_n_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU)
{
	
	double nn = nb*yn; 
	double np = nb*yp; 
	double mu_np;
	double dU, Qprime, E, R;
	double em, ab;

	if (use_dU == 1) {
		dU = deltaU;
	} else {
		dU = 0.;
	}
	
	Qprime = Q + dU;
	E = omega + Qprime; 
	mu_np = mu_hat - dU;

	if (use_WM_ab == 1) {
		R = WM_nue_abs(omega);
	} else {
		R = 1.;
	}

	if (E-lep_mass > 0.) {
		ab = R * g1 * g2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
		     pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(mu_l,E,T)); //*theta(E-lep_mass);

		em = R * c * g1 * g2 * eta_pn(np,nn,mu_np,T) * pow(E,2)*\
		     pow((1-pow(lep_mass/E,2)),0.5) * Fermi(mu_l,E,T); //*theta(E-lep_mass);
	} else {
		ab = 0;
		em = 0;
	}

	//double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/T); //to be checked

	return std::make_tuple(em,em+ab); //(em,ab)
}

//electron antineutrino absorption on protons
std::tuple<double,double> nu_p_abs(const double omega, const double nb, const double T, const double lep_mass, const double yp, const double yn, const double mu_l, const double mu_hat, const double deltaU)
{
	double nn = nb*yn;
	double np = nb*yp;
	double mu_np;
	double dU, Qprime, E, Rbar;
	double ab, em;

	if (use_dU == 1) {
		dU = deltaU;
	} else {
		dU = 0;
	}

	Qprime = Q + dU;
	E = omega - Qprime;
	mu_np = mu_hat - dU;

	if (use_WM_ab == 1) {
		Rbar = WM_anue_abs(omega);
	} else {
		Rbar = 1.;
	}
	
	if (E-lep_mass > 0.) {
		ab = Rbar * g1 * g2 * eta_pn(np,nn,mu_np,T) * pow(E,2) *\
		     pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(-mu_l,E,T)); // * theta(E-lep_mass);
		
		em = Rbar * c * g1 * g2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
		    pow((1-pow(lep_mass/E,2)),0.5) * Fermi(-mu_l,E,T); // * theta(E-lep_mass);
	} else {
		ab = 0.;
		em = 0.;
	}

	//double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/T); //to be checked

	return std::make_tuple(em,em+ab); //(em,ab)
}
