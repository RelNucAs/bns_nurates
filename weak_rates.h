//Opacities

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates

using namespace corrections;
using namespace constants;

namespace weakrates
{
	//electron neutrino absorption on neutrons
	double nu_n_abs(double nn, double np, double temp, double ye, double e_nu, double mu_e, double mu_np)
	{
		double sigma_nun_abs = ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_np(np,nn,mu_np,temp)*\
		                       (gV*gV+3*gA*gA)*(pow(e_nu+delta_np,2))*\
		                       pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*\
		                       blocking_factor_nu(e_nu, mu_e, temp); //*weak_magnetism(e_nu)
		return sigma_nun_abs; 
	}

	//electron neutrino absorption on protonss
	double nu_p_abs(double nn, double np, double temp, double ye, double e_nu_bar, double mu_e, double mu_np)
	{
		double sigma_nup_abs = ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_pn(np,nn,mu_np,temp)*\
		                       (gV*gV+3*gA*gA)*(pow(e_nu_bar-delta_np,2))*\
		                       pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*\
		                       blocking_factor_nu_bar(e_nu_bar, -mu_e, temp)*theta(e_nu_bar); //*weak_magnetism_bar(e_nu_bar)
		return sigma_nup_abs;
	}
}
