//Opacities

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates

using namespace corrections;
using namespace constants;

namespace weakrates
{
	//electron neutrino absorption on neutrons
	double nu_n_abs(double nb, double temp, double ye, double e_nu, double mu_e, double mu_np)
	{
		double sigma_nun_abs = (Gs/pi)*eta_np(nb*ye,nb*(1-ye),mu_np,temp)*\
		                       ((gV*gV+3*gA*gA)/4)*(pow((e_nu+delta_np),2))*\
		                       pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*\
		                       B_nu(e_nu, mu_e, temp); //*Wm(e_nu)
		return sigma_nun_abs; 
	}

	//electron neutrino absorption on protonss
	double nu_p_abs(double nb, double temp, double ye, double e_nu_bar, double mu_e, double mu_np)
	{
		double sigma_nup_abs = (Gs/pi)*eta_np(nb*(1-ye),nb*ye,mu_np,temp)*\
		                       ((gV*gV+3*gA*gA)/4)*(pow((e_nu_bar-delta_np),2))*\
		                       pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*\
		                       B_nu_bar(e_nu_bar, mu_e, temp)*theta(e_nu_bar); //*Wm_bar(e_nu_bar)
		return sigma_nup_abs;
	}
}
