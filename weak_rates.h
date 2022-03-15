//Opacities

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates

using namespace corrections;
using namespace constants;

namespace weakrates
{
	//electron neutrino absorption on neutrons
	double nu_n_abs(double nb, double temp, double ye, double e_nu, double e_nu_bar, double mu_e, double mu_np,  double ne, double pe, double ee)
	{
		double sigma_nun_abs = nb*sigma_o*((1+3*gA*gA)/4)*((e_nu+delta_np)/e_rm)*pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*Wm(e_nu)*B_nu(e_nu, e_nu_bar, mu_e, mu_np, temp, ne, pe, ee);
		return sigma_nun_abs; 
	}

	//electron neutrino absorption on protonss
	double nu_p_abs(double nb, double temp, double ye, double e_nu, double e_nu_bar, double mu_e, double mu_np,  double ne, double pe, double ee)
	{
		double sigma_nup_abs = nb*sigma_o*((1+3*gA*gA)/4)*((e_nu_bar-delta_np)/e_rm)*pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*Wm_bar(e_nu_bar)*B_nu_bar(e_nu, e_nu_bar, mu_e, mu_np, temp, ne, pe, ee);
		return sigma_nup_abs;
	}
}
