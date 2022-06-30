//Opacities

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "corrections.h" //Header file containing all relevant corrections to rates
#include "parameters.h" // Header file containing WM parameter

using namespace corrections;
using namespace constants;
using namespace parameters;

namespace weakrates
{
	//electron neutrino absorption on neutrons
	std::tuple<double,double> nu_n_abs(double rho, double temp, double ye, double yp, double yn, double e_nu, double mu_e, double mu_np)
	{
		double nn = rho*yn/mb; 
		double np = rho*yp/mb; 
		double R;

		if (use_WM_ab == 1) {
			R = WM_nue_abs(e_nu);
		} else {
			R = 1.;
		}

		double ab = R * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_np(np,nn,mu_np,temp)*\
		                 (gV*gV+3*gA*gA)*(pow(e_nu+delta_np,2))*\
		                pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*\
                                (1.-Fermi(mu_e,e_nu+delta_np,temp));

		//double em = ab * c * exp(-(e_nu-(mu_e-mu_np-delta_np))/temp);
		double em = R * c * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_pn(np,nn,mu_np,temp)*\
                                 (gV*gV+3*gA*gA)*(pow(e_nu+delta_np,2))*\
                                pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*\
                                Fermi(mu_e,e_nu+delta_np,temp);

		//printf("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",rho,temp,ye,yp,yn,B_nu(e_nu, mu_e, temp),eta_np(np,nn,mu_np,temp));
		return std::make_tuple(em,ab); 
	}

	//electron antineutrino absorption on protons
	std::tuple<double,double> nu_p_abs(double rho, double temp, double ye, double yp, double yn, double e_nu_bar, double mu_e, double mu_np)
	{
		double nn = rho*yn/mb;
		double np = rho*yp/mb;
		double Rbar;

		if (use_WM_ab == 1) {
                        Rbar = WM_nue_abs(e_nu_bar);
                } else {
                        Rbar = 1.;
		}

		double ab = Rbar * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_pn(np,nn,mu_np,temp)*\
		                    (gV*gV+3*gA*gA)*(pow(e_nu_bar-delta_np,2))*\
		                   pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*\
		                   (1.-Fermi(-mu_e, e_nu_bar-delta_np,temp))*theta(e_nu_bar);
		
		//double em = ab * c * exp(-(e_nu_bar-(mu_np+delta_np-mu_e))/temp);
		double em = Rbar * c * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_np(np,nn,mu_np,temp)*\
                                    (gV*gV+3*gA*gA)*(pow(e_nu_bar-delta_np,2))*\
                                   pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*\
                                   Fermi(-mu_e, e_nu_bar-delta_np,temp)*theta(e_nu_bar);
		return std::make_tuple(em,ab);
	}
}
