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
	std::tuple<double,double> nu_n_abs(double rho, double T, double Ml, double yp, double yn, double omega, double mu_l, double mu_hat)
	{
		double nn = rho*yn/mb; 
		double np = rho*yp/mb; 
		double R;

		if (use_WM_ab == 1) {
			R = WM_nue_abs(omega);
		} else {
			R = 1.;
		}

		double ab = R * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_np(np,nn,mu_hat,T)*\
		                 (gV*gV+3*gA*gA)*(pow(omega+delta_np,2))*\
		                pow((1-pow((Ml/(omega+delta_np)),2)),0.5)*\
                                (1.-Fermi(mu_l,omega+delta_np,T));

		//double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/T);
		double em = R * c * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_pn(np,nn,mu_hat,T)*\
                                 (gV*gV+3*gA*gA)*(pow(omega+delta_np,2))*\
                                pow((1-pow((Ml/(omega+delta_np)),2)),0.5)*\
                                Fermi(mu_l,omega+delta_np,T);

		return std::make_tuple(em,ab); 
	}

	//electron antineutrino absorption on protons
	std::tuple<double,double> nu_p_abs(double rho, double T, double Ml, double yp, double yn, double omega, double mu_l, double mu_hat)
	{
		double nn = rho*yn/mb;
		double np = rho*yp/mb;
		double Rbar;

		if (use_WM_ab == 1) {
                        Rbar = WM_nue_abs(omega);
                } else {
                        Rbar = 1.;
		}

		double ab = Rbar * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_pn(np,nn,mu_hat,T)*\
		                    (gV*gV+3*gA*gA)*(pow(omega-delta_np,2))*\
		                   pow((1-pow((Ml/(omega-delta_np)),2)),0.5)*\
		                   (1.-Fermi(-mu_l, omega-delta_np,T))*theta(omega);
		
		//double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/T);
		double em = Rbar * c * ((GF*GF/pi)/pow(h*c/(2*pi),4))*eta_np(np,nn,mu_hat,T)*\
                                    (gV*gV+3*gA*gA)*(pow(omega-delta_np,2))*\
                                   pow((1-pow((Ml/(omega-delta_np)),2)),0.5)*\
                                   Fermi(-mu_l, omega-delta_np,T)*theta(omega);
		return std::make_tuple(em,ab);
	}
}
