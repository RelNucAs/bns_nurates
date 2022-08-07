//Opacities

#pragma once //compile only once

#include "constants.h" //Header file for physical constants
#include "fermi.h" //Header file for Fermi distribution
#include "parameters.h" // Header file for code parameters
#include "nucfrmfac.h" //Header file for nucleon form factors
#include "weak_magnetism.h" //Header file for weak magnetism corrections

using namespace weakmag;
using namespace constants;
using namespace parameters;
using namespace fermi;

namespace nuabsem
{

	double eta_np(double np, double nn, double mu_np, double T, double dU)
	{
		if (nn == 0.){
			return 0; //enforce zero rates if no neutrons available
		}

		if (fabs(mu_np) < 1.e-2) {
			return nn;
		}

		double etanp = (np-nn)/(exp(-mu_np/T)-1);

		if (etanp < 0.){
			etanp = nn;
		}

			return etanp;
	  }

      	double eta_pn(double np, double nn, double mu_np, double T, double dU)
	{
		if (np == 0.){
			return 0; //enforce zero rates if no protons available
		}

		if (fabs(mu_np) < 1.e-2) {
			return np;
		}

		double etapn = (nn-np)/(exp(mu_np/T)-1);

		if (etapn < 0.) {
			etapn = np;
		}

		return etapn;
	}

	double theta(double e)
        {
                if (e-delta_np-me < 0.)
                        return 0.;
                else
                        return 1.;
        }


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
