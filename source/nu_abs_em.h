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

	double c1 = (GF*GF/pi) / pow(h*c/(2*pi),4);
	double c2 = gV*gV+3*gA*gA;

	double eta_np(double np, double nn, double mu_hat, double T)
	{

		if (nn == 0.){
			return 0; //enforce zero rates if no neutrons available
		}

		if (fabs(mu_hat) < 1.e-2) {
			return nn;
		}

		double etanp = (np-nn)/(exp(-mu_hat/T)-1);

		if (etanp < 0.){
			etanp = nn;
		}

			return etanp;
	  }

      	double eta_pn(double np, double nn, double mu_hat, double T)
	{

		if (np == 0.){
			return 0; //enforce zero rates if no protons available
		}

		if (fabs(mu_hat) < 1.e-2) {
			return np;
		}

		double etapn = (nn-np)/(exp(mu_hat/T)-1);

		if (etapn < 0.) {
			etapn = np;
		}

		return etapn;
	}

	double theta(double x)
        {
                if (x < 0.)
                        return 0.;
                else
                        return 1.;
        }

	
	//electron neutrino absorption on neutrons
	std::tuple<double,double> nu_n_abs(double omega, double nb, double T, double lep_mass, double yp, double yn, double mu_l, double mu_hat, double deltaU)
	{
		
		double nn = nb*yn; 
		double np = nb*yp; 
		double mu_np;
		double dU, E, R;

		if (use_dU == 1) {
			dU = deltaU;
		} else {
			dU = 0;
		}
		
		E = omega + Q + dU; 
		mu_np = mu_hat - dU;

		if (use_WM_ab == 1) {
			R = WM_nue_abs(omega);
		} else {
			R = 1.;
		}

		double ab = R * c1 * c2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
		                pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(mu_l,E,T));

		//double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/T); //to be checked

		double em = R * c * c1 * c2 * eta_pn(np,nn,mu_np,T) * pow(E,2)*\
                                pow((1-pow(lep_mass/E,2)),0.5) * Fermi(mu_l,E,T);

		return std::make_tuple(em,ab); 
	}

	//electron antineutrino absorption on protons
	std::tuple<double,double> nu_p_abs(double omega, double nb, double T, double lep_mass, double yp, double yn, double mu_l, double mu_hat, double deltaU)
	{
		double nn = nb*yn;
		double np = nb*yp;
		double mu_np;
		double dU, E, Rbar;

                if (use_dU == 1) {
                        dU = deltaU;
                } else {
                        dU = 0;
                }

                E = omega - Q - dU;
                mu_np = mu_hat - dU;

		if (use_WM_ab == 1) {
                        Rbar = WM_nue_abs(omega);
                } else {
                        Rbar = 1.;
		}

		double ab = Rbar * c1 * c2 * eta_pn(np,nn,mu_np,T) * pow(E,2) *\
		                   pow((1-pow(lep_mass/E,2)),0.5) * (1.-Fermi(-mu_l,E,T)) * theta(E-lep_mass);
		
		//double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/T); //to be checked

		double em = Rbar * c * c1 * c2 * eta_np(np,nn,mu_np,T) * pow(E,2) *\
                                   pow((1-pow(lep_mass/E,2)),0.5) * Fermi(-mu_l,E,T) * theta(E-lep_mass);
		
		return std::make_tuple(em,ab);
	}

}
