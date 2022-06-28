//Compute scattering kernel for neutrino elastic scattering on nucleons

#pragma once //compile only once

#include <cmath> //Needed for sqrt computation
#include "constants.h" //Header file containing all relevant constants
#include "parameters.h" //Header file containing all parameters
#include "nucfrmfac.h" //Header file for computation of nucleon form factors

using namespace constants;
using namespace parameters;
using namespace test;

namespace elastic_scatt
{
	// neutral nucleon current form factors
	const double hnv = -0.5;
	const double hna = -0.5*gA;
	const double hpv =  0.5-2.*sinsqthetaw;
	const double hpa =  0.5*gA;
	const double hbar = 0.5*h/pi; //[MeV*s]	
	
	double c0 = (4.*pi*GF*GF)/h;
	double c1 = (4./3.*pi*GF*GF)/h;
	
	std::tuple<double,double> WM_scatt(double Enu, int reacflag){
		
		double R, Rbar;
		
		if (use_WM_sc != 0) {
			double cv, ca, F2;
			std::tie(cv,ca,F2) = nucfrmfac(Enu,reacflag); //form factors
        		double x = 0.; //assume an average angle x=0

          		double ehor = Enu* MeV/(mb*c*c);
			double tmp = (4.*ca*(cv+F2)) / (cv*cv*(1.+x) + ca*ca*(3.-x));
		        R    = (1.+(tmp-3.) *ehor*(1.-x));
          		Rbar = (1.+(-tmp-3.)*ehor*(1.-x));
		} else {
			R = 1.;
			Rbar = 1;
		}
		//printf("%.6e, %.6e\n", R, Rbar);
		return  std::make_tuple(R,Rbar);
	}


	//double eta_nn(double rho, double T, double Ye, double yn){
		//dmudrho = ....; //[MeV*cm^3/g]
		//etann = yn / (T*mu*dmudrho); //[1/cm^3]
		//return etann;
		
		
	double nu_n_scatt(double Enu, double rho, double T, double yn){
		double nn = yn*rho/mb;
		double eFn, tmp, etann;

		double phi0, phi1;
		double B_IS_n;
		
		// degeneracy parameter eta_nn
		if (nn < 0.) {
			etann = 0.;
		} else {
			eFn = 0.5*hbar*hbar * pow(3.*pi*pi*nn,2./3.) / mb * MeV; //from erg to MeV

			tmp = 1.5*T/eFn;
			etann = nn*tmp/sqrt(1.+tmp*tmp);
		}

		phi0 = c0*etann*(hnv*hnv + 3.*hna*hna);
		phi1 = c1*etann*(hnv*hnv -    hna*hna);
		
		B_IS_n = (2.*pi/c)/pow(2.*pi*hbar*c,3.) * Enu*Enu * (phi1-phi0);
		
		return B_IS_n;			
	}


	double nu_p_scatt(double Enu, double rho, double T, double yp){
		double np = yp*rho/mb;
		double eFp, tmp, etapp;

		double phi0, phi1;
		double B_IS_p;
		
		//double E_
		// degeneracy parameter eta_pp
		if (np < 0.) {
			etapp = 0.;
		} else {
			eFp = 0.5*hbar*hbar * pow(3.*pi*pi*np,2./3.) / mb * MeV;
			tmp = 1.5*T/eFp;
			etapp = np*tmp/sqrt(1.+tmp*tmp);
		}

		phi0 = c0*etapp*(hpv*hpv + 3.*hpa*hpa);
		phi1 = c1*etapp*(hpv*hpv -    hpa*hpa);
		
		B_IS_p = (2.*pi/c)/pow(2.*pi*hbar*c,3.) * Enu*Enu * (phi1-phi0);
		
		return B_IS_p;			
	}

	double nu_N_scatt_tot(double Enu, double rho, double T, double yn, double yp){
		double Rp, Rn, Rbarp, Rbarn;
		std::tie(Rp,Rbarp) = WM_scatt(Enu, 1);
		std::tie(Rn,Rbarn) = WM_scatt(Enu, 2);
		double B_IS_tot;

		B_IS_tot = Rn*nu_n_scatt(Enu,rho,T,yn) + Rp*nu_p_scatt(Enu,rho,T,yp);

		return B_IS_tot;
	}
	
	double anu_N_scatt_tot(double Enu, double rho, double T, double yn, double yp){
		double Rp, Rn, Rbarp, Rbarn;
		std::tie(Rp,Rbarp) = WM_scatt(Enu, 1);
		std::tie(Rn,Rbarn) = WM_scatt(Enu, 2);
		double B_IS_tot;

		B_IS_tot = Rbarn*nu_n_scatt(Enu,rho,T,yn) + Rbarp*nu_p_scatt(Enu,rho,T,yp);

		return B_IS_tot;
	}
}
