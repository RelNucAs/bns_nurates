//corrections
//Fermi distributions

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "nucfrmfac.h" //Header file for computation of nucleon form factors

using namespace constants;
using namespace test;


namespace corrections
{
  //number densities

  double eta_np(double np, double nn, double mu_np, double T)
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

  double eta_pn(double np, double nn, double mu_np, double T)
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


  //Fermi distributions

  double Fermi(double mu, double E, double T)
	{
  		return 1/(exp((E-mu)/T)+1);
	}


  //blocking factor and stimulated absorption 

  double blocking_factor_nu(double e_nu, double mu_e, double T)
	{
		return 1 - Fermi(mu_e, e_nu+delta_np,T);
	} 

  double blocking_factor_nu_bar(double e_nu_bar, double mu_e, double T)
	{
		return 1 - Fermi(mu_e, e_nu_bar-delta_np,T);
	}

  //weak magnetism

  double WM_nue_abs(double e_nu)
	{
		double cv, ca, F2;
		std::tie(cv,ca,F2) = nucfrmfac(e_nu,3); //form factors

		double ehor = e_nu * MeV/(mb*c*c);
		double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
		double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
		double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
                return (tmp1+tmp2)/tmp3;
        }


  double WM_anue_abs(double e_nu)
        {
                double cv, ca, F2;
                std::tie(cv,ca,F2) = nucfrmfac(e_nu,3); //form factors

                double ehor = e_nu * MeV/(mb*c*c);
                double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
                double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
                double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);

                return (tmp1-tmp2)/tmp3;
        }


  //Miscillenous

  double theta(double e)
	{
		if (e-delta_np-me < 0)
			return 0;
		else
			return 1;
	}

}
