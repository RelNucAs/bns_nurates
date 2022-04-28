//corrections
//Fermi distributions

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants

using namespace constants;


namespace corrections
{
  //number densities

  double eta_np(double np, double nn, double mu_np, double temp)
  {
  	return (np-nn)/(exp(-mu_np/(kb*temp))-1);
  }

  double eta_pn(double np, double nn, double mu_np, double temp)
  {
  	return (nn-np)/(exp(mu_np/(kb*temp))-1);
  }

  //Fermi distributions

  double Fermi(double mu, double E, double temp)
  {
  	return 1/(exp((E-mu)/(kb*temp))+1);
  }


	//blocking factor and stimulated absorption 

	double B_nu(double e_nu, double mu_e, double temp)
	{
		return 1 - Fermi(mu_e, e_nu+delta_np,temp);
	} 

	double B_nu_bar(double e_nu_bar, double mu_e, double temp)
	{
		return 1 - Fermi(mu_e, e_nu_bar-delta_np,temp);
	}

	//weak magnetism

	double Wm(double e_nu)
	{
		return (1+1.1*e_nu)/n_rm;
	}

	double Wm_bar(double e_nu_bar)
	{
		return (1-7.1*e_nu_bar)/n_rm;
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
