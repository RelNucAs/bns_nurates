//corrections
//Fermi distributions

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants

using namespace constants;


namespace corrections
{
  //Fermi distributions

  double Fn(double mu_np, double ne)
  {
  	return 1/(exp(ne-mu_np)+1);
  }

  double Fp(double mu_np, double pe)
  {
  	return 1/(exp(pe-mu_np)+1);
  }

  double Fe(double mu_e, double ee)
  {
  	return 1/(exp(ee-mu_e)+1);
  }

  double Fnu(double mu_np, double mu_e, double e_nu)
  {
  	return 1/(exp(1-mu_np-mu_e-e_nu));
  }

  double F_tilde_nu(double mu_np, double mu_e, double e_nu_bar)
  {
  	return 1/(exp(1-mu_np-mu_e-e_nu_bar));
  }

	double Feqnu(double e_nu, double mu_e, double mu_np, double temp)
	{
		return pow(exp((e_nu-mu_e+mu_np)/(kb*temp))+1,-1);  //beta equilibrium assumed (mu_nu = mu_e-mu_np)
	}

	double Feqnu_bar(double e_nu_bar, double mu_e, double mu_np, double temp)
	{
		return pow(exp((e_nu_bar+mu_e-mu_np)/(kb*temp))+1,-1);  //beta equilibrium assumed (mu_nu = mu_e-mu_np)
	}

	//blocking factor and stimulated absorption
	double B_nu(double e_nu, double e_nu_bar, double mu_e, double mu_np, double temp, double ne, double pe, double ee)
	{
		return (Fn(mu_np, ne)*(1-Fe(mu_e, ee))*(1-Fp(mu_np, pe))*(Feqnu(e_nu, mu_e, mu_np, temp)-Fnu(mu_np, mu_e, e_nu)))/(1-F_tilde_nu(mu_np, mu_e, e_nu_bar));
	} 

	double B_nu_bar(double e_nu, double e_nu_bar, double mu_e, double mu_np, double temp, double ne, double pe, double ee)
	{
		return (Fn(mu_np, ne)*(1-Fe(mu_e, ee))*(1-Fp(mu_np, pe))*(Feqnu_bar(e_nu_bar, mu_e, mu_np, temp)-Fnu(mu_np, mu_e, e_nu)))/(1-F_tilde_nu(mu_np, mu_e, e_nu_bar));
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
}
