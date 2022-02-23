//corrections
//Fermi distributions
#include "constants.h" //Header file containing all relevant constants

using namespace constants;


namespace corrections
{

	double Feqnu(double e_nu, double mu_e, double mu_np, double temp)
	{
		return pow(exp((e_nu-mu_e+mu_np)/(kb*temp))+1,-1)  //beta equilibrium assumed (mu_nu = mu_e-mu_np)
	}

	double Feqnu_bar(double e_nu_bar, double mu_e, double mu_np, double temp)
	{
		return pow(exp((e_nu_bar+mu_e-mu_np)/(kb*temp))+1,-1)  //beta equilibrium assumed (mu_nu = mu_e-mu_np)
	}

	//blocking factor and stimulated absorption
	double B_nu(double e_nu, double mu_e, double mu_np, double temp)
	{
		Bf = (Fn*(1-Fe)*(1-Fp)*(Feqnu(e_nu, mu_e, mu_np, temp)-Fnu))/(1-F_tilde_nu)
		return Bf
	} 

	double B_nu_bar(double e_nu_bar, double mu_e, double mu_np, double temp)
	{
		Bf = (Fn*(1-Fe)*(1-Fp)*(Feqnu_bar(e_nu_bar, mu_e, mu_np, temp)-Fnu))/(1-F_tilde_nu)
		return Bf
	}

	//weak magnetism
	double Wm(double e_nu)
	{
		return (1+1.1*e_nu)/n_rm
	}

	double Wm_bar(double e_nu_bar)
	{
		return (1-7.1*e_nu_bar)/n_rm
	}
}