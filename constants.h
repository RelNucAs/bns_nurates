//Unit conversions
#pragma once //compile only once

#include <cmath> //math library for basic operations, pi

namespace constants
{
	const double eV = 1.6*10e-12; //convert eV to CGS
	const double cm = 10e13; //convert fm to cm 

	//Global Constants (MeV)

	const double h_bar = 6.582*(10e-22);
	const double c = 29979245800;
	const double Gf = 1.436*(10e-49);
	const double me = 0.5109;
	const double mn = 939.565;
	const double mp = 938.272;
	const double pi = M_PI;
	const double kb = 1; //in eV/eV 

	//Simplying constants
	const double e_rm = me;
	const double n_rm = mn;
	const double p_rm = mp;
	const double Gs = 5.88*(10e-18); //in Mev-2fm2 //(4*(Gf*Gf)*(e_rm*e_rm))/(pi); //(pi*pow(h_bar*c,4))
	const double gA = 1.23;
	const double gV = 1;
	const double delta_np = 1.29332;
}
