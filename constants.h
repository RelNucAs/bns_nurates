//Unit conversions
#pragma once //compile only once

#include <cmath> //math library for basic operations, pi

namespace constants
{
	const double MeV = 1.602176634e-6; //convert MeV to CGS (erg)
	const double cm = 10e13; //convert fm to cm 

	//Global Constants (MeV)

	const double h_bar = 6.582*(10e-22);
	const double h = 4.1356943e-21; //[MeV*s]
 	const double c = 2.997924562e+10;
 	const double GF = 8.957e-44; //[MeV*cm^3]
 	const double me = 0.510998928; // MeV
	const double mn = 939.565; // MeV
	const double mp = 938.272; // MeV
	const double mb = 1.674e-24; // g
	const double pi = 3.1415926535898;
	const double kb = 1; //in eV/eV 

	//Simplying constants
	const double e_rm = me;
	const double n_rm = mn;
	const double p_rm = mp;
	const double gA = 1.23;
	const double gV = 1.;
	const double gS = 0.;
	const double delta_np = 1.2935; // MeV
	const double sinsqthetaw = 0.2325;
}
