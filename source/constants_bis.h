//Unit conversions
#pragma once //compile only once

#include <cmath> //math library for basic operations, pi

namespace constants
{
	// Units conversion
	const double MeV = 1.602176634e-6; //convert MeV to CGS (erg)
	const double cm = 10e13; //convert fm to cm 

 	
	// Electron mass
	const double me = 0.51099895; // MeV
	
	// Muon mass
	const double mmu = 105.65837; // MeV
	
	// Neutron mass
	//const double mn = 939.565; // MeV
	const double mn = 939.5654133; // MeV
	
	// Proton mass 
	//const double mp = 938.272; // MeV
	const double mp = 938.2720813; // MeV
	const double Q = 1.2935; // MeV
	
	// Atomic mass unit
	const double mb = 1.674e-24; // g
        const double mu = 1.66054e-24; //g
	
	// Physical constants
	const double h_bar = 6.582*(10e-22); // MeV*s
	//const double h = 4.1356943e-21; // MeV*s
 	const double c = 2.997924562e+10; // cm/s
	const double h = 1.23984193e-10/c;
 	const double GF = 8.957e-44; // MeV*cm^3
	const double pi = acos((double) -1.);
        const double kB = 8.617333262145e-11; // MeV/K

	// Coupling constants
	const double gA = 1.23;
	const double gV = 1.;
	const double gS = 0.;
	const double sinsqthetaw = 0.2325;

        // Neutral current nucleon form factors (Q^2=0)
        const double hnv = -0.5;
        const double hna = -0.5*gA;
        const double hpv =  0.5-2.*sinsqthetaw;
        const double hpa =  0.5*gA;

}
