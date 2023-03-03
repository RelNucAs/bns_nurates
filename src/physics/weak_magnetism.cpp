// !\file weak_magnetism.cpp
// \brief Evaluation of phase space/recoil/weak magnetism correction for (anti)neutrino
//        emission/absorption on nucleons and elastic scattering on nucleons
//        Ref: Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)


#include <tuple>

#include "constants.hpp"              //Header file for physical constants
#include "parameters.hpp"             //Header file for code parameters
#include "physics/nucfrmfac.hpp"      //Header file for nucleon form factors
#include "physics/weak_magnetism.hpp" //Header file for weak magnetism correction

using namespace constants;
using namespace parameters;
using namespace formfactors;

// Correction for electron neutrino absorption on neutron (nue + n -> e- + p): reacflag = 3
// Input: neutrino energy [MeV]
double WM_nue_abs(const double e_nu) { 
	double cv, ca, F2;
	std::tie(cv,ca,F2) = nucfrmfac(e_nu,3); //nuclear form factors

	const double ehor = e_nu * MeV/(mb*c*c); //Eq.(4)
	const double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
	const double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
	//const double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
	const double tmp3 = (gV*gV+3.0*gA*gA)*pow(1.+2.*ehor,3.);

	return (tmp1+tmp2)/tmp3; //Eq.(22)
}


// Correction for electron antineutrino absorption on proton (anue + p -> e+ + n): reacflag = 3
// Input: neutrino energy [MeV]
double WM_anue_abs(const double e_nu) {
	double cv, ca, F2;
	std::tie(cv,ca,F2) = nucfrmfac(e_nu,3); //nuclear form factors

	const double ehor = e_nu * MeV/(mb*c*c); //Eq.(4)
	const double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
	const double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
	//const double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
	const double tmp3 = (gV*gV+3.0*gA*gA)*pow(1.+2.*ehor,3.);
	
	return (tmp1-tmp2)/tmp3; //Eq.(22)
}

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N): reacflag = 1 | 2
// Input: neutrino energy [MeV]
// Output: correction to zeroth (R0) and first Legendre (R1) coefficients of scattering kernel
std::tuple<double,double> WM_scatt(const double Enu, const int reacflag) {
	double cv, ca, F2;
	
	std::tie(cv,ca,F2) = nucfrmfac(Enu,reacflag); //nuclear form factors
	//std::tie(cv_0,ca_0,F2_0) = nucfrmfac(0.,reacflag); //nuclear form factors at Q^2=0
	 
	if (reacflag == 1) {
		const double h0 = hpv*hpv + 3.*hpa*hpa;
		const double h1 = hpv*hpv -    hpa*hpa;
	} else if (reacflag == 2) {
		const double h0 = hnv*hnv + 3.*hna*hna;
		const double h1 = hnv*hnv -    hna*hna;
	}

	const double ehor = Enu* MeV/(mb*c*c);
	
	/* Low-energy limit derived from Eq.(12) */
	const double R0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0; //correction to zeroth coefficient
	const double R1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1; //correction to first coefficient

	return  std::make_tuple(R0,R1);
}
