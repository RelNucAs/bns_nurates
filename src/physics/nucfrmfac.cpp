// Calculation of single nucleon form factors

#include <cmath>
#include <tuple>

#include "constants.hpp" //Header file containing all relevant constants
#include "nucfrmfac.hpp"

using namespace constants;
using namespace formfactors;

std::tuple<double,double,double> nucfrmfac(const double E, const int reacflag) {
	const double ehor = E * MeV/(mb*c*c);

	const double tau = 0.5*ehor*ehor/(1.+ehor);
	const double eta = 1./(1.+5.6*tau);
	const double G = 1./pow(1.+4.97*tau,2);
	const double Fp1 = (1.+tau*(1.+lamp))*G/(1.+tau);
	const double Fp2 = lamp*G/(1.+tau);
	const double Fn1 = tau*lamn*(1.-eta)*G/(1.+tau);
	const double Fn2 = lamn*(1.+tau*eta)*G/(1.+tau);

	double cv,ca,F2;

	//form factors
	if (reacflag == 1) {
		cv  = (0.5-2.*sinsqthetaw)*Fp1 - 0.5*Fn1;
		ca  = 0.5*(gA-gS)/pow(1.+3.53*tau,2);
		F2  = (0.5-2.*sinsqthetaw)*Fp2 - 0.5*Fn2;
	} else if (reacflag == 2) {
		cv = (0.5-2.*sinsqthetaw)*Fn1 - 0.5*Fp1;
		ca = -0.5*(gA+gS)/pow(1.+3.53*tau,2);
		F2 = (0.5-2.*sinsqthetaw)*Fn2 - 0.5*Fp2;
	} else if (reacflag == 3) {
		cv = Fp1 - Fn1;
		ca = gA/pow(1.+3.53*tau,2);
		F2 = Fp2 - Fn2;
	} else {
		printf ("Error: reacflag out of range");
		exit (EXIT_FAILURE);
	}
return std::make_tuple(cv,ca,F2);
}	
