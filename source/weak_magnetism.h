// Calculation of weak magnetism correction for neutrino emission/absorption and scattering

#pragma once //compile only once

#include "constants.h" //Header file containing all relevant constants
#include "parameters.h" //Header file containing all relevant constants
#include "nucfrmfac.h" //Header file for computation of nucleon form factors

using namespace constants;
using namespace parameters;
using namespace formfactors;

namespace weakmag
{

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


  std::tuple<double,double> WM_scatt(double Enu, int reacflag){

                double R, Rbar;
		double cv, ca, F2;
		std::tie(cv,ca,F2) = nucfrmfac(Enu,reacflag); //form factors
		double x = 0.; //assume an average angle x=0

		double ehor = Enu* MeV/(mb*c*c);
		double tmp = (4.*ca*(cv+F2)) / (cv*cv*(1.+x) + ca*ca*(3.-x));
		R    = (1.+(tmp-3.) *ehor*(1.-x));
		Rbar = (1.+(-tmp-3.)*ehor*(1.-x));

                return  std::make_tuple(R,Rbar);

        }

}
