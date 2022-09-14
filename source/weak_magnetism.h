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
		//double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
		double tmp3 = (gV*gV+3.0*gA*gA)*pow(1.+2.*ehor,3.);
                return (tmp1+tmp2)/tmp3;
        }


  double WM_anue_abs(double e_nu)
        {
                double cv, ca, F2;
                std::tie(cv,ca,F2) = nucfrmfac(e_nu,3); //form factors

                double ehor = e_nu * MeV/(mb*c*c);
                double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
                double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
                //double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
		double tmp3 = (gV*gV+3.0*gA*gA)*pow(1.+2.*ehor,3.);

                return (tmp1-tmp2)/tmp3;
	}


  std::tuple<double,double> WM_scatt(double Enu, int reacflag){

                double R0, R1;
		double cv, ca, F2;
		double h0, h1;
		
		std::tie(cv,ca,F2) = nucfrmfac(Enu,reacflag); //form factors
		//std::tie(cv_0,ca_0,F2_0) = nucfrmfac(0.,reacflag); //form factors at Q^2=0
		  
    		if (reacflag == 1) {
                        h0 = hpv*hpv + 3.*hpa*hpa;
                        h1 = hpv*hpv -    hpa*hpa;
                } else if (reacflag == 2) {
                        h0 = hnv*hnv + 3.*hna*hna;
                        h1 = hnv*hnv -    hna*hna;
                }

		double ehor = Enu* MeV/(mb*c*c);
		
		R0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0;
		R1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1;

                return  std::make_tuple(R0,R1);

        }

}
