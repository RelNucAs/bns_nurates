#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#include "../../constants.h"
#include "weak_magnetism.h"

//! \file nucfrmfac.c
//  \brief Calculation of single nucleon form factors as in C.J. Horowitz, 2002
//         (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001).
//         These are needed to compute the recoil and weak magnetism correction
//         for (anti)neutrino absorption on nucleons and elastic scattering on
//         nucleons.

// Reactions are distinguished using the following indices:
//   reacflag = 0: (anti)neutrino scattering on proton  (nu p -> nu p)
//   reacflag = 1: (anti)neutrino scattering on neutron (nu n -> nu n)
//   reacflag = 2: (anti)neutrino absorption on nucleon (nue n -> e- p, anue p -> e+ n)

// @TODO: change names of variables and functions following Google C style


// Nucleon constants
const double lamp =  1.793; //  proton magnetic moment?
const double lamn = -1.913; // neutron magnetic moment

// Computation of single nucleon form factors for reaction reacflag,
// given the (anti)neutrino energy

/*
 * Input:
 * 	- E : (anti)neutrino energy [MeV]
 * 	- reacflag : index defining the reaction (see above)
 * 
 * Output:
 * 	- cv : vector form factor
 * 	- ca : axial vector form factor
 * 	- F2 : tensor/Pauli form factor
 *
*/

void nucfrmfac(const double E, double* cv, double* ca, double* F2, const int reacflag) {
  /* (Anti)neutrino energy rescaled by the nucleon mass */
  const double ehor = E * kMeV/(kMb*kClight*kClight); //Eq.(4), dimensionless

  const double tau = 0.5*ehor*ehor/(1.+ehor);       //Eq.(B10) 
  const double eta = 1./(1.+5.6*tau);               //Eq.(B16)
  const double G   = 1./pow(1.+4.97*tau,2.);        //Eq.(B17)
  const double Fp1 = (1.+tau*(1.+lamp))*G/(1.+tau); //Eq.(B11) 
  const double Fp2 = lamp*G/(1.+tau);               //Eq.(B12)
  const double Fn1 = tau*lamn*(1.-eta)*G/(1.+tau);  //Eq.(B13)
  const double Fn2 = lamn*(1.+tau*eta)*G/(1.+tau);  //Eq.(B14)
  
  double frm1, frm2, frm3;

  /* Different parametrization depending on the reaction */
  if (reacflag == 1) {
    frm1  = (0.5-2.*kSinsqthetaw)*Fp1 - 0.5*Fn1; //Eq.(B1)
    frm2  = 0.5*(kGa-kGs)/pow(1.+3.53*tau,2.);    //Eq.(B2)
    frm3  = (0.5-2.*kSinsqthetaw)*Fp2 - 0.5*Fn2; //Eq.(B3)
  } else if (reacflag == 2) {
    frm1 = (0.5-2.*kSinsqthetaw)*Fn1 - 0.5*Fp1;  //Eq.(B4)
    frm2 = -0.5*(kGa+kGs)/pow(1.+3.53*tau,2.);    //Eq.(B5)
    frm3 = (0.5-2.*kSinsqthetaw)*Fn2 - 0.5*Fp2;  //Eq.(B6)
  } else if (reacflag == 3) {
    frm1 = Fp1 - Fn1;                           //Eq.(B7)
    frm2 = kGa/pow(1.+3.53*tau,2.);              //Eq.(B8)
    frm3 = Fp2 - Fn2;                           //Eq.(B9)
  } else {
    printf("Error: reacflag out of range in nucfrmfac\n");
    exit(EXIT_FAILURE);
  }

  cv = &frm1;
  ca = &frm2;
  F2 = &frm3;

  return;
}	
