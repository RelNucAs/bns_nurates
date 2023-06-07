#include <math.h>

//#include "../../constants.h"
#include "weak_magnetism.h"
//#include "parameters.hpp"             //Header file for code parameters
//#include "physics/nucfrmfac.hpp"      //Header file for nucleon form factors
//#include "physics/weak_magnetism.hpp" //Header file for weak magnetism correction

// !\file weak_magnetism.c
// \brief Evaluation of phase space/recoil/weak magnetism correction for (anti)neutrino
//        emission/absorption on nucleons and elastic scattering on nucleons
//        Ref: Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)


// TODO: change names of variables and functions following Google C style



// TODO: function that computes WM both for abs on n and p
//       (useful only if there exists a single function for the two reactions) 

// Correction for electron neutrino absorption on neutron (nue + n -> e- + p): reacflag = 3
// Input: neutrino energy [kMeV]
double WM_nue_abs(const double e_nu) { 
  double cv, ca, F2;
  nucfrmfac(e_nu, &cv, &ca, &F2, 3); //nuclear form factors
  
  const double ehor = e_nu * kMeV/(kMb*kClight*kClight); //Eq.(4)
  const double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
  const double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
  //const double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
  const double tmp3 = (kGv*kGv+3.0*kGa*kGa)*pow(1.+2.*ehor,3.);
  
  return (tmp1+tmp2)/tmp3; //Eq.(22)
}


// Correction for electron antineutrino absorption on proton (anue + p -> e+ + n): reacflag = 3
// Input: neutrino energy [kMeV]
double WM_anue_abs(const double e_nu) {
  double cv, ca, F2;
  nucfrmfac(e_nu, &cv, &ca, &F2, 3); //nuclear form factors
  
  const double ehor = e_nu * kMeV/(kMb*kClight*kClight); //Eq.(4)
  const double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
  const double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
  //const double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
  const double tmp3 = (kGv*kGv+3.0*kGa*kGa)*pow(1.+2.*ehor,3.);
  
  return (tmp1-tmp2)/tmp3; //Eq.(22)
}

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N): reacflag = 1 | 2
// Input: neutrino energy [kMeV]
// Output: correction to zeroth (R0) and first Legendre (R1) coefficients of scattering kernel
void WM_scatt(const double e_nu, double* R0, double* R1, const int reacflag) {
  double cv, ca, F2;
  double h0, h1;
  double wm_0, wm_1;

  nucfrmfac(e_nu, &cv, &ca, &F2, reacflag); //nuclear form factors
  //nucfrmfac(0., &cv_0, &ca_0, &F2_0, reacflag); //nuclear form factors at Q^2=0
 
  // TODO: evaluate this at compile time  
  if (reacflag == 1) {
    h0 = kHpv*kHpv + 3.*kHpa*kHpa;
    h1 = kHpv*kHpv -    kHpa*kHpa;
  } else if (reacflag == 2) {
    h0 = kHnv*kHnv + 3.*kHna*kHna;
    h1 = kHnv*kHnv -    kHna*kHna;
  }
  
  const double ehor = e_nu * kMeV/(kMb*kClight*kClight);
  
  /* Low-energy limit derived from Eq.(12) */
  wm_0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0; //correction to zeroth coefficient
  wm_1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1; //correction to first coefficient

  R0 = &wm_0;
  R1 = &wm_1;

  return;
}
