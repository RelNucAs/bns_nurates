#include <math.h>

#include "../../constants.h"
#include "weak_magnetism.h"

// !\file weak_magnetism.c
// \brief Evaluation of phase space/recoil/weak magnetism correction for (anti)neutrino
//        emission/absorption on nucleons and elastic scattering on nucleons
//        Ref: Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)


// R    -> Correction for electron     neutrino absorption on neutron (nu_l + n -> l- + p)
// Rbar -> Correction for electron antineutrino absorption on proton (anu_l + p -> l+ + n)
// reacflag = 3 (for nuclear form factors)
// Input: omega -> neutrino energy [MeV]
void WMAbsEm(const double omega, double* R, double* Rbar) { 
  double cv, ca, F2;
  double R_nu, R_anu;

  NucFrmFac(omega, &cv, &ca, &F2, 3); //nuclear form factors
  
  const double ehor = omega * kMeV/(kMb*kClight*kClight); //Eq.(4)
  const double tmp1 = cv*cv*(1.+4.*ehor+16./3.*ehor*ehor) + 3.*ca*ca*pow(1.+4./3.*ehor,2.) + 8./3.*cv*F2*ehor*ehor + 5./3.*ehor*ehor*(1.+2./5.*ehor)*F2*F2;
  const double tmp2 = 4.*(cv+F2)*ca*ehor*(1.+4./3.*ehor);
  //const double tmp3 = (cv*cv+3.0*ca*ca)*pow(1.+2.*ehor,3.);
  const double tmp3 = (kGv*kGv+3.0*kGa*kGa)*pow(1.+2.*ehor,3.);
  
  R_nu  = (tmp1+tmp2)/tmp3; // Eq.(22) 
  R_anu = (tmp1-tmp2)/tmp3; //Eq.(22)

  R     = &R_nu;
  Rbar = &R_anu;
  
  return;
}

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N): reacflag = 1 | 2
// Input: omega -> neutrino energy [MeV]
// Output: correction to zeroth (R0) and first Legendre (R1) coefficients of scattering kernel
void WMScatt(const double omega, double* R0, double* R1, const int reacflag) {
  double cv, ca, F2;
  double h0, h1;
  double wm_0, wm_1;

  NucFrmFac(omega, &cv, &ca, &F2, reacflag); //nuclear form factors
  //NucFrmFac(0., &cv_0, &ca_0, &F2_0, reacflag); //nuclear form factors at Q^2=0
 
  // @TODO: evaluate this at compile time  
  if (reacflag == 1) {
    h0 = kHpv*kHpv + 3.*kHpa*kHpa;
    h1 = kHpv*kHpv -    kHpa*kHpa;
  } else if (reacflag == 2) {
    h0 = kHnv*kHnv + 3.*kHna*kHna;
    h1 = kHnv*kHnv -    kHna*kHna;
  }
  
  const double ehor = omega * kMeV/(kMb*kClight*kClight);
  
  /* Low-energy limit derived from Eq.(12) */
  wm_0 = (cv*cv + 3.*ca*ca + 1.5*ehor*ehor*F2*F2 + 2.*ehor*ehor*cv*F2) / h0; //correction to zeroth coefficient
  wm_1 = (cv*cv -    ca*ca -  2.*ehor*ehor*F2*F2 - 4.*ehor*ehor*cv*F2) / h1; //correction to first coefficient

  R0 = &wm_0;
  R1 = &wm_1;

  return;
}
