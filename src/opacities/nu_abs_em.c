// \file nu_abs_em.c
// \brief Computation of emissivity and absorptivity for neutrino absorption on neutron
//        (and inverse reaction) and for antineutrino absorption on proton (and inverse reaction)
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
//        
// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001) 
//
// Possible inclusion of nucleon interation correction as in 
// Hempel, 2015 (https://journals.aps.org/prc/abstract/10.1103/PhysRevC.91.055807)

#include <math.h>

#include "opacities.h"
#include "../constants.h"
#include "../functions/functions.h"
#include "weak_magnetism/weak_magnetism.h"

// @TODO: change names of variables and functions following Google C style


// Definition of parameters
// @TODO: decide how to switch on/off corrections to the rates
const int use_WM_ab = 0; // flag for activating weak magnetism (and related) corrections
const int use_dU = 1;    // flag for activating dU correction
// 0: not active, 1: active


// Definition of constants
const double g0 = 0.5 * kH * kClight / kPi;
const double g1 = (kGf*kGf/kPi) / (g0 * g0 * g0 * g0); // (GF*GF/pi) / pow(hbar*c,4.),
                                                    // constant in front of
                                                    // Eq.(C13,C15,C19,C20)
const double g2 = kGv*kGv+3.*kGa*kGa; // constant in Eq.(C13,C15,C19,C20)
const double g3 = g1 * g2;
const double mu_thres = 1.E-02;   // mu_hat threshold value in eta_NN evaluation

/* Inputs:
 * 	omega    [MeV] : (anti)neutrino energy
 * 	nb      [cm-3] : baryon number density
 * 	temp     [MeV] : temperature
 * 	lep_mass [MeV] : lepton mass (electron/muon) 
 * 	yp             : proton fraction
 * 	yn             : neutron fraction
 * 	mu_l     [MeV] : lepton chemical potential, rest mass included
 * 	mu_hat   [MeV] : neutron-proton chemical potential difference, rest mass NOT included
 * 	dU       [MeV] : nuclear interaction correction on mu_hat
*/


// Generic function for nucleon phase space integration
double eta_NN_abs(const double n_in, const double n_out,
	          const double mu_hat,
		  const double temp) {
  if (n_in == 0.) return 0.; // enforce zero rates if no nucleons 
                             // available in the initial state

  if (fabs(mu_hat) < mu_thres) return n_in; // backup solution if mu_hat
                                            // too small

  const double etanp = (n_out-n_in) / (exp(-mu_hat/temp)-1.); // Eq.(C14), [cm-3]

  if (etanp < 0.) return n_in; //backup if etanp is negative
	
  return etanp;
}

// Nucleon phase space integration for X + n -> X + p 
double eta_np(const double nn, const double np,
	      const double mu_hat,
	      const double temp) {
  return eta_NN_abs(nn, np, mu_hat, temp);
}

// Nucleon phase space integration for X + p -> X + n 
double eta_pn(const double nn, const double np,
              const double mu_hat,
	      const double temp) {
  return eta_NN_abs(np, nn, -mu_hat, temp);
}

// @TODO: optimized function which computes abs on p and n at the same time

// Neutrino absorption on neutron (nul + n -> l- + p)
MyOpacity nu_n_abs(const double omega,
                   const double nb, const double temp,
		   const double lep_mass,
		   const double yp, const double yn,
		   const double mu_l, const double mu_hat,
		   const double deltaU) {
  const double nn = nb*yn; // neutron number density [cm-3]
  const double np = nb*yp; // proton number density  [cm-3]
  double Qprime, E, mu_np;
  double dU = 0.;
  double R = 1.;
  double tmp, fd;

  MyOpacity out = {0.0, 0.0}; // {em, ab}

  // Nucleon interaction correction to chemical potential
  if (use_dU == 1) dU = deltaU; // [MeV]
 
  Qprime = kQ + dU;     // [MeV], Eq.(79) in Hempel
  E = omega + Qprime;  // [MeV]
  mu_np = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

  // Phase space, recoil and weak magnetism correction
  if (use_WM_ab == 1) R = WM_nue_abs(omega);

  // Check kinematics constraint for the reaction
  if (E-lep_mass > 0.) {
    tmp = R * kClight * g3 * E * E * sqrt(1.-pow(lep_mass/E,2.));
    fd  = FermiDistr(mu_l,E,temp);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)
    
    // Absoprtivity [s-1]
    out.ab = eta_np(np,nn,mu_np,temp) * tmp * (1 - fd); // Eq.(C13)
		
    // Emissivity [s-1]
    out.em = eta_pn(np,nn,mu_np,temp) * tmp * fd;       // Eq.(C15)
  }

  /* Emissivity from detailed balance (NOT TESTED) */
  //double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/temp);

  return out;
}

// Antineutrino absorption on neutron (anul + p -> l+ + n)
MyOpacity nu_p_abs(const double omega,
                   const double nb, const double temp,
		   const double lep_mass,
		   const double yp, const double yn,
		   const double mu_l, const double mu_hat,
		   const double deltaU) {
  const double nn = nb*yn; //neutron number density [cm-3]
  const double np = nb*yp; //proton number density  [cm-3]
  double Qprime, E, mu_np;
  double dU = 0.;
  double Rbar = 1.;
  double tmp, fd;

  MyOpacity out = {0.0, 0.0}; // {em, ab}

  // Nucleon interaction correction to chemical potential
  if (use_dU == 1) dU = deltaU; // [MeV]
	
  Qprime = kQ + dU;     // [MeV], Eq.(79) in Hempel
  E = omega - Qprime;  // [MeV]
  mu_np = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

  // Phase space, recoil and weak magnetism correction
  if (use_WM_ab == 1)  Rbar = WM_anue_abs(omega);
	
  // Check kinematics constraint for the reaction
  if (E-lep_mass > 0.) {
    tmp = Rbar * kClight * g3 * E * E * sqrt(1.-pow(lep_mass/E,2.));
    fd  = FermiDistr(-mu_l,E,temp);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)

    // Absorptivity [s-1]
    out.ab = eta_pn(np,nn,mu_np,temp) * tmp * (1. - fd); // Eq.(C19)

    // Emissivity [s-1]
    out.em = eta_np(np,nn,mu_np,temp) * tmp * fd;          // Eq.(C20)
  }

  /* Emissivity from detailed balance (NOT TESTED) */
  //double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/temp);

  return out;
}


// Stimulated absoption versions
MyOpacity nu_n_abs_stim(const double omega,
                        const double nb, const double temp,
		       	const double lep_mass,
		       	const double yp, const double yn,
		       	const double mu_l, const double mu_hat,
		       	const double deltaU) {
  MyOpacity in = nu_n_abs(omega, nb, temp, lep_mass, yp, yn, mu_l, mu_hat, deltaU);
  MyOpacity out = {in.em, in.em + in.ab}; // stimulated absorption

  return out;
}

MyOpacity nu_p_abs_stim(const double omega,
                        const double nb, const double temp,
                        const double lep_mass,
                        const double yp, const double yn,
                        const double mu_l, const double mu_hat,
                        const double deltaU) {
  MyOpacity in = nu_p_abs(omega, nb, temp, lep_mass, yp, yn, mu_l, mu_hat, deltaU);
  MyOpacity out = {in.em, in.em + in.ab}; // stimulated absorption

  return out;
}


// @TODO: move the following somewhere else
// Theta step function 
double theta(const double x) {
  if (x < 0.) return 0.;
  else        return 1.;
}
