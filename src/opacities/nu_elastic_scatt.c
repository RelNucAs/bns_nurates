// \file nu_elastic_scatt.c
// \brief Computation of first two Legendre coefficients for neutrino elastic
//        scattering on nucleons
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)

#include <math.h>

#include "opacities.h"
#include "../constants.h"
#include "weak_magnetism/weak_magnetism.h"

// TODO: change names of variables and functions following Google C style

// Definition of parameters
// TODO: decide how to switch on/off corrections to the rates
const int use_WM_sc = 0; // flag for activating weak magnetism (and related) corrections
// 0: not active, 1: active

// Definition of constants
const double hbar = 0.5 * kH / kPi; //[MeV*s]       
const double c0 = (2. * kPi * kGf * kGf) / kH;
const double h0_p = c0 * (kHpv*kHpv + 3.*kHpa*kHpa); // 0th Leg coeff, protons
const double h1_p = c0 * (kHpv*kHpv -    kHpa*kHpa); // 1st Leg coeff, protons
const double h0_n = c0 * (kHnv*kHnv + 3.*kHna*kHna); // 0th Leg coeff, neutrons
const double h1_n = c0 * (kHnv*kHnv -    kHna*kHna); // 1st Leg coeff, neutrons

/* Inputs:
 *      omega    [MeV] : (anti)neutrino energy
 *      mu_ang         : cosine of neutrino polar angle
 *      nb      [cm-3] : baryon number density
 *      temp     [MeV] : temperature
 *      yp             : proton fraction
 *      yn             : neutron fraction
*/


/* Computation of degeneracy parameter eta_NN */
double eta_NN_sc(const double nb, const double temp, const double yN) {
  const double nN = yN*nb;

  if (nN <= 0.) return 0.;  // Enforce zero rates if no nucleons are present
  
  // Linear interpolation between degenerate and nondegnerate limit in Eq.(C37)
  const double eFN = 0.5*hbar*hbar * pow(3.*kPi*kPi*nN,2./3.) / kMb * kMeV; // [MeV]
  const double tmp = 1.5*temp/eFN;
  return nN*tmp/sqrt(1.+tmp*tmp); // [cm-3]
  
  /* Alternative computation: evaluate via numerical derivative (some quick tests showed it's unstable) */
  //dmudrho = ....; //[MeV*cm^3/g]
  //etann = yn / (T*mu*dmudrho); //[1/cm^3]
}

// TODO: elastic scattering is kind of a middle-case between closed-form
//       opacity and opacity integrated from kernel, since it requires
//       only angular intergration. Specific implementation is possibly required
// Legendre expansion of scattering kernel up to l=1
double nu_N_scatt_kern(const double omega, const double mu_ang,
		       const double nb, const double temp, const double yN,
		       const int reacflag) {
  double R0 = 1., R1 = 1.;
  double ker;
	
  const double etaNN = eta_NN_sc(nb, temp, yN); //degeneracy parameter eta_NN, Eq.(C37)

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (use_WM_sc != 0)  WM_scatt(omega, &R0, &R1, reacflag);

  if (reacflag == 1) {
    // Scattering on proton
    ker = h0_p * R0 + h1_p * R1 * mu_ang;
  } else if (reacflag == 2) {
    // Scattering on neutron
    ker = h0_n * R0 + h1_n * R1 * mu_ang;
  }

  return omega * omega * ker; // scattering kernel * enu**2, Eq.(C36)

  // TODO: understand where to add the constant in front
  //       (2.*kPi/kClight)/pow(2.*kPi*hbar*kClight,3.)
}


// Scattering kernel for neutrino-proton scattering
double nu_p_scatt_kern(const double omega, const double mu_ang,
                       const double nb, const double temp, const double yp) {
  return nu_N_scatt_kern(omega, mu_ang, nb, temp, yp, 1);
}

// Scattering kernel for neutrino-neutron scattering
double nu_n_scatt_kern(const double omega, const double mu_ang,
                       const double nb, const double temp, const double yn) {
  return nu_N_scatt_kern(omega, mu_ang, nb, temp, yn, 2);
}

// Sum of scattering kernels for neutrino-nucleon scattering (protons + neutrons)
double nu_N_scatt_kern_tot(const double omega, const double mu_ang,
		           const double nb, const double temp,
			   const double yp, const double yn) {
  return nu_p_scatt_kern(omega, mu_ang, nb, temp, yp) + // proton  contribution
	 nu_n_scatt_kern(omega, mu_ang, nb, temp, yn);  // neutron contribution
}
