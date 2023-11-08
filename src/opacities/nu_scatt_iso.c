// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
/**
 * @file nu_scatt_iso.c
 * @brief Implementation of iso-energetic scattering of neutrinos on nucleons.
 *
 * This file implements iso-energetic scattering on nucleons from [Bruenn 1985](https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
 * with optional inclusion of phase space, recoil and weak magnetism corrections from [Horowitz 2002](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)
 *
 * @note The iso-energetic scattering kernel is same for production and absorption and for all neutrino species.
 */

#include <math.h>
#include <assert.h>

#include "opacities.h"
#include "constants.h"
#include "weak_magnetism.h"
#include "distribution.h"
#include "integration.h"

// Constants for the isoenergetic scattering on nucleons
#define kThreePiSquared2TwoThird 9.570780000627304
//#define kIsoKer (2. * kPi *  kGf *  kGf) / kHbar / (kHbarClight * kHbarClight * kHbarClight * kHbarClight * kHbarClight * kHbarClight) // 2 pi Gf^2 / hbar [MeV cm^6 s^-1] from Eqn. (C36) of Bruenn

/**
 * @fn double EtaNNSc(const double nb, const double temp, const double yN)
 * @brief Computes degeneracy parameter \f$\eta_{NN}\f$ from Eqn. (C37) of Bruenn.
 *
 * @param nb    baryon number density \f$[cm^{-3}]\f$
 * @param temp  temperature \f$[MeV]\f$
 * @param yN    proton/neutron fraction
 * @return      The degeneracy parameter \f$\eta_{NN} [cm^-3]\f$
 */
double EtaNNSc(const double nb, const double temp, const double yN) {
  static const double eF_const = 0.5 * kHbar * kHbar / kMb * kMeV; // constant in Fermi energy calculation

  const double nN = yN * nb; // nucleon (neutron/proton) number density [cm^-3]

  if (nN <= 0.) return 0.;  // Enforce zero rates if no nucleons are present

  // Linear interpolation between degenerate and non-degnerate limit in Eq.(C37)
  const double eFN = eF_const * kThreePiSquared2TwoThird * pow(nN, 2. / 3.); // Fermi energy computation [MeV]
  const double aux = 1.5 * temp / eFN;
  return nN * aux / sqrt(1. + aux * aux); // [cm-3]

}

/**
 * @fn double IsoScattNucleon(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const double yN, const int reacflag)
 * @brief Computes Spectral scattering opacity for iso-energetic scattering on nucleons (protons/neutrons)
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure containing opacity parameters
 * @param eos_pars      structure containing equation of state parameters (needs baryon number density \f$[cm^{-3}]\f$ and temperature \f$[MeV]\f$)
 * @param yN            proton/neutron fraction
 * @param reacflag      choice of nucleon (1: proton scattering 2: neutron scattering)
 * @return              "Eq.(A41)" \f$[MeV cm^{3} s^{-1}]\f$
 */
double IsoScattNucleon(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const double yN, const int reacflag) {
  static const double kIsoKer = (2. * kPi *  kGf *  kGf) / kHbar / (kHClight * kHClight * kHClight);
    
  static const double h0_p = kHpv * kHpv + 3. * kHpa * kHpa;   
  static const double h1_p = kHpv * kHpv - kHpa * kHpa;        
  static const double h0_n = kHnv * kHnv + 3. * kHna * kHna;   
  static const double h1_n = kHnv * kHnv - kHna * kHna;        

  static const double c0_p = kIsoKer * h0_p;  // 0th Legendre coefficient (protons)
  static const double c1_p = kIsoKer * h1_p;  // 1st Legendre coefficient (protons)
  static const double c0_n = kIsoKer * h0_n;  // 0th Legendre coefficient (neutrons)
  static const double c1_n = kIsoKer * h1_n;  // 1st Legendre coefficient (neutrons)

  double R0 = 1., R1 = 1.;
  double leg_0, leg_1;

  const double nb = eos_pars->nb;     // Number baryon density [cm^-3]
  const double temp = eos_pars->temp; // Temperature [MeV]

  const double etaNN = EtaNNSc(nb, temp, yN); // degeneracy parameter eta_NN [cm^-3] from Eqn. (C37) of Bruenn

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (opacity_pars->use_WM_sc) { WMScatt(omega, &R0, &R1, reacflag); }

  if (reacflag == 1) {
    // Scattering on proton
    leg_0 = c0_p * R0;
    leg_1 = c1_p * R1; // [MeV cm^6 s^-1]
  } else if (reacflag == 2) {
    // Scattering on neutron
    leg_0 = c0_n * R0;
    leg_1 = c1_n * R1; // [MeV cm^6 s^-1]
  }

  return etaNN * (leg_0 - leg_1 / 3.); // "Eq.(A41)" [MeV cm^3 s-1]
}

/**
 * @fn double IsoScattProton(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos on protons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV cm{^3} s^{-1}]\f$
 */
double IsoScattProton(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yp, 1);
}

/**
 * @fn double IsoScattNeutron(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos on neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV cm{^3} s^{-1}]\f$
 */
double IsoScattNeutron(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yn, 2);
}

/**
 * @fn double IsoScattTotal(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the total spectral scattering opacity for scattering of neutrinos on protons and neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV cm{^3} s^{-1}]\f$
 */
double IsoScattTotal(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  const double iso_nu_p = IsoScattProton(omega, opacity_pars, eos_pars); // proton contribution
  const double iso_nu_n = IsoScattNeutron(omega, opacity_pars, eos_pars); // neutron contribution

  return iso_nu_p + iso_nu_n;
}


/**
 * @fn double IsoScattLegCoeff(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const int l)
 * @brief Computes Spectral scattering opacity for iso-energetic scattering on nucleons (protons/neutrons)
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure containing opacity parameters
 * @param eos_pars      structure containing equation of state parameters (needs baryon number density \f$[cm^{-3}]\f$ and temperature \f$[MeV]\f$)
 * @param l             order of Legendre coefficient 
 * @return              Legendre coefficient of order l of the isoscattering kernel
 */
double IsoScattLegCoeff(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const int l) {
  assert(l >= 0 && l <= 1); // Legendre order must be either zero or one

  static const double kIsoKer = (2. * kPi *  kGf *  kGf) / kHbar / (kHClight * kHClight * kHClight);
    
  static const double h0_p = kHpv * kHpv + 3. * kHpa * kHpa;   
  static const double h1_p = kHpv * kHpv - kHpa * kHpa;        
  static const double h0_n = kHnv * kHnv + 3. * kHna * kHna;   
  static const double h1_n = kHnv * kHnv - kHna * kHna;        

  static const double c0_p = kIsoKer * h0_p;  // 0th Legendre coefficient (protons)
  static const double c1_p = kIsoKer * h1_p;  // 1st Legendre coefficient (protons)
  static const double c0_n = kIsoKer * h0_n;  // 0th Legendre coefficient (neutrons)
  static const double c1_n = kIsoKer * h1_n;  // 1st Legendre coefficient (neutrons)

  double R0_n = 1., R1_n = 1.;
  double R0_p = 1., R1_p = 1.;
  double leg;

  const double nb = eos_pars->nb;     // Number baryon density [cm^-3]
  const double temp = eos_pars->temp; // Temperature [MeV]
  const double yp = eos_pars->yp;     // Proton fraction
  const double yn = eos_pars->yn;     // Neutron fraction

  // degeneracy parameter from Eqn. (C37) of Bruenn
  const double eta_pp = EtaNNSc(nb, temp, yp); // scattering on protons [cm^-3]
  const double eta_nn = EtaNNSc(nb, temp, yn); // scattering on neutrons [cm^-3]

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (opacity_pars->use_WM_sc) {
     WMScatt(omega, &R0_p, &R1_p, 1);
     WMScatt(omega, &R0_n, &R1_n, 2);
  }

  if (l == 0) {
    leg = eta_pp * c0_p * R0_p + eta_nn * c0_n * R0_n; // [MeV cm^6 s^-1]
  } else {
    leg = eta_pp * c1_p * R1_p + eta_nn * c1_n * R1_n; // [MeV cm^6 s^-1]
  }

  return leg;
}