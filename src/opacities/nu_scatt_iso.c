// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
/**
 * @file nu_scatt_iso.c
 * @brief Implementation of iso-energetic scattering of neutrinos with nucleons.
 *
 * This file implements iso-energetic scattering on nucleons from [Bruenn 1985](https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
 * with optional inclusion of phase space, recoil and weak magnetism corrections from [Horowitz 2002](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)
 *
 * @note The iso-energetic scattering kernel is same for production and absorption and for all neutrino species.
 */

#include <math.h>

#include "opacities.h"
#include "weak_magnetism/weak_magnetism.h"
#include "../distribution/distribution.h"
#include "../integration/integration.h"

// Definition of local constants
#define kIsoKer (2. * kPi *  kGf * kGf) / kHbar // 2 pi Gf^2 / hbar [MeV cm^6 s^-1] from Eqn. (C36) of Bruenn

#define h0_p  kHpv*kHpv + 3.*kHpa*kHpa          // 0th Legendre coefficient (protons)
#define h1_p  kHpv*kHpv -    kHpa*kHpa          // 1st Legedndre coefficient (protons)
#define h0_n  kHnv*kHnv + 3.*kHna*kHna          // 0th Legendre coefficient (neutrons)
#define h1_n  kHnv*kHnv -    kHna*kHna          // 1st Legendre coefficient (neutrons)

#define c0_p  kIsoKer * h0_p                    // 0th Legendre coefficient (protons)
#define c1_p  kIsoKer * h1_p                    // 1st Legendre coefficient (protons)
#define c0_n  kIsoKer * h0_n                    // 0th Legendre coefficient (neutrons)
#define c1_n  kIsoKer * h1_n                    // 1st Legendre coefficient (neutrons)

/**
 * @fn double EtaNNsc(const double nb, const double temp, const double yN)
 * @brief Computes degeneracy parameter \f$\eta_{NN}\f$ from Eqn. (C37) of Bruenn.
 *
 * @param nb    baryon number density \f$[cm^{-3}]\f$
 * @param temp  temperature \f$[MeV]\f$
 * @param yN    proton/neutron fraction
 * @return      The degeneracy parameter \f$\eta_{NN} [cm^-3]\f$
 */
double EtaNNsc(const double nb, const double temp, const double yN) {

  static const double three_pi_sqr = 3. * kPi * kPi;
  static const double eF_const = 0.5 * kHbar * kHbar / kMb * kMeV; // constant in Fermi energy calculation

  const double nN = yN * nb;

  if (nN <= 0.) return 0.;  // Enforce zero rates if no nucleons are present

  // Linear interpolation between degenerate and nondegnerate limit in Eq.(C37)
  const double eFN = eF_const * pow(three_pi_sqr * nN, 2. / 3.); // Fermi energy computation [MeV]
  const double tmp = 1.5 * temp / eFN;
  return nN * tmp / sqrt(1. + tmp * tmp); // [cm-3]

}

/**
 * @fn double IsoScattNucleon(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const double yN, const int reacflag)
 * @brief Computes Spectral scattering opacity for iso-energetic scattering on nucleons (protons/neutrons)
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure containing opacity parameters
 * @param eos_pars      structure containing equation of state parameters (needs baryon number density \f$[cm^{-3}]\f$ and temperature \f$[MeV]\f$)
 * @param yN            proton/neutron fraction
 * @param reacflag      choice of nucleon (1: proton scattering 2: neutron scattering)
 * @return              "Eq.(A41)" \f$[MeV^{3} cm{^3} s^{-1}]\f$
 */
double IsoScattNucleon(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const double yN, const int reacflag) {

  double R0 = 1., R1 = 1.;
  double leg_0, leg_1;

  const double nb = eos_pars->nb;   // Number baryon density [cm-3]
  const double temp = eos_pars->temp; // Temperature [MeV]

  const double etaNN = EtaNNsc(nb, temp, yN); // degeneracy parameter eta_NN [cm^-3] from Eqn. (C37) of Bruenn

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (opacity_pars->use_WM_sc) { WMScatt(omega, &R0, &R1, reacflag); }

  if (reacflag == 1) {
    // Scattering on proton
    leg_0 = c0_p * R0;
    leg_1 = c1_p * R1; // [MeV cm^3 s-1]
  } else if (reacflag == 2) {
    // Scattering on neutron
    leg_0 = c0_n * R0;
    leg_1 = c1_n * R1; // [MeV cm^3 s-1]
  }

  return omega * omega * etaNN * (leg_1 / 3. - leg_0); // "Eq.(A41)" [MeV^3 cm^3 s-1]
}

/**
 * @fn double IsoScattProton(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos on protons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV^{3} cm{^3} s^{-1}]\f$
 */
double IsoScattProton(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yp, 1);
}

/**
 * @fn double IsoScattNeutron(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the spectral scattering opacity for scattering of neutrinos on neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV^{3} cm{^3} s^{-1}]\f$
 */
double IsoScattNeutron(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yn, 2);
}

/**
 * @fn double IsoScattTotal(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars)
 * @brief Computes the total spectral scattering opacity for scattering of neutrinos on protons and neutrons
 * @param omega         neutrino energy \f$[MeV]\f$
 * @param opacity_pars  structure for opacity parameters
 * @param eos_pars      structure for equation of state parameters
 * @return              "Eq.(A41)" \f$[MeV^{3} cm{^3} s^{-1}]\f$
 */
double IsoScattTotal(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  const double iso_nu_p = IsoScattProton(omega, opacity_pars, eos_pars); // proton contribution
  const double iso_nu_n = IsoScattNeutron(omega, opacity_pars, eos_pars); // neutron contribution

  return iso_nu_p + iso_nu_n;
}


// @TODO: generalize all the following functions with structures containing all nu species

// NuEnergyScatteringIntegrand function
// integrand for the computation of the energy opacity
// scattering coefficient
// (constants are added after the integration)
MyQuadratureIntegrand NuEnergyScatteringIntegrand(double *x, void *p) {
  GreyOpacityParams *grey_pars = (GreyOpacityParams *) p;

  const double iso_scatt = IsoScattTotal(x[0], &grey_pars->opacity_pars, &grey_pars->eos_pars);

  MyQuadratureIntegrand result;

  result.integrand[0] = x[0] * x[0] * x[0] * iso_scatt * TotalNuF(x[0], &grey_pars->distr_pars, 0);
  result.integrand[1] = x[0] * x[0] * x[0] * iso_scatt * TotalNuF(x[0], &grey_pars->distr_pars, 1);
  result.integrand[2] = x[0] * x[0] * x[0] * iso_scatt * TotalNuF(x[0], &grey_pars->distr_pars, 2);

  return result;
}

// Computation of scattering coefficients
SourceCoeffs IsoScattCoeffs(GreyOpacityParams *grey_pars) {
  MyQuadratureIntegrand iso;

  SourceCoeffs out;

  MyFunctionMultiD integrand;

  integrand.dim = 1;
  integrand.params = grey_pars;
  integrand.my_quadrature_integrand.n = 3;

  MyQuadrature quad = quadrature_default; //{.type=kGauleg, .dim=1, .nx=32, .ny=1, .nz=1, .alpha=0., .x1=0., .x2=1., .y1=-42., .y2=-42., .z1=-42., .z2=-42.};

  GaussLegendreMultiD(&quad);

  double s = 2. * grey_pars->distr_pars.eta_t[0] * grey_pars->distr_pars.temp_t[0]; //@TODO: we are choose the same s for all species

  out.R_nue = 0.;
  out.R_anue = 0.;

  //out.R_num  = 0.;
  //out.R_anum = 0.;

  out.R_nux = 0.;

  integrand.function = &NuEnergyScatteringIntegrand;
  iso = GaussLegendreIntegrate1D(&quad, &integrand, s);

  out.Q_nue = 4. * kPi * iso.integrand[0] / grey_pars->m1_pars.J[0] / kClight;
  out.Q_anue = 4. * kPi * iso.integrand[1] / grey_pars->m1_pars.J[1] / kClight;

  //out.Q_num  = iso;
  //out.Q_anum = iso;

  out.Q_nux = 4. * kPi * iso.integrand[2] / grey_pars->m1_pars.J[2] / kClight;

  return out;
}

