// \file nu_scatt_iso.c
// \brief Implementation of isoenergetic neutrino scattering on nucleons
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)

// The isoenergetic scattering kernel is equal for production and absorption and is the same for all neutrino types

#include <math.h>

#include "opacities.h"
#include "weak_magnetism/weak_magnetism.h"
//#include "../constants.h"
#include "../bns_nurates.h"
#include "../distribution/distribution.h"
//#include "../functions/functions.h"
#include "../integration/integration.h"

// Definition of constants
#define kIsoKer (2. * kPi *  kGf * kGf) / kHbar // constant in (C36) [MeV cm^6 s-1], 
                                                  // Gf^2/hbar = 1.2E-65 MeV cm^6 s^-1 from (C40)
#define h0_p  kHpv*kHpv + 3.*kHpa*kHpa  // 0th Leg coeff, protons
#define h1_p  kHpv*kHpv -    kHpa*kHpa  // 1st Leg coeff, protons
#define h0_n  kHnv*kHnv + 3.*kHna*kHna  // 0th Leg coeff, neutrons
#define h1_n  kHnv*kHnv -    kHna*kHna  // 1st Leg coeff, neutrons

#define c0_p  kIsoKer * h0_p  // 0th Leg coeff, protons
#define c1_p  kIsoKer * h1_p  // 1st Leg coeff, protons
#define c0_n  kIsoKer * h0_n  // 0th Leg coeff, neutrons
#define c1_n  kIsoKer * h1_n  // 1st Leg coeff, neutrons

/* Inputs:
 *      omega    [MeV] : (anti)neutrino energy (elastic process)
 *      nb      [cm-3] : baryon number density
 *      temp     [MeV] : temperature
 *      yp             : proton fraction
 *      yn             : neutron fraction
*/


/* Computation of degeneracy parameter eta_NN */
double EtaNNsc(const double nb, const double temp, const double yN) {
  static const double three_pi_sqr = 3. * kPi * kPi;
  static const double eF_const = 0.5 * kHbar * kHbar / kMb * kMeV; // constant in Fermi energy calculation

  const double nN = yN*nb;

  if (nN <= 0.) return 0.;  // Enforce zero rates if no nucleons are present
  
  // Linear interpolation between degenerate and nondegnerate limit in Eq.(C37)
  const double eFN = eF_const * pow(three_pi_sqr * nN, 2./3.); // Fermi energy computation [MeV]
  const double tmp = 1.5 * temp / eFN;
  return nN * tmp / sqrt(1. + tmp*tmp); // [cm-3]
  
  /* Alternative computation: evaluate via numerical derivative (some quick tests showed it's unstable) */
  //dmudrho = ....; //[MeV*cm^3/g]
  //etann = yn / (T*mu*dmudrho); //[1/cm^3]
}

// Spectral scattering opacity for isoenergetic scattering on nucleons
double IsoScattNucleon(double omega, OpacityParams *opacity_pars,
                       MyEOSParams *eos_pars,
                       const double yN, const int reacflag) {
  double R0 = 1., R1 = 1.;
  double leg_0, leg_1;
	   
  const double nb   = eos_pars->nb;   // Number baryon density [cm-3]
  const double temp = eos_pars->temp; // Temperature [MeV]

  const double etaNN = EtaNNsc(nb, temp, yN); // degeneracy parameter eta_NN, Eq.(C37)

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (opacity_pars->use_WM_sc)  WMScatt(omega, &R0, &R1, reacflag);

  if (reacflag == 1) {
    // Scattering on proton
    leg_0 = h0_p * R0;
    leg_1 = h1_p * R1; // [MeV cm^3 s-1]
  } else if (reacflag == 2) {
    // Scattering on neutron
    leg_0 = h0_n * R0;
    leg_1 = h1_n * R1; // [MeV cm^3 s-1]
  }
  
  return omega * omega * etaNN * (leg_1 / 3. - leg_0); // "Eq.(A41)" [MeV^3 cm^3 s-1]
}

// Neutrino-proton contribution
double IsoScattProton(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yp, 1);
}

// Neutrino-neutron contribution
double IsoScattNeutron(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  return IsoScattNucleon(omega, opacity_pars, eos_pars, eos_pars->yn, 2);
}

// Sum of the two contributions for neutrino-nucleon scattering (protons + neutrons)
double IsoScattTotal(double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  const double iso_nu_p = IsoScattProton(omega , opacity_pars, eos_pars); // proton contribution
  const double iso_nu_n = IsoScattNeutron(omega, opacity_pars, eos_pars); // neutron contribution

  //MyKernel elastic_kernel = {.absorption_e = nu_p_ker.absorption_e + nu_n_ker.absorption_e,
  //                           .production_e = nu_p_ker.production_e + nu_n_ker.production_e,
  //                           .absorption_x = nu_p_ker.absorption_x + nu_n_ker.absorption_x,
  //                           .production_x = nu_p_ker.production_x + nu_n_ker.production_x};

  return iso_nu_p + iso_nu_n;
}

// NuNumberScatteringIntegrand function
// integrand for the computation of the number opacity
// scattering coefficient
// (constants are added after the integration)
double NuNumberScatteringIntegrand(double *x, void *p) {
  GreyOpacityParams *grey_pars = (GreyOpacityParams *) p;
  
  const double iso_scatt = IsoScattTotal(x[0], &grey_pars->opacity_pars, &grey_pars->eos_pars);

  return x[0] * x[0] * iso_scatt * TotalNuF(x[0], &grey_pars->distr_pars);
}

// NuEnergyScatteringIntegrand function
// integrand for the computation of the energy opacity
// scattering coefficient
// (constants are added after the integration)
double NuEnergyScatteringIntegrand(double *x, void *p) {
  return x[0] * NuNumberScatteringIntegrand(x, p); 
}

// Computation of scattering coefficients
SourceCoeffs IsoScattCoeffs(GreyOpacityParams *grey_pars) {
  double iso_n, iso_j;
  
  SourceCoeffs out;

  MyFunction integrand;

  integrand.dim = 1;
  integrand.params = grey_pars;

  // @TODO: compute quadrature points once for all
  const int np = 32;

  //static MyQuadrature quadrature_default = {.type=kGauleg, .dim=1, .nx=32, .ny=1, .nz=1, .alpha=0., .x1=0., .x2=1., .y1=-42., .y2=-42., .z1=-42., .z2=-42.};
  MyQuadrature quad = quadrature_default;

  quad.nx = np;
  quad.type = kGauleg;
  quad.dim = 1;
  quad.x1 = 0.;
  quad.x2 = 1.;

  GaussLegendreMultiD(&quad);

  double n = NuNumber(&grey_pars->distr_pars);
  double J = NuEnergy(&grey_pars->distr_pars);

  double s = 2. * grey_pars->distr_pars.eta_t * grey_pars->distr_pars.temp_t;

  integrand.function = &NuNumberScatteringIntegrand;
  iso_n = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &integrand, s) / n / pow(kH * kClight, 3.);
  
  out.R_nue  = iso_n;
  out.R_anue = iso_n;
  
  //out.R_num  = iso_n;
  //out.R_anum = iso_n;
  
  out.R_nux = iso_n;

  integrand.function = &NuEnergyScatteringIntegrand;
  iso_j = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &integrand, s) / J / pow(kH * kClight, 3.);

  out.Q_nue  = iso_j;
  out.Q_anue = iso_j;
  
  //out.Q_num  = iso_j;
  //out.Q_anum = iso_j;
  
  out.Q_nux = iso_j;
  
  return out;
}

