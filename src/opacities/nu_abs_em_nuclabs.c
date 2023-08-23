// \file nu_abs_em_nuclabs.c
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
#include "weak_magnetism/weak_magnetism.h"
#include "../constants.h"
#include "../functions/functions.h"
#include "../integration/integration.h"


// @TODO: align with constant definitions for opacities coming from kernel integration
// Definition of constants
const double g0 = 0.5 * kH * kClight / kPi;
const double g1 = (kGf*kGf/kPi) / (g0 * g0 * g0 * g0); // (GF*GF/pi) / pow(hbar*c,4.),
                                                    // constant in front of
                                                    // Eq.(C13,C15,C19,C20)
const double g2 = kGv*kGv+3.*kGa*kGa; // constant in Eq.(C13,C15,C19,C20)
const double g3 = g1 * g2;
const double mu_thres = 1.E-02;   // mu_hat threshold value in EtaNNAbs evaluation

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
double EtaNNAbs(const double n_in, const double n_out,
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
double EtaNP(const double nn, const double np,
	      const double mu_hat,
	      const double temp) {
  return EtaNNAbs(nn, np, mu_hat, temp);
}

// Nucleon phase space integration for X + p -> X + n 
double EtaPN(const double nn, const double np,
              const double mu_hat,
	      const double temp) {
  return EtaNNAbs(np, nn, -mu_hat, temp);
}

// Neutrino absorption on neutron (nul + n -> l- + p)
// Antineutrino absorption on neutron (anul + p -> l+ + n)
void AbsOpacitySingleLep(OpacityParams* opacity_pars,
                         MyEOSParams *eos_pars,
                         const double mLep, const double muLep,
                         double *out) {   // out[0] --> nu_abs
                                          // out[1] --> nu_em
                                          // out[2] --> anu_abs
                                          // out[3] --> anu_em
  double Qprime, mu_np;
  double etanp, etapn;
  double E_e, E_p;
  double tmp, fd_e, fd_p;
  double R = 1., Rbar = 1.;
  double dU = 0.;

  const double omega  = opacity_pars->omega; // (Anti)neutrino energy [MeV]

  const double nb     = eos_pars->nb;      // Number baryon density [cm-3]
  const double temp   = eos_pars->temp;    // Temperature [MeV]
  const double yp     = eos_pars->yp;      // Proton fraction
  const double yn     = eos_pars->yn;      // Neutron fraction
  const double mu_p   = eos_pars->mu_p;    // Proton chemical potential [MeV]
  const double mu_n   = eos_pars->mu_n;    // Neutron chemical potential [MeV]

  const double nn = nb*yn; // Neutron number density [cm-3]
  const double np = nb*yp; // Proton number density  [cm-3]

  // Neutron minus proton chem. potentials (corrected for the mass difference)
  const double mu_hat = mu_n - mu_p - kQ; // [MeV] //kQ needed if mu_p and mu_n are relativistic chem. potentials
  
  // Nucleon interaction correction to chemical potential
  if (opacity_pars->use_dU) dU = eos_pars->dU; // [MeV]

  Qprime = kQ + dU;     // [MeV], Eq.(79) in Hempel
  mu_np  = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel
                                           
  etanp = EtaNP(nn,np,mu_np,temp); // Eq. (C14)
  etapn = EtaPN(nn,np,mu_np,temp);

  E_e = omega + Qprime;  // Electron energy [MeV]
  E_p = omega - Qprime;  // Positron energy [MeV]

  // Phase space, recoil and weak magnetism correction
  if (opacity_pars->use_WM_ab) WMAbsEm(omega, &R, &Rbar);

  // @TODO: Leonardo -> look for smartest choice to check kinematics constraints
  // Check kinematics constraint for neutrino absorption
  if (E_e-mLep > 0.) {
    tmp  = R * kClight * g3 * E_e * E_e * sqrt(1.-pow(mLep/E_e,2.)); // remove c to get output in cm-1
    fd_e = FermiDistr(E_e,temp,muLep);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)
    out[0] = etanp * tmp * (1 - fd_e); // Neutrino absopritvity [s-1], Eq.(C13)
    out[1] = etapn * tmp * fd_e;       // Neutrino emissivity   [s-1], Eq.(C15)
  }

  // Check kinematics constraint for antineutrino absorption
  if (E_p-mLep > 0.) {
    tmp  = Rbar * kClight * g3 * E_p * E_p * sqrt(1.-pow(mLep/E_p,2.)); // remove c to get output in cm-1
    fd_p = FermiDistr(E_p,temp,-muLep);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)
    out[2] = etapn * tmp * (1. - fd_p); // Antineutrino absopritvity [s-1], Eq.(C19)
    out[3] = etanp * tmp * fd_p;        // Antineutrino emissivity   [s-1], Eq.(C20)
  }
  /* Emissivity from detailed balance (NOT TESTED) */
  //double em = ab * c * exp(-(omega-(mu_l-mu_hat-delta_np))/temp);
  //double em = ab * c * exp(-(omega_bar-(mu_hat+delta_np-mu_l))/temp);

  return;
}

MyOpacity AbsOpacity(OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  MyOpacity MyOut = {0.0}; // initialize to zero

  // Electron (anti)neutrino
  double el_out[4] = {0.0};
  AbsOpacitySingleLep(opacity_pars, eos_pars, kMe, eos_pars->mu_e, el_out);

  MyOut.ab_nue   = el_out[0];
  MyOut.em_nue   = el_out[1];
  MyOut.ab_anue  = el_out[2];
  MyOut.em_anue  = el_out[3];

  // Uncomment the following when considering also muons 
  // // Muon (anti)neutrino
  //double mu_out[4] = {0.0};
  //AbsOpacitySingleLep(opacity_pars, eos_pars, kMmu, eos_pars->mu_mu, mu_out);

  //MyOut.ab_num   = mu_out[0];
  //MyOut.em_num   = mu_out[1];
  //MyOut.ab_anum  = mu_out[2];
  //MyOut.em_anum  = mu_out[3];
  
  return MyOut;
}

// Stimulated absorption versions
MyOpacity StimAbsOpacity(OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  MyOpacity abs_opacity = AbsOpacity(opacity_pars, eos_pars);
  MyOpacity stim_abs_opacity = {.ab_nue  = abs_opacity.ab_nue  + abs_opacity.em_nue,
                                .em_nue  = abs_opacity.em_nue,
                                .ab_anue = abs_opacity.ab_anue + abs_opacity.em_anue,
                                .em_anue = abs_opacity.em_anue,
                                //.ab_num  = abs_opacity.ab_num  + abs_opacity.em_num,
                                //.em_num  = abs_opacity.em_num,
                                //.ab_anum = abs_opacity.ab_anum + abs_opacity.em_anum,
                                //.em_anum = abs_opacity.em_anum,
                                .ab_nux  = 0.,
                                .em_nux  = 0.};
                                
  return stim_abs_opacity;
}

// NueAbsNumberEmissivityIntegrand function
// integrand for the computation of the number emission
// coefficient for electron-type neutrino 
// (v**3 * j(nu), constants are added after the
// integration)
double NueAbsNumberEmissivityIntegrand(double *x, void *p) {
  MyEOSParams *eos_pars = (MyEOSParams *) p;
  
  // @TODO: modify way in which corrections can be included
  OpacityParams opacity_pars = {.omega = x[0],
                                .use_dU = 0,
                                .use_WM_ab = 0};
                                
  MyOpacity out = AbsOpacity(&opacity_pars, eos_pars);

  return x[0] * x[0] * out.em_nue;
}

// NueAbsEnergyEmissivityIntegrand function
// integrand for the computation of the energy emission
// coefficient for electron-type neutrino 
// (v**3 * j(nu), constants are added after the
// integration)
double NueAbsEnergyEmissivityIntegrand(double *x, void *p) {
  return x[0] * NueAbsNumberEmissivityIntegrand(x, p); 
  //MyEOSParams *eos_pars = (MyEOSParams *) p;
  
  // @TODO: modify way in which corrections can be included
  //OpacityParams opacity_pars = {.omega = x[0],
  //                              .use_dU = 0,
  //                              .use_WM_ab = 0};
                                
  //MyOpacity out = AbsOpacity(&opacity_pars, eos_pars);

  //return x[0] * x[0] * x[0] * out.em_nue;
}


// ANueAbsNumberEmissivityIntegrand function
// integrand for the computation of the number emission
// coefficient for electron-type antineutrino 
// (v**2 * j(nu), constants are added after the
// integration)
double ANueAbsNumberEmissivityIntegrand(double *x, void *p) {
  MyEOSParams *eos_pars = (MyEOSParams *) p;
  
  // @TODO: modify way in which corrections can be included
  OpacityParams opacity_pars = {.omega = x[0],
                                .use_dU = 0,
                                .use_WM_ab = 0};
                                
  MyOpacity out = AbsOpacity(&opacity_pars, eos_pars);

  return x[0] * x[0] * out.em_anue;
}

// ANueAbsEnergyEmissivityIntegrand function
// integrand for the computation of the energy emission
// coefficient for electron-type antineutrino 
// (v**3 * j(nu), constants are added after the
// integration)
double ANueAbsEnergyEmissivityIntegrand(double *x, void *p) {
  return x[0] * ANueAbsNumberEmissivityIntegrand(x, p); 
}

SourceCoeffs NuclAbsEmissionCoeffs(MyEOSParams *eos_pars) {
  SourceCoeffs out;

  MyFunction em_integrand;

  em_integrand.dim = 1;
  em_integrand.params = eos_pars;

  // @TODO: compute quadrature points once for all
  const int n = 32;

  //static MyQuadrature quadrature_default = {.type=kGauleg, .dim=1, .nx=32, .ny=1, .nz=1, .alpha=0., .x1=0., .x2=1., .y1=-42., .y2=-42., .z1=-42., .z2=-42.};
  MyQuadrature quad = quadrature_default;

  quad.nx = n;
  quad.type = kGauleg;
  quad.dim = 1;
  quad.x1 = 0.;
  quad.x2 = 1.;

  GaussLegendreMultiD(&quad);

  double s = eos_pars->mu_e - kQ;

  em_integrand.function = &NueAbsNumberEmissivityIntegrand;
  out.R_nue  = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &em_integrand, s) / pow(kH * kClight, 3.);

  em_integrand.function = &ANueAbsNumberEmissivityIntegrand;
  out.R_anue = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &em_integrand, s) / pow(kH * kClight, 3.);
  
  //out.R_num = ...
  //out.R_anum = ...
  
  out.R_nux = 0.;

  em_integrand.function = &NueAbsEnergyEmissivityIntegrand;
  out.Q_nue  = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &em_integrand, s) / pow(kH * kClight, 3.);

  em_integrand.function = &ANueAbsEnergyEmissivityIntegrand;
  out.Q_anue = 4. * kPi * GaussLegendreIntegrateZeroInf(&quad, &em_integrand, s) / pow(kH * kClight, 3.);
  
  //out.Q_num  = ...
  //out.Q_anum = ...
 
  out.Q_nux = 0;
  
  return out;
}

// @TODO: move the following somewhere else
// Theta step function 
double theta(const double x) {
  if (x < 0.) return 0.;
  else        return 1.;
}
