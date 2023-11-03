//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
// \file nu_abs_em_beta.c
// \brief Computes emissivity and absorptivity for neutrino absorption on neutron
//        and for antineutrino absorption on proton and their inverse reactions
//
// Reference: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)
//
// Also include:
//           Possible inclusion of phase space, recoil, weak magnetism correction as in
//           Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)
//
//           Possible inclusion of nucleon interation correction as in
//           Hempel, 2015 (https://journals.aps.org/prc/abstract/10.1103/PhysRevC.91.055807)

#include <math.h>

#include "opacities.h"
#include "weak_magnetism.h"
#include "constants.h"
#include "distribution.h"
#include "functions.h"
#include "integration.h"

// Constants for beta reactions
#define mu_thres 1.0E-02  // mu_hat threshold value in EtaNNAbs evaluation

// ----------------------------------------------------------------------------------

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

/* Generic function for nucleon phase space integration
 * Computes absorptivity from Eqn. (C14) of Bruenn
 *
 * Note:
 *
 * mu_hat [MeV] is neutron-proton chemical potential difference, rest mass NOT included
 * This function is never called directly.
 *
 *  Output: j [cm^-3]
 */
double EtaNNAbs(const double n_in, const double n_out, const double mu_hat, const double temp) {

  // if no nucleons, enforce zero rates
  if (n_in == 0.) return 0.;

  // if mu_hat too small, neglect nucleon degeneracy as backup
  if (fabs(mu_hat) < mu_thres) return n_in;

  // Eq.(C14), [cm-3]
  const double etanp = (n_out - n_in) / (SafeExp(-mu_hat / temp) - 1.);

  // backup if etanp is negative
  if (etanp < 0.) return n_in;

  return etanp;
}

/* Perform phase space integration for:
 *    X + n -> X + p
 * Inputs:
 *      nn [MeV], np [MeV], mu_hat [MeV], temp [MeV]
 * Output:
 *      Eta_np from Bruenn (C14) [cm^-3]
 */
double EtaNP(const double nn, const double np, const double mu_hat, const double temp) {
  return EtaNNAbs(nn, np, mu_hat, temp);
}

/* Perform phase space integration for:
 *    X + p -> X + n
 * Inputs:
 *      nn [MeV], np [MeV], mu_hat [MeV], temp [MeV]
 * Output:
 *      Eta_pn [cm^-3]
 */
double EtaPN(const double nn, const double np, const double mu_hat, const double temp) {
  return EtaNNAbs(np, nn, -mu_hat, temp);
}


// @TODO: add effective mass correction to the rates

/* Compute opacities for:
 *
 * 1. Neutrino absorption on neutron (nu + n -> l- + p)
 * 2. Antineutrino absorption on neutron (anul + p -> l+ + n)
 *
 * Inputs:
 *    omega [MeV], mLep [MeV], muLep [MeV]
 * Outputs: j_x and 1/lambda_x for electron neutrino and electron antineutrino
 *    out[0]: j_nue [s^-1], out[1]: 1/lamda_nue [s^-1], out[2]: j_anue [s^-1], out[3]: 1/lamda_anue [s^-1]
 */
void AbsOpacitySingleLep(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars, const double mLep, const double muLep, double *out) {
  static const double k1 = kClight * kGf * kGf / kPi / (kHbarClight * kHbarClight * kHbarClight * kHbarClight);
  static const double k2 = kGv * kGv + 3. * kGa * kGa; 
  static const double kAbsEmConst = k1 * k2; // constant in Eq.(C13,C15,C19,C20)

  double Qprime, mu_np;
  double etanp, etapn;
  double E_e, E_p;
  double aux, fd_e, fd_p;
  double R = 1., Rbar = 1.;
  double dU = 0.;

  const double nb = eos_pars->nb;         // Number baryon density [cm-3]
  const double temp = eos_pars->temp;     // Temperature [MeV]
  const double yp = eos_pars->yp;         // Proton fraction
  const double yn = eos_pars->yn;         // Neutron fraction
  const double mu_p = eos_pars->mu_p;     // Proton chemical potential [MeV]
  const double mu_n = eos_pars->mu_n;     // Neutron chemical potential [MeV]

  const double nn = nb * yn;              // Neutron number density [cm-3]
  const double np = nb * yp;              // Proton number density  [cm-3]

  // Neutron minus proton chem. potentials (corrected for the mass difference)
  const double mu_hat = mu_n - mu_p - kQ; // [MeV] //kQ needed if mu_p and mu_n are relativistic chem. potentials

  // Nucleon interaction correction to chemical potential
  if (opacity_pars->use_dU) dU = eos_pars->dU; // [MeV]

  Qprime = kQ + dU;     // [MeV], Eq.(79) in Hempel
  mu_np = mu_hat - dU; // [MeV], Eq.(80,86) in Hempel

  etanp = EtaNP(nn, np, mu_np, temp); // Eq. (C14)
  etapn = EtaPN(nn, np, mu_np, temp);

  E_e = omega + Qprime;  // Electron energy [MeV]
  E_p = omega - Qprime;  // Positron energy [MeV]

  // Phase space, recoil and weak magnetism correction
  if (opacity_pars->use_WM_ab) WMAbsEm(omega, &R, &Rbar);

  // @TODO: Leonardo -> look for smartest choice to check kinematics constraints
  // Check kinematics constraint for neutrino absorption
  if (E_e - mLep > 0.) {
    aux = R *  kAbsEmConst * E_e * E_e * sqrt(1. - mLep * mLep / (E_e * E_e)); // remove c to get output in cm-1
    fd_e = FermiDistr(E_e, temp, muLep);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)
    out[1] = etapn * aux * fd_e;       // Neutrino emissivity   [s-1], Eq.(C15)
    out[0] = out[1] * SafeExp((omega - (mu_p + muLep - mu_n)) / temp);

    // without detailed balance
    //out[0] = etanp * aux * (1 - fd_e); // Neutrino absorptivity [s-1], Eq.(C13)
  }

  // Check kinematics constraint for antineutrino absorption
  if (E_p - mLep > 0.) {
    aux = Rbar * kClight * kAbsEmConst * E_p * E_p * sqrt(1. - mLep * mLep / (E_p * E_p)); // remove c to get output in cm-1
    fd_p = FermiDistr(E_p, temp, -muLep);
    // @TODO: eventually think about a specifically designed function for (1-FermiDistr)
    out[3] = etanp * aux * fd_p;        // Antineutrino emissivity   [s-1], Eq.(C20)
    out[2] = out[3] * SafeExp((omega - (mu_n - mu_p - muLep)) / temp);
    
    // without detailed balance
    //out[2] = etapn * aux * (1. - fd_p); // Antineutrino absorptivity [s-1], Eq.(C19)
  }
}

/* Compute the absortivity and inverse mean free path
 *
 * Outputs:
 *      em_nue = j_nu [s^-1] for nue
 *      ab_nue = 1/lambda_nu [s^-1] for nue
 *      em_anue = j_nu [s^-1] for anue
 *      ab_anue = 1/lambda_nu [s^-1] for anue
 *
 * @TODO: add support for muons
 */
MyOpacity AbsOpacity(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  MyOpacity MyOut = {0.0}; // initialize to zero

  // Electron (anti)neutrino
  double el_out[4] = {0.0};
  AbsOpacitySingleLep(omega, opacity_pars, eos_pars, kMe, eos_pars->mu_e, el_out);

  MyOut.abs[id_nue] = el_out[0];
  MyOut.em[id_nue] = el_out[1];
  MyOut.abs[id_anue] = el_out[2];
  MyOut.em[id_anue] = el_out[3];

  // Uncomment the following when considering also muons 
  // // Muon (anti)neutrino
  //double mu_out[4] = {0.0};
  //AbsOpacitySingleLep(omega, opacity_pars, eos_pars, kMmu, eos_pars->mu_mu, mu_out);

  //MyOut.abs[id_num] = mu_out[0];
  //MyOut.em[id_num] = mu_out[1];
  //MyOut.abs[id_anum] = mu_out[2];
  //MyOut.em[id_anum] = mu_out[3];

  return MyOut;
}

/* Compute the stimulated absorption opacities
 *
 * For all species:
 *      em_nu = j_nu [s^-1]
 *      ab_nu = j_nu + 1/lambda_nu [s^-1]
 *
 * Both opacities are energy dependent
 * @TODO: add support for muons
 */
MyOpacity StimAbsOpacity(const double omega, OpacityParams *opacity_pars, MyEOSParams *eos_pars) {
  MyOpacity abs_opacity = AbsOpacity(omega, opacity_pars, eos_pars);
  MyOpacity stim_abs_opacity = {
      .em[id_nue]  = abs_opacity.em[id_nue],
      .abs[id_nue]  = abs_opacity.abs[id_nue] + abs_opacity.em[id_nue],
      .em[id_anue] = abs_opacity.em[id_anue],
      .abs[id_anue] = abs_opacity.abs[id_anue] + abs_opacity.em[id_anue],
      .abs[id_nux] = 0.,
      .em[id_nux] = 0.
  };

  return stim_abs_opacity;
}
