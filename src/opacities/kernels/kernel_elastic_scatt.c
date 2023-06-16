// \file nu_elastic_scatt.c
// \brief Computation of first two Legendre coefficients for neutrino elastic
//        scattering on nucleons
//        Ref: Bruenn, 1985 (https://articles.adsabs.harvard.edu/pdf/1985ApJS...58..771B)

// Possible inclusion of phase space, recoil, weak magnetism correction as in
// Horowitz, 2002 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)

// The elastic scattering kernel is equal for production and absorption and is the same for all neutrino types
// This reaction only contributes for non-isotropic neutrino distribution functions!
// To get the corresponding emissivity/opacity, it requires only a 1D integration over neutrino polar angle


#include <math.h>

#include "kernels.h"
#include "../../constants.h"
#include "../weak_magnetism/weak_magnetism.h"

// @TODO: change names of variables and functions following Google C style

// Definition of constants
const double c0 = (2. * kPi *  kGf * kGf) / kHbar; // constant in (C36) [MeV cm^6 s-1], 
                                                   // Gf^2/hbar = 1.2E-65 MeV cm^6 s^-1 from (C40)
const double c1 = 2. * kPi * c0;  // angular integration over dpi      
const double c2 = 3.*kPi*kPi;                       // constant in Fermi energy calculation
const double c3 = 0.5 * kHbar * kHbar / kMb * kMeV; // constant in Fermi energy calculation

const double h0_p = c1 * (kHpv*kHpv + 3.*kHpa*kHpa); // 0th Leg coeff, protons
const double h1_p = c1 * (kHpv*kHpv -    kHpa*kHpa); // 1st Leg coeff, protons
const double h0_n = c1 * (kHnv*kHnv + 3.*kHna*kHna); // 0th Leg coeff, neutrons
const double h1_n = c1 * (kHnv*kHnv -    kHna*kHna); // 1st Leg coeff, neutrons

/* Inputs:
 *      omega    [MeV] : (anti)neutrino energy (elastic process)
 *      mu             : cosine of the polar angle of nu
 *      mu_prime       : cosine of the polar angle of nu_prime
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
  const double eFN = c3 * pow(c2*nN,2./3.); // Fermi energy computation [MeV]
  const double tmp = 1.5*temp/eFN;
  return nN*tmp/sqrt(1.+tmp*tmp); // [cm-3]
  
  /* Alternative computation: evaluate via numerical derivative (some quick tests showed it's unstable) */
  //dmudrho = ....; //[MeV*cm^3/g]
  //etann = yn / (T*mu*dmudrho); //[1/cm^3]
}

// Legendre expansion of scattering kernel up to l=1
MyKernel nu_N_scatt_kern(ElasticScattParams *kernel_pars,
                         MyEOSParams *eos_pars,
                         const double yN, const int reacflag) {
  double R0 = 1., R1 = 1.;
  double tmp, ker;
	
  const double omega    = kernel_pars->omega;    // (Anti)neutrino energy [MeV]
  const double mu       = kernel_pars->mu;       // cosine of polar angle for nu
  const double mu_prime = kernel_pars->mu_prime; // cosine of polar angle for nu'
  
  const double nb   = eos_pars->nb;   // Number baryon density [cm-3]
  const double temp = eos_pars->temp; // Temperature [MeV]

  const double etaNN = eta_NN_sc(nb, temp, yN); // degeneracy parameter eta_NN, Eq.(C37)

  // Phase space, recoil and weak magnetism corrections
  // R0 (R1) is the correction to the zeroth (first) Legendre coefficient
  if (kernel_pars->use_WM_sc)  WM_scatt(omega, &R0, &R1, reacflag);

  if (reacflag == 1) {
    // Scattering on proton
    tmp = h0_p * R0 + h1_p * R1 * mu * mu_prime; // [MeV cm^3 s-1]
  } else if (reacflag == 2) {
    // Scattering on neutron
    tmp = h0_n * R0 + h1_n * R1 * mu * mu_prime; // [MeV cm^3 s-1]
  }
  
  ker = omega * omega * tmp; // scattering kernel * enu**2 [MeV^3 cm^3 s-1], Eq.(C36)

  // Reaction rate is the same for all neutrino types
  MyKernel elastic_kernel =
      {.absorption_e = ker, .production_e = ker, .absorption_x = ker, .production_x = ker};

  return elastic_kernel;

}


// Scattering kernel for neutrino-proton scattering
MyKernel nu_p_scatt_kern(ElasticScattParams *kernel_params, MyEOSParams *eos_pars) {
  return nu_N_scatt_kern(kernel_params, eos_pars, eos_pars->yp, 1);
}

// Scattering kernel for neutrino-neutron scattering
MyKernel nu_n_scatt_kern(ElasticScattParams *kernel_params, MyEOSParams *eos_pars) {
  return nu_N_scatt_kern(kernel_params, eos_pars, eos_pars->yn, 2);
}

// Sum of scattering kernels for neutrino-nucleon scattering (protons + neutrons)
MyKernel nu_N_scatt_kern_tot(ElasticScattParams *kernel_params, MyEOSParams *eos_pars) {
  MyKernel nu_p_ker = nu_p_scatt_kern(kernel_params, eos_pars); // proton contribution
  MyKernel nu_n_ker = nu_n_scatt_kern(kernel_params, eos_pars); // neutron contribution

  MyKernel elastic_kernel = {.absorption_e = nu_p_ker.absorption_e + nu_p_ker.absorption_e,
                             .production_e = nu_p_ker.production_e + nu_p_ker.production_e,
                             .absorption_x = nu_p_ker.absorption_x + nu_p_ker.absorption_x,
                             .production_x = nu_p_ker.production_x + nu_p_ker.production_x};

  return elastic_kernel;
}
