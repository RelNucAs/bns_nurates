// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_total.c
//  \brief functions for constructing the total neutrino distribution function by combining the optically thick and thin regimes

#include <math.h>

#include "distribution.h"
#include "../constants.h"
#include "../integration/integration.h"

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 * species:     species of neutrino
 */
double TotalNuF(double omega, NuDistributionParams *distr_pars, int species) {
  double w_t = distr_pars->w_t[species];
  double w_f = distr_pars->w_f[species];

  double f_thick = NuFThick(omega, distr_pars, species);
  double f_thin = NuFThin(omega, distr_pars, species);

  return w_t * f_thick + w_f * f_thin;
}

/* Calculate distribution function parameters in the thick and thin regime from M1 quantities
 *
 * M1_params:   M1 quantities
 * eos_params:  parameters from EOS
 */
NuDistributionParams CalculateDistrParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars) {
  NuDistributionParams out_distr_pars;

  CalculateThickParamsFromM1(M1_pars, eos_pars, &out_distr_pars);
  CalculateThinParamsFromM1(M1_pars, &out_distr_pars);

  return out_distr_pars;
}

/* Integrand for computing neutrino number density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuNumberIntegrand(double *x, void *p) {
  NuDistributionParams *distr_pars = (NuDistributionParams *) p;

  MyQuadratureIntegrand result;

  result.integrand[0] = x[0] * x[0] * TotalNuF(x[0], distr_pars, 0);
  result.integrand[1] = x[0] * x[0] * TotalNuF(x[0], distr_pars, 1);
  result.integrand[2] = x[0] * x[0] * TotalNuF(x[0], distr_pars, 2);

  return result;
}

/* Compute neutrino number density (besides 4*pi factor)
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuNumber(NuDistributionParams *distr_pars) {
  MyFunctionMultiD integrand;

  integrand.dim = 1;
  integrand.params = distr_pars;
  integrand.my_quadrature_integrand.n = 3;

  MyQuadrature quad = quadrature_default;

  GaussLegendreMultiD(&quad);

  double s = distr_pars->temp_t[0] * distr_pars->eta_t[0]; // @TODO: currently breaking the integral at the same point for all species

  integrand.function = &NuNumberIntegrand;
  MyQuadratureIntegrand result = GaussLegendreIntegrate1D(&quad, &integrand, s);
  const double result_factor = 4. * kPi / pow(kH * kClight, 3.);

  for (int species = 0; species < 3; species++) {
    result.integrand[species] = result.integrand[species] * result_factor;
  }

  return result;
}


/* Integrand for computing the neutrino energy density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuEnergyIntegrand(double *x, void *p) {
  MyQuadratureIntegrand result = NuNumberIntegrand(x, p);

  for (int species = 0; species < 3; species++) {
    result.integrand[species] = x[0] * result.integrand[species];
  }

  return result;
}

/* Compute neutrino energy density (besides 4*pi factor)
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuEnergy(NuDistributionParams *distr_pars) {
  MyFunctionMultiD integrand;

  integrand.dim = 1;
  integrand.params = distr_pars;
  integrand.my_quadrature_integrand.n = 3;

  MyQuadrature quad = quadrature_default;

  GaussLegendreMultiD(&quad);

  double s = distr_pars->temp_t[0] * distr_pars->eta_t[0]; // @TODO: currently breaking the integral at the same point for all species

  integrand.function = &NuEnergyIntegrand;
  MyQuadratureIntegrand result = GaussLegendreIntegrate1D(&quad, &integrand, s);
  const double result_factor = 4. * kPi / pow(kH * kClight, 3.);

  for (int species = 0; species < 3; species++) {
    result.integrand[species] = result.integrand[species] * result_factor;
  }

  return result;
}
