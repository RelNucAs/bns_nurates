// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_total.h
//  \brief functions for constructing the total neutrino distribution function by combining the optically thick and thin regimes

#include "distribution.h"
#include "../integration/integration.h"

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 */
double TotalNuF(double omega, NuDistributionParams *distr_pars) {
  double w_t = distr_pars->w_t;
  double w_f = distr_pars->w_f;

  double f_thick = NuFThick(omega, distr_pars);
  double f_thin = NuFThin(omega, distr_pars);

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

// Integrand for computation of neutrino number density
double NuNumberIntegrand(double *x, void *p) {
  NuDistributionParams *distr_pars = (NuDistributionParams *) p;

  return x[0] * x[0] * TotalNuF(x[0], distr_pars);
}

// Neutrino number density (besides 4*pi factor)
double NuNumber(NuDistributionParams *distr_pars) {
  MyFunction integrand;

  integrand.dim = 1;
  integrand.params = distr_pars;

  // @TODO: compute quadrature points once for all
  const int n = 32;

  MyQuadrature quad = quadrature_default;

  quad.nx = n;
  quad.type = kGauleg;
  quad.dim = 1;
  quad.x1 = 0.;
  quad.x2 = 1.;

  GaussLegendreMultiD(&quad);

  double s = distr_pars->temp_t * distr_pars->eta_t;

  integrand.function = &NuNumberIntegrand;
  return GaussLegendreIntegrateZeroInf(&quad, &integrand, s);
}

// Integrand for computation of neutrino energy density
double NuEnergyIntegrand(double *x, void *p) {
  return x[0] * NuNumberIntegrand(x, p);
}

// Neutrino energy density (besides 4*pi factor)
double NuEnergy(NuDistributionParams *distr_pars) {
  MyFunction integrand;

  integrand.dim = 1;
  integrand.params = distr_pars;

  // @TODO: compute quadrature points once for all
  const int n = 32;

  MyQuadrature quad = quadrature_default;

  quad.nx = n;
  quad.type = kGauleg;
  quad.dim = 1;
  quad.x1 = 0.;
  quad.x2 = 1.;

  GaussLegendreMultiD(&quad);

  double s = distr_pars->temp_t * distr_pars->eta_t;

  integrand.function = &NuEnergyIntegrand;
  return GaussLegendreIntegrateZeroInf(&quad, &integrand, s);
}
