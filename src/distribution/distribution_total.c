// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_total.c
//  \brief functions for constructing the total neutrino distribution function
//  by combining the optically thick and thin regimes

#include <math.h>

#include "distribution.h"
#include "constants.h"
#include "integration.h"
#include "functions.h"

#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))


/* Function for evaluating parameters of neutrino distribution function at
 * equilibrium
 */
NuDistributionParams NuEquilibriumParams(MyEOSParams* eos_pars)
{
    NuDistributionParams out_distr;

    const double temp = eos_pars->temp; // [MeV]
    const double mu_e = eos_pars->mu_e; // [MeV]
    const double mu_p = eos_pars->mu_p; // [MeV]
    const double mu_n = eos_pars->mu_n; // [MeV]

    for (int idx = 0; idx < total_num_species; idx++)
    {
        out_distr.w_t[idx]    = 1.;
        out_distr.temp_t[idx] = temp;

        out_distr.w_f[idx]    = 0.;
        out_distr.c_f[idx]    = 1.; //
        out_distr.temp_f[idx] = temp;
        out_distr.beta_f[idx] = 1.;
    }

    out_distr.eta_t[id_nue]  = (mu_e - mu_n + mu_p) / temp;
    out_distr.eta_t[id_anue] = -out_distr.eta_t[id_nue];
    out_distr.eta_t[id_nux]  = 0.;
    out_distr.eta_t[id_anux] = 0.;

    return out_distr;
}

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 * species:     species of neutrino
 */
double TotalNuF(double omega, NuDistributionParams* distr_pars, int species)
{
    double w_t = distr_pars->w_t[species];
    double w_f = distr_pars->w_f[species];

    double f_thick = NuFThick(omega, distr_pars, species);
    double f_thin  = NuFThin(omega, distr_pars, species);

    return w_t * f_thick + w_f * f_thin;
}

/* Calculate distribution function parameters in the thick and thin regime from
 * M1 quantities
 *
 * M1_params:   M1 quantities
 * eos_params:  parameters from EOS
 */
NuDistributionParams CalculateDistrParamsFromM1(M1Quantities* M1_pars,
                                                MyEOSParams* eos_pars)
{
    NuDistributionParams out_distr_pars;

    CalculateThickParamsFromM1(M1_pars, eos_pars, &out_distr_pars);
    CalculateThinParamsFromM1(M1_pars, &out_distr_pars);

    return out_distr_pars;
}

/* Integrand for computing neutrino number density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuNumberIntegrand(double* x, void* p)
{
    NuDistributionParams* distr_pars = (NuDistributionParams*)p;

    MyQuadratureIntegrand result;

    result.n = total_num_species;
    for (int idx = 0; idx < total_num_species; idx++)
    {
        result.integrand[idx] = x[0] * x[0] * TotalNuF(x[0], distr_pars, idx);
    }

    return result;
}

/* Compute neutrino number density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuNumber(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    double s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuNumberIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const double result_factor = 4. * kPi / POW3(kHClight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}

/* Integrand for computing the neutrino energy density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuEnergyIntegrand(double* x, void* p)
{
    MyQuadratureIntegrand result = NuNumberIntegrand(x, p);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = x[0] * result.integrand[species];
    }

    return result;
}

/* Compute neutrino energy density
 *
 * Computes this for three neutrino species
 */
MyQuadratureIntegrand NuEnergy(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    double s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuEnergyIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const double result_factor = 4. * kPi / POW3(kHClight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}

void ComputeM1DensitiesEq(MyEOSParams* eos_params,
                          NuDistributionParams* nu_distribution_params,
                          M1Quantities* m1_pars)
{

    const double four_pi_hc3 = 4. * kPi / POW3(kHClight);
    const double n_prefactor = four_pi_hc3 * POW3(eos_params->temp);
    const double j_prefactor = n_prefactor * eos_params->temp;

    for (int idx = 0; idx < total_num_species; idx++)
    {
        m1_pars->n[idx] =
            n_prefactor * FDI_p2(nu_distribution_params->eta_t[idx]);
        m1_pars->J[idx] =
            j_prefactor * FDI_p3(nu_distribution_params->eta_t[idx]);
    }

    return;
}
