// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution.h
//  \brief header file for distribution function reconstruction from M1
//  parameters
//         supports three different species: electron neutrino, electron
//         anti-neutrino and mu/tau neutrino/antineutrino

#ifndef BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
#define BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_

#include <math.h>
#include <Kokkos_Core.hpp>

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"
#include "integration.hpp"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))
#define CONST_C_F 0.6

/* ===========================================================================
 * Functions for the optically thick distribution function
 * ===========================================================================
 */

/* Neutrino distribution function in optically thick regime: Fermi-Dirac
distribution
 *
 * omega:       neutrino energy
 * distr_pars:  uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 * species:     species of neutrino
 */
KOKKOS_INLINE_FUNCTION
double NuFThick(double omega, NuDistributionParams* distr_pars, int species)
{
    const double temp = distr_pars->temp_t[species];
    const double mu = distr_pars->temp_t[species] * distr_pars->eta_t[species];
    return FermiDistr(omega, temp, mu);
}

/* Recover distribution function parameters for optically thick regime from M1
quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy
density)
 * eos_pars:        uses fluid temperature
 * out_distr_pars:  computes trapped neutrino temperature and degeneracy
parameter
 */
KOKKOS_INLINE_FUNCTION
void CalculateThickParamsFromM1(M1Quantities* M1_pars, MyEOSParams* eos_pars,
                                NuDistributionParams* out_distr_pars)
{
    (void)eos_pars;

    for (int species = 0; species < total_num_species; ++species)
    {
        // n and J are assumed to be in cgs units
        const double n = M1_pars->n[species];
        const double j =
            M1_pars->J[species] / kMeV; // conversion from erg cm-3 to MeV cm-3
        const double chi = M1_pars->chi[species];

        if (n == 0.)
        {
            out_distr_pars->w_t[species]    = 0.;
            out_distr_pars->eta_t[species]  = 1.;
            out_distr_pars->temp_t[species] = 1.;

            continue;
        }

        out_distr_pars->w_t[species] = 1.5 * (1. - chi);

        const double y = POW4(n) * POW3(kHClight) / (4. * kPi * POW3(j));

        if (y < 0.005)
        {
            out_distr_pars->eta_t[species] =
                log((y * (y * (y * (y * (y * (y + 0.19926987701997) +
                                         38865.0149220478) +
                                    14364.6331099737) +
                               5750.1878570758) +
                          1120.71610972194) +
                     1.60356108438235e-8) /
                    (y * (y * (y * (y * (y * (y + 1.0) + 38840.0743174942) -
                                    99.9009680656931) -
                               171.874844843596) +
                          75.7101579899442) +
                     83.0160130941424));
        }
        else if (y <= 0.7)
        {
            out_distr_pars->eta_t[species] =
                (y * (y * (y * (y * (y * (41.3836568203438 * y +
                                          32.5515666786612) -
                                     157.774993512235) +
                                66.5726772426253) +
                           14.4883415579211) +
                      0.315360380575709) +
                 0.000660414331285249) /
                    (y * (y * (y * (y * (y * (y + 1.8888797407042) -
                                         5.35488690539183) +
                                    1.94673781342617) +
                               0.483128792035557) +
                          0.0113386564109086) +
                     2.64160073447322e-5) -
                30.;
        }
        else if (y > 0.7 && y < 0.7901234567745267)
        {
            out_distr_pars->eta_t[species] =
                exp((y * (y * (y * (y * (y * (3852.81416018959 * y -
                                              5316.18895799799) +
                                         1102.91561586553) +
                                    1.54082262710661e-6) +
                               1732.89925128741) -
                          1769.59868329086) +
                     586.406885304906) /
                    (y * (y * (y * (y * (y * (y + 255.936658313629) +
                                         9.42360945627147e-5) -
                                    81.2467063138386) -
                               180.100197053091) -
                          89.0343496217014) +
                     143.849128123195));
        }
        else
        {
            out_distr_pars->eta_t[species] = 20.;
        }

        out_distr_pars->eta_t[species] =
            MIN(out_distr_pars->eta_t[species], 20.);

        out_distr_pars->temp_t[species] =
            FDI_p2(out_distr_pars->eta_t[species]) * j /
            (FDI_p3(out_distr_pars->eta_t[species]) * n);
    }
}

/* ===========================================================================
 * Functions for the optically thin distribution function
 * ===========================================================================
 */

/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann
 * distribution
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 * species:     species of neutrino
 */
KOKKOS_INLINE_FUNCTION
double NuFThin(double omega, NuDistributionParams* distr_pars, int species)
{
    const double temp_f = distr_pars->temp_f[species];
    const double c_f    = distr_pars->c_f[species];
    const double beta_f = distr_pars->beta_f[species];

    return beta_f * pow(omega, c_f) * exp(-omega / temp_f);
}

/* Recover distribution function parameters for optically thin regime from M1
 * quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy
 * density) out_distr_pars:  computes free neutrino temperature and c_f
 */
KOKKOS_INLINE_FUNCTION
void CalculateThinParamsFromM1(M1Quantities* M1_pars,
                               NuDistributionParams* out_distr_pars)
{
    for (int species = 0; species < total_num_species; species++)
    {
        const double n = M1_pars->n[species];
        const double j =
            M1_pars->J[species] / kMeV; // conversion from erg cm-3 to MeV cm-3
        const double chi = M1_pars->chi[species];

        if (n == 0.)
        {
            out_distr_pars->w_f[species]    = 0.;
            out_distr_pars->beta_f[species] = 0.;
            out_distr_pars->c_f[species]    = CONST_C_F;
            out_distr_pars->temp_f[species] = 1.;

            continue;
        }

        out_distr_pars->w_f[species] = 0.5 * (3. * chi - 1.);

        out_distr_pars->c_f[species] = CONST_C_F;

        out_distr_pars->temp_f[species] = j / (n * (CONST_C_F + 3.));

        // n and J are assumed to be in cgs units
        out_distr_pars->beta_f[species] =
            n * POW3(kHClight) /
            (4. * kPi * GammaStirling(CONST_C_F + 3.) *
             pow(out_distr_pars->temp_f[species], CONST_C_F + 3.));
    }
}

/* ===========================================================================
 * Functions for constructing the total distribution function
 * ===========================================================================
 */

/* Function for evaluating parameters of neutrino distribution function at
 * equilibrium
 */
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
inline MyQuadratureIntegrand NuNumber(NuDistributionParams* distr_pars)
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
inline MyQuadratureIntegrand NuEnergyIntegrand(double* x, void* p)
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
inline MyQuadratureIntegrand NuEnergy(NuDistributionParams* distr_pars)
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

KOKKOS_INLINE_FUNCTION
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

#endif // BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
