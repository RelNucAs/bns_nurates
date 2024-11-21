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

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"

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
bs_real NuFThick(const bs_real omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{
    const bs_real T_t  = distr_pars->temp_t[nuid];
    const bs_real mu_t = distr_pars->temp_t[nuid] * distr_pars->eta_t[nuid];

    return FermiDistr(omega, T_t, mu_t);
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
void CalculateThickParamsFromM1(const M1Quantities* M1_pars,
                                NuDistributionParams* out_distribution_pars)
{

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        const bs_real n = M1_pars->n[nuid]; // [nm^-3]
        // Convert J from g nm-1 s^-1 to MeV nm-3
        const bs_real J   = M1_pars->J[nuid] / kBS_MeV;
        const bs_real chi = M1_pars->chi[nuid];

        BS_ASSERT(
            n > 0.,
            "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
            nuid, nuid, n);
        BS_ASSERT(chi >= 1. / 3. && chi <= 1.,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_t[nuid] = 1.5 * (1. - chi);

        const bs_real y = n * POW3(n / J) / kBS_FourPi_hc3;

        if (y < 0.005)
        {
            out_distribution_pars->eta_t[nuid] =
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
            out_distribution_pars->eta_t[nuid] =
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
            out_distribution_pars->eta_t[nuid] =
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
            out_distribution_pars->eta_t[nuid] = 20.;
        }

        out_distribution_pars->eta_t[nuid] =
            fmin(out_distribution_pars->eta_t[nuid], 20.);

        out_distribution_pars->temp_t[nuid] =
            FDI_p2(out_distribution_pars->eta_t[nuid]) * J /
            (FDI_p3(out_distribution_pars->eta_t[nuid]) * n);
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
bs_real NuFThin(const bs_real omega, const NuDistributionParams* distr_pars,
                const int nuid)
{
    const bs_real T_f    = distr_pars->temp_f[nuid];
    const bs_real c_f    = distr_pars->c_f[nuid];
    const bs_real beta_f = distr_pars->beta_f[nuid];

    return beta_f * pow(omega, c_f) * exp(-omega / T_f);
}

/* Recover distribution function parameters for optically thin regime from M1
 * quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy
 * density) out_distr_pars:  computes free neutrino temperature and c_f
 */
KOKKOS_INLINE_FUNCTION
void CalculateThinParamsFromM1(const M1Quantities* M1_pars,
                               NuDistributionParams* out_distribution_pars)
{
    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        const bs_real n = M1_pars->n[nuid]; // [nm^-3]
        // Convert J from g nm-1 s^-1 to MeV nm-3
        const bs_real J   = M1_pars->J[nuid] / kBS_MeV;
        const bs_real chi = M1_pars->chi[nuid];

        BS_ASSERT(
            n > 0.,
            "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
            nuid, nuid, n);
        BS_ASSERT(chi >= 1. / 3. && chi <= 1.,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_f[nuid] = 0.5 * (3. * chi - 1.);

        out_distribution_pars->c_f[nuid] = CONST_C_F;

        const bs_real Tnu                   = J / (n * (CONST_C_F + 3.));
        out_distribution_pars->temp_f[nuid] = Tnu;

        out_distribution_pars->beta_f[nuid] =
            n / (kBS_FourPi_hc3 * GammaStirling(CONST_C_F + 3.) *
                 pow(Tnu, CONST_C_F + 3.));
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
NuDistributionParams NuEquilibriumParams(const MyEOSParams* eos_pars)
{
    NuDistributionParams out;

    const bs_real T    = eos_pars->temp; // [MeV]
    const bs_real mu_e = eos_pars->mu_e; // [MeV]
    const bs_real mu_p = eos_pars->mu_p; // [MeV]
    const bs_real mu_n = eos_pars->mu_n; // [MeV]

    BS_ASSERT(T > 0. && T < 1e3,
              "Given temperature is either negative or beyond 1 GeV.");

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        out.w_t[idx]    = 1.;
        out.temp_t[idx] = T;

        out.w_f[idx]    = 0.;
        out.c_f[idx]    = CONST_C_F;
        out.temp_f[idx] = T;
        out.beta_f[idx] = 1.;
    }

    out.eta_t[id_nue]  = (mu_e - mu_n + mu_p) / T;
    out.eta_t[id_anue] = -out.eta_t[id_nue];
    out.eta_t[id_nux]  = 0.;
    out.eta_t[id_anux] = 0.;

    return out;
}

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 * species:     species of neutrino
 */
KOKKOS_INLINE_FUNCTION
bs_real TotalNuF(const bs_real omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{

    BS_ASSERT(omega >= 0., "Neutrino energy is negative.");
    BS_ASSERT(nuid >= 0 && nuid < total_num_species,
              "Invalid neutrino species ID.");

    const bs_real w_t = distr_pars->w_t[nuid];
    const bs_real w_f = distr_pars->w_f[nuid];

    const bs_real f_thick = NuFThick(omega, distr_pars, nuid);
    const bs_real f_thin  = NuFThin(omega, distr_pars, nuid);

    return w_t * f_thick + w_f * f_thin;
}

/* Calculate distribution function parameters in the thick and thin regime from
 * M1 quantities
 *
 * M1_params:   M1 quantities
 * eos_params:  parameters from EOS
 */
KOKKOS_INLINE_FUNCTION
NuDistributionParams CalculateDistrParamsFromM1(const M1Quantities* M1_pars,
                                                const MyEOSParams* eos_pars)
{
    NuDistributionParams out;

    CalculateThickParamsFromM1(M1_pars, &out);
    CalculateThinParamsFromM1(M1_pars, &out);

    return out;
}

/* Integrand for computing neutrino number density
 *
 * Computes this for three neutrino species
 */
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand NuNumberIntegrand(bs_real* x, void* p)
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
/*
inline MyQuadratureIntegrand NuNumber(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    bs_real s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuNumberIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const bs_real result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}
*/

/* Integrand for computing the neutrino energy density
 *
 * Computes this for three neutrino species
 */
inline MyQuadratureIntegrand NuEnergyIntegrand(bs_real* x, void* p)
{
    MyQuadratureIntegrand result = NuNumberIntegrand(x, p);

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        result.integrand[nuid] = x[0] * result.integrand[nuid];
    }

    return result;
}

/* Compute neutrino energy density
 *
 * Computes this for three neutrino species
 */
/*
inline MyQuadratureIntegrand NuEnergy(NuDistributionParams* distr_pars)
{
    MyFunctionMultiD integrand;

    integrand.dim                       = 1;
    integrand.params                    = distr_pars;
    integrand.my_quadrature_integrand.n = total_num_species;

    MyQuadrature quad = quadrature_default;

    GaussLegendreMultiD(&quad);

    bs_real s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuEnergyIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const bs_real result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

    for (int species = 0; species < total_num_species; ++species)
    {
        result.integrand[species] = result.integrand[species] * result_factor;
    }

    return result;
}
*/

KOKKOS_INLINE_FUNCTION
void ComputeM1DensitiesEq(const MyEOSParams* eos_pars,
                          const NuDistributionParams* nu_distribution_params,
                          M1Quantities* m1_pars)
{
    const bs_real n_prefactor = kBS_FourPi_hc3 * POW3(eos_pars->temp);
    const bs_real j_prefactor = n_prefactor * eos_pars->temp;

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        m1_pars->n[nuid] =
            n_prefactor * FDI_p2(nu_distribution_params->eta_t[nuid]);
        m1_pars->J[nuid] =
            j_prefactor * FDI_p3(nu_distribution_params->eta_t[nuid]);

        BS_ASSERT(isfinite(m1_pars->n[nuid]),
                  "Neutrino (species=%d) number density computed assuming "
                  "equilibrium is not finite (n=%e nm^-3).",
                  nuid, m1_pars->n[nuid]);
        BS_ASSERT(isfinite(m1_pars->J[nuid]),
                  "Neutrino (species=%d) energy density computed assuming "
                  "equilibrium is not finite (J=%e MeV nm^-3).",
                  nuid, m1_pars->J[nuid]);
    }

    return;
}

#endif // BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
