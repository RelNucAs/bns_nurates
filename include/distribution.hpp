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
BS_REAL NuFThick(const BS_REAL omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{
    const BS_REAL T_t  = distr_pars->temp_t[nuid];
    const BS_REAL mu_t = distr_pars->temp_t[nuid] * distr_pars->eta_t[nuid];

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
    constexpr BS_REAL zero         = 0;
    constexpr BS_REAL one          = 1;
    constexpr BS_REAL half         = 0.5;
    constexpr BS_REAL three_halves = 1.5;
    constexpr BS_REAL one_third    = 1. / 3.;

    constexpr BS_REAL twenty = 20;
    constexpr BS_REAL thirty = 30;

    constexpr BS_REAL y1 = 0.005;
    constexpr BS_REAL y2 = 0.7;
    constexpr BS_REAL y3 = 0.7901234567745267;

    constexpr BS_REAL kBS_ThickDistr_num_1[6] = {
        0.19926987701997, 38865.0149220478, 14364.6331099737,
        5750.1878570758,  1120.71610972194, 1.60356108438235e-8};
    constexpr BS_REAL kBS_ThickDistr_num_2[7] = {
        41.3836568203438,    32.5515666786612, 157.774993512235,
        66.5726772426253,    14.4883415579211, 0.315360380575709,
        0.000660414331285249};
    constexpr BS_REAL kBS_ThickDistr_num_3[7] = {
        3852.81416018959,    5316.18895799799, 1102.91561586553,
        1.54082262710661e-6, 1732.89925128741, 1769.59868329086,
        586.406885304906};

    constexpr BS_REAL kBS_ThickDistr_den_1[6] = {1.0,
                                                 38840.0743174942,
                                                 99.9009680656931,
                                                 171.874844843596,
                                                 75.7101579899442,
                                                 83.0160130941424};
    constexpr BS_REAL kBS_ThickDistr_den_2[6] = {
        1.8888797407042,   5.35488690539183,   1.94673781342617,
        0.483128792035557, 0.0113386564109086, 2.64160073447322e-5};
    constexpr BS_REAL kBS_ThickDistr_den_3[6] = {
        255.936658313629, 9.42360945627147e-5, 81.2467063138386,
        180.100197053091, 89.0343496217014,    143.849128123195};

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        const BS_REAL n = fmax(1e-100, M1_pars->n[nuid]); // [nm^-3]
        // Convert J from g nm-1 s^-1 to MeV nm-3
        const BS_REAL J   = fmax(1e-100, M1_pars->J[nuid]) / kBS_MeV;
        const BS_REAL chi = M1_pars->chi[nuid];

        // BS_ASSERT(
            // n > zero,
            // "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
            // nuid, nuid, n);
        BS_ASSERT(chi >= one_third && chi <= one,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_t[nuid] = three_halves * (one - chi);

        const BS_REAL y = n * POW3(n / J) / kBS_FourPi_hc3;

        if (y < y1)
        {
            out_distribution_pars->eta_t[nuid] =
                log((y * (y * (y * (y * (y * (y + kBS_ThickDistr_num_1[0]) +
                                         kBS_ThickDistr_num_1[1]) +
                                    kBS_ThickDistr_num_1[2]) +
                               kBS_ThickDistr_num_1[3]) +
                          kBS_ThickDistr_num_1[4]) +
                     kBS_ThickDistr_num_1[5]) /
                    (y * (y * (y * (y * (y * (y + kBS_ThickDistr_den_1[0]) +
                                         kBS_ThickDistr_den_1[1]) -
                                    kBS_ThickDistr_den_1[2]) -
                               kBS_ThickDistr_den_1[3]) +
                          kBS_ThickDistr_den_1[4]) +
                     kBS_ThickDistr_den_1[5]));
        }
        else if (y <= y2)
        {
            out_distribution_pars->eta_t[nuid] =
                (y * (y * (y * (y * (y * (kBS_ThickDistr_num_2[0] * y +
                                          kBS_ThickDistr_num_2[1]) -
                                     kBS_ThickDistr_num_2[2]) +
                                kBS_ThickDistr_num_2[3]) +
                           kBS_ThickDistr_num_2[4]) +
                      kBS_ThickDistr_num_2[5]) +
                 kBS_ThickDistr_num_2[6]) /
                    (y * (y * (y * (y * (y * (y + kBS_ThickDistr_den_2[0]) -
                                         kBS_ThickDistr_den_2[1]) +
                                    kBS_ThickDistr_den_2[2]) +
                               kBS_ThickDistr_den_2[3]) +
                          kBS_ThickDistr_den_2[4]) +
                     kBS_ThickDistr_den_2[5]) -
                thirty;
        }
        else if (y > y2 && y < y3)
        {
            out_distribution_pars->eta_t[nuid] =
                exp((y * (y * (y * (y * (y * (kBS_ThickDistr_num_3[0] * y -
                                              kBS_ThickDistr_num_3[1]) +
                                         kBS_ThickDistr_num_3[2]) +
                                    kBS_ThickDistr_num_3[3]) +
                               kBS_ThickDistr_num_3[4]) -
                          kBS_ThickDistr_num_3[5]) +
                     kBS_ThickDistr_num_3[6]) /
                    (y * (y * (y * (y * (y * (y + kBS_ThickDistr_den_3[0]) +
                                         kBS_ThickDistr_den_3[1]) -
                                    kBS_ThickDistr_den_3[2]) -
                               kBS_ThickDistr_den_3[3]) -
                          kBS_ThickDistr_den_3[4]) +
                     kBS_ThickDistr_den_3[5]));
        }
        else
        {
            out_distribution_pars->eta_t[nuid] = twenty;
        }

        out_distribution_pars->eta_t[nuid] =
            fmin(out_distribution_pars->eta_t[nuid], twenty);

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
BS_REAL NuFThin(const BS_REAL omega, const NuDistributionParams* distr_pars,
                const int nuid)
{
    const BS_REAL T_f    = distr_pars->temp_f[nuid];
    const BS_REAL c_f    = distr_pars->c_f[nuid];
    const BS_REAL beta_f = distr_pars->beta_f[nuid];

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
    constexpr BS_REAL zero      = 0;
    constexpr BS_REAL one       = 1;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL half      = 0.5;
    constexpr BS_REAL one_third = 1. / 3.;

    constexpr BS_REAL c_f = CONST_C_F;

    for (int nuid = 0; nuid < total_num_species; ++nuid)
    {
        const BS_REAL n = fmax(1e-100, M1_pars->n[nuid]); // [nm^-3]
        // Convert J from g nm-1 s^-1 to MeV nm-3
        const BS_REAL J   = fmax(1e-100, M1_pars->J[nuid]) / kBS_MeV;
        const BS_REAL chi = M1_pars->chi[nuid];

        // BS_ASSERT(
            // n > zero,
            // "Neutrino (species=%d) number density is non-positive (n[%d]=%e).",
            // nuid, nuid, n);
        BS_ASSERT(chi >= one_third && chi <= one,
                  "Invalid Eddington factor (chi[%d]=%e).", nuid, chi);

        out_distribution_pars->w_f[nuid] = half * (three * chi - one);

        out_distribution_pars->c_f[nuid] = c_f;

        const BS_REAL Tnu                   = J / (n * (c_f + three));
        out_distribution_pars->temp_f[nuid] = Tnu;

        out_distribution_pars->beta_f[nuid] =
            n / (kBS_FourPi_hc3 * GammaStirling(c_f + three) *
                 pow(Tnu, c_f + three));
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
    constexpr BS_REAL zero     = 0;
    constexpr BS_REAL one      = 1;
    constexpr BS_REAL thousand = 1e3;

    constexpr BS_REAL c_f = CONST_C_F;

    NuDistributionParams out;

    const BS_REAL T    = eos_pars->temp; // [MeV]
    const BS_REAL mu_e = eos_pars->mu_e; // [MeV]
    const BS_REAL mu_p = eos_pars->mu_p; // [MeV]
    const BS_REAL mu_n = eos_pars->mu_n; // [MeV]

    BS_ASSERT(T > zero && T < thousand,
              "Given temperature is either negative or beyond 1 GeV (T=%e).",
              T);

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        out.w_t[idx]    = one;
        out.temp_t[idx] = T;

        out.w_f[idx]    = zero;
        out.c_f[idx]    = c_f;
        out.temp_f[idx] = T;
        out.beta_f[idx] = one;
    }

    out.eta_t[id_nue]  = (mu_e - mu_n + mu_p) / T;
    out.eta_t[id_anue] = -out.eta_t[id_nue];
    out.eta_t[id_nux]  = zero;
    out.eta_t[id_anux] = zero;

    return out;
}

/* Total neutrino distribution combining optically thick and thin regimes
 *
 * omega:       neutrino energy
 * distr_pars:  neutrino distribution parameters for thick and thin regimes
 * species:     species of neutrino
 */
KOKKOS_INLINE_FUNCTION
BS_REAL TotalNuF(const BS_REAL omega, const NuDistributionParams* distr_pars,
                 const int nuid)
{
    constexpr BS_REAL zero = 0;

    BS_ASSERT(omega >= zero, "Neutrino energy is negative.");
    BS_ASSERT(nuid >= 0 && nuid < total_num_species,
              "Invalid neutrino species ID.");

    const BS_REAL w_t = distr_pars->w_t[nuid];
    const BS_REAL w_f = distr_pars->w_f[nuid];

    const BS_REAL f_thick = NuFThick(omega, distr_pars, nuid);
    const BS_REAL f_thin  = NuFThin(omega, distr_pars, nuid);

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
MyQuadratureIntegrand NuNumberIntegrand(BS_REAL* x, void* p)
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

    BS_REAL s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuNumberIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const BS_REAL result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

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
inline MyQuadratureIntegrand NuEnergyIntegrand(BS_REAL* x, void* p)
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

    BS_REAL s[total_num_species];

    s[id_nue]  = fabs(distr_pars->temp_t[id_nue] * distr_pars->eta_t[id_nue]);
    s[id_anue] = fabs(distr_pars->temp_t[id_anue] * distr_pars->eta_t[id_anue]);
    s[id_nux]  = s[id_nue]; // @TODO: cannot be equal to zero
    s[id_anux] = s[id_nue]; // @TODO: cannot be equal to zero

    integrand.function = &NuEnergyIntegrand;
    MyQuadratureIntegrand result =
        GaussLegendreIntegrate1D(&quad, &integrand, s);
    const BS_REAL result_factor = 4. * kBS_Pi / POW3(kBS_H * kBS_Clight);

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
    const BS_REAL n_prefactor = kBS_FourPi_hc3 * POW3(eos_pars->temp);
    const BS_REAL j_prefactor = n_prefactor * eos_pars->temp;

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
