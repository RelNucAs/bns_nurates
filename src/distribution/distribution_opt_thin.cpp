// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_opt_thin.c
//  \brief compute neutrino distribution function & M1 parameters for optically
//  thin regime

#include <math.h>

#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))

#define CONST_C_F 0.6


/* Neutrino distribution function for optically thin regime: Maxwell-Boltzmann
 * distribution
 *
 * omega:       neutrino energy
 * distr_pars:  optically thin parameters
 * species:     species of neutrino
 */
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
