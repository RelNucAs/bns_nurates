// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================

//! \file distribution_opt_thick.h
//  \brief compute neutrino distribution function & M1 parameters for optically
//  thick regime

#include <math.h>

#include "distribution.h"
#include "bns_nurates.h"
#include "constants.h"
#include "functions.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))


/* Neutrino distribution function in optically thick regime: Fermi-Dirac
distribution
 *
 * omega:       neutrino energy
 * distr_pars:  uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 * species:     species of neutrino
 */
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
void CalculateThickParamsFromM1(M1Quantities* M1_pars, MyEOSParams* eos_pars,
                                NuDistributionParams* out_distr_pars)
{
    (void)eos_pars;

    for (int species = 0; species < total_num_species; ++species)
    {
        // n and J are assumed to be in cgs units
        const double n   = M1_pars->n[species];
        const double j   = M1_pars->J[species] / kMeV; // conversion from erg cm-3 to MeV cm-3
        const double chi = M1_pars->chi[species];

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
