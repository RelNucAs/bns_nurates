// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file distribution_opt_thick.h
//  \brief compute neutrino distribution function & M1 parameters for optically thick regime

#include <math.h>

#include "distribution.h"
#include "bns_nurates.h"
#include "constants.h"
#include "functions.h"

/* Neutrino distribution function in optically thick regime: Fermi-Dirac distribution
 *
 * omega:       neutrino energy
 * distr_pars:  uses temp_t (fluid temperature) and eta_t (degeneracy parameter)
 * species:     species of neutrino
 */
double NuFThick(double omega, NuDistributionParams *distr_pars, int species) {
  double temp = distr_pars->temp_t[species];
  double mu = distr_pars->temp_t[species] * distr_pars->eta_t[species]; // @TODO: simplify if needed
  return FermiDistr(omega, temp, mu);
}

/* Recover distribution function parameters for optically thick regime from M1 quantities
 *
 * M1_params:       uses n (neutrino number density) and J (neutrino energy density)
 * eos_pars:        uses fluid temperature
 * out_distr_pars:  computes trapped neutrino temperature and degeneracy parameter
 */


void CalculateThickParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars, NuDistributionParams *out_distr_pars) {
  (void)eos_pars;

  // set degeneracy parameter for different neutrino species

  //out_distr_pars->eta_t[0] = (eos_pars->mu_p + eos_pars->mu_e - eos_pars->mu_n) / eos_pars->temp;
  //out_distr_pars->eta_t[1] = -(eos_pars->mu_p + eos_pars->mu_e - eos_pars->mu_n) / eos_pars->temp;
  //out_distr_pars->eta_t[2] = 0.;

  for (int species = 0; species < total_num_species; species++) {

    out_distr_pars->w_t[species] = 1.5 * (1. - M1_pars->chi[species]);
    //out_distr_pars->temp_t[species] = eos_pars->temp;

    // @TODO: Currently disabling alternative. What do we do with this ?

    double y = M1_pars->n[species] * M1_pars->n[species] * M1_pars->n[species] * M1_pars->n[species] * kHClight * kHClight * kHClight
        / (4. * M_PI * M1_pars->J[species] * M1_pars->J[species] * M1_pars->J[species]);

    if (y < 0.04) {
      out_distr_pars->eta_t[species] =
          (y * (y * (y * (y * (193601090.674965 - 1108185464.38267 * y) - 5417182.68352186) - 132141.667235385) - 103.882788946162) - 0.0014855034293057) /
              (y * (y * (y * (1.0 * y + 5426593.58422611) + 35266.6610309076) + 15.7105593937989) + 7.27351598778796e-5);
    } else if (y <= 0.7) {
      out_distr_pars->eta_t[species] =
          (y * (y * (y * (y * (2.98735454268239 * y - 4.88726196424146) - 0.901791627658626) + 2.4819363543555) - 0.00548761969223768) - 0.00769760353422157) /
              (y * (y * (y * (1.0 * y - 1.58869835679685) + 0.441658509378521) + 0.14883216362587) + 0.0024946879233069);
    } else if (y > 0.7 && y < 0.7901234567745267) {
      out_distr_pars->eta_t[species] =
          (y * (y * (y * (y * (2827.84724452959 * y - 3707.82829755322) - 2.49562579572847) + 1254.79286964601) - 258.239389707494) - 3.58793204563292) /
              (y * (y * (y * (1.0 * y + 85.5940548187225) - 68.8227897049312) - 55.1337376497705) + 43.9181682873458);
    } else {
      out_distr_pars->eta_t[species] = 1000.;
    }

    out_distr_pars->temp_t[species] = FDI_p2(out_distr_pars->eta_t[species]) * M1_pars->J[species] / (FDI_p3(out_distr_pars->eta_t[species]) * M1_pars->n[species]);
  }
}

