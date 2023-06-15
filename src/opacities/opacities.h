//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_OPACITIES_H_

#include "../bns_nurates.h"

/*===========================================================================*/

// nu_abs_em.h

// Neutrino absorption on neutron (nul + n -> l- + p)
MyOpacity nu_n_abs(const double omega,
                   const double nb, const double temp,
                   const double lep_mass,
                   const double yp, const double yn,
                   const double mu_l, const double mu_hat,
                   const double deltaU);

// Antineutrino absorption on neutron (anul + p -> l+ + n)
MyOpacity nu_p_abs(const double omega,
                   const double nb, const double temp,
                   const double lep_mass,
                   const double yp, const double yn,
                   const double mu_l, const double mu_hat,
                   const double deltaU);


/* Stimulated absoption versions */
// Neutrino absorption on neutron (nul + n -> l- + p)
MyOpacity nu_n_abs_stim(const double omega,
                        const double nb, const double temp,
                        const double lep_mass,
                        const double yp, const double yn,
                        const double mu_l, const double mu_hat,
                        const double deltaU);

// Antineutrino absorption on neutron (anul + p -> l+ + n)
MyOpacity nu_p_abs_stim(const double omega,
                        const double nb, const double temp,
                        const double lep_mass,
                        const double yp, const double yn,
                        const double mu_l, const double mu_hat,
                        const double deltaU);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_OPACITIES_OPACITIES_H_
