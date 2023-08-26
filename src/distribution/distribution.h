//
// Created by leonardo on 8/25/23.
//

#ifndef BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
#define BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_

#include "../bns_nurates.h"

// NuDistributionParams struct
// structure needed for storing the parameters of the 
// neutrino distribution function
struct NuDistributionParams {
  // optically thick
  double w_t;
  double temp_t;
  double eta_t;
  // optically thin
  double w_f;
  double temp_f;
  double c_f;
};
typedef struct NuDistributionParams NuDistributionParams;

/*===========================================================================*/

// opt_thick.c

// Neutrino distribution function in optically thick regime
double NuFThick(double omega, NuDistributionParams *distr_pars);

// Recover parameters of thick distribution function from M1 quantities
void CalculateThickParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars,
                                NuDistributionParams *out_distr_pars);

/*===========================================================================*/

// opt_thin.c

// Neutrino distribution function in optically thin regime
double NuFThin(double omega, NuDistributionParams *distr_pars);

// Recover parameters of thin distribution function from M1 quantities
void ThinFromM1(M1Quantities *M1_pars, NuDistributionParams *out_distr_pars);

/*===========================================================================*/

// total_distr_func.c

// Total neutrino distribution combining optically thick and thin regimes
double TotalNuF(double omega, NuDistributionParams *distr_pars);

// Recover parameters of thick and thin distribution function from M1 quantities
NuDistributionParams DistrParamsFromM1(M1Quantities *M1_pars, MyEOSParams *eos_pars);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_DISTRIBUTION_DISTRIBUTION_H_
