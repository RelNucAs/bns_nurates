// \file nu_gray_total.c
// \brief Computation of total grey coefficients

#include "opacities.h"
#include "../bns_nurates.h"

// Total emission coefficients
SourceCoeffs EmissionCoeffs(MyEOSParams *eos_pars) {
  SourceCoeffs em_beta = NuclAbsEmissionCoeffs(eos_pars);
  
  // @TODO: substitute with the proper function
  SourceCoeffs em_pair = {.R_nue = 0., .R_anue = 0., .R_nux = 0.,
                          .Q_nue = 0., .Q_anue = 0., .Q_nux = 0.};

  // @TODO: substitute with the proper function
  SourceCoeffs em_brem = {.R_nue = 0., .R_anue = 0., .R_nux = 0.,
                          .Q_nue = 0., .Q_anue = 0., .Q_nux = 0.};

  SourceCoeffs em_tot = {.R_nue  = em_beta.R_nue  + em_pair.R_nue  + em_brem.R_nue,
                         .R_anue = em_beta.R_anue + em_pair.R_anue + em_brem.R_anue,
                         //.R_num  = em_beta.R_num  + em_pair.R_num  + em_brem.R_num,
                         //.R_anum = em_beta.R_anum + em_pair.R_anum + em_brem.R_anum,
                         .R_nux  = em_beta.R_nux  + em_pair.R_nux  + em_brem.R_nux,
                         .Q_nue  = em_beta.Q_nue  + em_pair.Q_nue  + em_brem.Q_nue,
                         .Q_anue = em_beta.Q_anue + em_pair.Q_anue + em_brem.Q_anue,
                         //.Q_num  = em_beta.Q_num  + em_pair.Q_num  + em_brem.Q_num,
                         //.Q_anum = em_beta.Q_anum + em_pair.Q_anum + em_brem.Q_anum,
                         .Q_nux  = em_beta.Q_nux  + em_pair.Q_nux  + em_brem.Q_nux};

  return em_tot;  
}


// Total opacity coefficients
SourceCoeffs OpacityCoeffs(GreyOpacityParams *grey_pars) {
  SourceCoeffs ab_beta = NuclAbsOpacityCoeffs(grey_pars);
  
  // @TODO: substitute with the proper function
  SourceCoeffs ab_pair = {.R_nue = 0., .R_anue = 0., .R_nux = 0.,
                          .Q_nue = 0., .Q_anue = 0., .Q_nux = 0.};

  // @TODO: substitute with the proper function
  SourceCoeffs ab_brem = {.R_nue = 0., .R_anue = 0., .R_nux = 0.,
                          .Q_nue = 0., .Q_anue = 0., .Q_nux = 0.};
             
  SourceCoeffs ab_tot = {.R_nue  = ab_beta.R_nue  + ab_pair.R_nue  + ab_brem.R_nue,
                         .R_anue = ab_beta.R_anue + ab_pair.R_anue + ab_brem.R_anue,
                         //.R_num  = ab_beta.R_num  + ab_pair.R_num  + ab_brem.R_num,
                         //.R_anum = ab_beta.R_anum + ab_pair.R_anum + ab_brem.R_anum,
                         .R_nux  = ab_beta.R_nux  + ab_pair.R_nux  + ab_brem.R_nux,
                         .Q_nue  = ab_beta.Q_nue  + ab_pair.Q_nue  + ab_brem.Q_nue,
                         .Q_anue = ab_beta.Q_anue + ab_pair.Q_anue + ab_brem.Q_anue,
                         //.Q_num  = ab_beta.Q_num  + ab_pair.Q_num  + ab_brem.Q_num,
                         //.Q_anum = ab_beta.Q_anum + ab_pair.Q_anum + ab_brem.Q_anum,
                         .Q_nux  = ab_beta.Q_nux  + ab_pair.Q_nux  + ab_brem.Q_nux};

  return ab_tot;  
}

// Total scattering coefficients
SourceCoeffs ScatteringCoeffs(GreyOpacityParams *grey_pars) {
  return IsoScattCoeffs(grey_pars);  
}