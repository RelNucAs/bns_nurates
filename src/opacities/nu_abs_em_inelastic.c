//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file nu_abs_em_inelastic_scatt.c
//  \brief calculate emissivity and inverse mean free path for the inelastic scattering on electrons/positrons

#include <math.h>
#include "bns_nurates.h"
#include "constants.h"
#include "kernels.h"
#include "functions.h"
#include "integration.h"

void InelasticOpacitiesTable(MyQuadrature *quad, GreyOpacityParams *grey_pars, double t, M1Matrix *out) {
  double nu, nu_bar;

  const int n = quad->nx;

  MyKernelOutput inel_1, inel_2;

  for (int i = 0; i < n; i++) {

    for (int j = i; j < n; j++) {
      
      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t * quad->points[j];
     
      // compute the pair kernels
      grey_pars->kernel_pars.inelastic_kernel_params.omega = nu;
      grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

      CrossedInelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars, &inel_1, &inel_2);

      for (int idx = 0; idx < total_num_species; idx++) {
        out->m1_mat_em[idx][i][j] = inel_1.em[idx];
        out->m1_mat_em[idx][j][i] = inel_2.em[idx];
      
        out->m1_mat_ab[idx][i][j] = inel_1.abs[idx];
        out->m1_mat_ab[idx][j][i] = inel_2.abs[idx];
      }

      // energies and parameters
      nu = t / quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      grey_pars->kernel_pars.inelastic_kernel_params.omega = nu;
      grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
      
      CrossedInelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars, &inel_1, &inel_2);

      for (int idx = 0; idx < total_num_species; idx++) {
        out->m1_mat_em[idx][n+i][n+j] = inel_1.em[idx];
        out->m1_mat_em[idx][n+j][n+i] = inel_2.em[idx];
      
        out->m1_mat_ab[idx][n+i][n+j] = inel_1.abs[idx];
        out->m1_mat_ab[idx][n+j][n+i] = inel_2.abs[idx];
      }   
    
    }

    for (int j = 0; j < n; j++) {

      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      grey_pars->kernel_pars.inelastic_kernel_params.omega = nu;
      grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
      
      CrossedInelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars, &inel_1, &inel_2);

      for (int idx = 0; idx < total_num_species; idx++) {
        out->m1_mat_em[idx][i][n+j] = inel_1.em[idx];
        out->m1_mat_em[idx][n+j][i] = inel_2.em[idx];
      
        out->m1_mat_ab[idx][i][n+j] = inel_1.abs[idx];
        out->m1_mat_ab[idx][n+j][i] = inel_2.abs[idx];
      }  
    }
  }

  return;
}
