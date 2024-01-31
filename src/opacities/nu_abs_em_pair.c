//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file nu_abs_em_pair.c
//  \brief calculate emissivity and inverse mean free path for the pair process

#include <math.h>
#include "bns_nurates.h"
#include "constants.h"
#include "kernels.h"
#include "functions.h"
#include "integration.h"

// Calculate the 4 integrands needed for computing pair opacities (assuming no dependence on mu' and phi')
// var[3]: [omega' mu' phi'] can in generally have 3 variables but the kernel we are considering has already
// been integrated over mu' and phi' and therefore only depends on the anti-neutrino energy (nu') and EOS/Kernel
// parameters.
MyKernelQuantity PairEmissivityAbsorptivityIntegrandFermi(double *var, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  my_kernel_params->pair_kernel_params.omega_prime = var[0];
  MyKernelQuantity pair_kernel = PairKernelsPhiMuIntegrated(my_eos_params, &my_kernel_params->pair_kernel_params);

  double fermi_e = FermiDistr(var[0], my_eos_params->temp, my_eos_params->mu_e / my_eos_params->temp);
  double fermi_x = FermiDistr(var[0], my_eos_params->temp, 0.);

  MyKernelQuantity result;

  result.em_e = var[0] * var[0] * (1. - fermi_e) * pair_kernel.em_e;
  result.em_x = var[0] * var[0] * (1. - fermi_x) * pair_kernel.em_x;

  result.abs_e = var[0] * var[0] * fermi_e * pair_kernel.abs_e;
  result.abs_x = var[0] * var[0] * fermi_x * pair_kernel.abs_x;

  return result;

}

MyKernelQuantity PairOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  double t = 1.5 * my_eos_params->temp;
  MyFunctionSpecial func;
  func.function = &PairEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyKernelQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}


void PairOpacitiesTable(MyQuadrature *quad, GreyOpacityParams *grey_pars, double t, M1Matrix *out) {
  double nu, nu_bar;

  const int n = quad->nx;

  MyKernelOutput pair_1, pair_2;

  for (int i = 0; i < quad->nx; i++) {

    for (int j = i; j < quad->nx; j++) {
      
      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t * quad->points[j];
     
      // compute the pair kernels
      grey_pars->kernel_pars.pair_kernel_params.omega = nu;
      grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_bar;
      grey_pars->kernel_pars.pair_kernel_params.cos_theta = 1.;
      grey_pars->kernel_pars.pair_kernel_params.filter = 0.;
      grey_pars->kernel_pars.pair_kernel_params.lmax = 0;
      grey_pars->kernel_pars.pair_kernel_params.mu = 1.;
      grey_pars->kernel_pars.pair_kernel_params.mu_prime = 1.;

      PairKernelsM1Test(&grey_pars->eos_pars, &grey_pars->kernel_pars.pair_kernel_params, &pair_1, &pair_2);
      
      for (int idx = 0; idx < total_num_species; idx ++) {
        out->m1_mat_em[idx][i][j] = pair_1.em[idx];
        out->m1_mat_em[idx][j][i] = pair_2.em[idx];

        out->m1_mat_ab[idx][i][j] = pair_1.abs[idx];
        out->m1_mat_ab[idx][j][i] = pair_2.abs[idx];
      }
        
 

      // energies and parameters
      nu = t / quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      grey_pars->kernel_pars.pair_kernel_params.omega = nu;
      grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_bar;
      
      PairKernelsM1Test(&grey_pars->eos_pars, &grey_pars->kernel_pars.pair_kernel_params, &pair_1, &pair_2);
      
      for (int idx = 0; idx < total_num_species; idx ++) {
        out->m1_mat_em[idx][n+i][n+j] = pair_1.em[idx];
        out->m1_mat_em[idx][n+j][n+i] = pair_2.em[idx];

        out->m1_mat_ab[idx][n+i][n+j] = pair_1.abs[idx];
        out->m1_mat_ab[idx][n+j][n+i] = pair_2.abs[idx];
      }
    }

    for (int j = 0; j < quad->nx; j++) {

      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      grey_pars->kernel_pars.pair_kernel_params.omega = nu;
      grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_bar;
      
      PairKernelsM1Test(&grey_pars->eos_pars, &grey_pars->kernel_pars.pair_kernel_params, &pair_1, &pair_2);
      
      for (int idx = 0; idx < total_num_species; idx ++) {
        out->m1_mat_em[idx][i][n+j] = pair_1.em[idx];
        out->m1_mat_em[idx][n+j][i] = pair_2.em[idx];

        out->m1_mat_ab[idx][i][n+j] = pair_1.abs[idx];
        out->m1_mat_ab[idx][n+j][i] = pair_2.abs[idx];
      }
    }
  }

  return;
}
