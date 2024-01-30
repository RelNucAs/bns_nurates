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


void PairOpacitiesTable(MyQuadrature *quad, MyEOSParams *eos_pars, MyKernelParams *kernel_pars, double t, M1Matrix *out) {
  double nu, nu_bar;

  const int n = quad->nx;

  MyKernelQuantity pair_1, pair_2;

  for (int idx = 0; idx < total_num_species; idx++) {
    out->m1_mat_ab[idx] = (double **) malloc(sizeof(double *) * 2 * quad->nx);
    out->m1_mat_em[idx] = (double **) malloc(sizeof(double *) * 2 * quad->nx);

    for (int i = 0; i < 2 * quad->nx; i++) {
      out->m1_mat_ab[idx][i] = (double *) malloc(sizeof(double) * 2 * quad->nx);
      out->m1_mat_em[idx][i] = (double *) malloc(sizeof(double) * 2 * quad->nx);
    }
  }

  for (int i = 0; i < quad->nx; i++) {

    for (int j = i; j < quad->nx; j++) {
      
      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t * quad->points[j];
     
      // compute the pair kernels
      kernel_pars->pair_kernel_params.omega = nu;
      kernel_pars->pair_kernel_params.omega_prime = nu_bar;
      kernel_pars->pair_kernel_params.cos_theta = 1.;
      kernel_pars->pair_kernel_params.filter = 0.;
      kernel_pars->pair_kernel_params.lmax = 0;
      kernel_pars->pair_kernel_params.mu = 1.;
      kernel_pars->pair_kernel_params.mu_prime = 1.;

      PairKernelsM1Test(eos_pars, &kernel_pars->pair_kernel_params, &pair_1, &pair_2);
      
      out->m1_mat_em[id_nue][i][j] = pair_1.em_e;
      out->m1_mat_em[id_nue][j][i] = pair_2.em_e;

      out->m1_mat_em[id_anue][i][j] = pair_1.em_e;
      out->m1_mat_em[id_anue][j][i] = pair_2.em_e;

      out->m1_mat_em[id_nux][i][j] = pair_1.em_x;
      out->m1_mat_em[id_nux][j][i] = pair_2.em_x;

      out->m1_mat_em[id_anux][i][j] = pair_1.em_x;
      out->m1_mat_em[id_anux][j][i] = pair_2.em_x;    


      out->m1_mat_ab[id_nue][i][j] = pair_1.abs_e;
      out->m1_mat_ab[id_nue][j][i] = pair_2.abs_e;

      out->m1_mat_ab[id_anue][i][j] = pair_1.abs_e;
      out->m1_mat_ab[id_anue][j][i] = pair_2.abs_e;

      out->m1_mat_ab[id_nux][i][j] = pair_1.abs_x;
      out->m1_mat_ab[id_nux][j][i] = pair_2.abs_x;

      out->m1_mat_ab[id_anux][i][j] = pair_1.abs_x;
      out->m1_mat_ab[id_anux][j][i] = pair_2.abs_x;
   

      // energies and parameters
      nu = t / quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      kernel_pars->pair_kernel_params.omega = nu;
      kernel_pars->pair_kernel_params.omega_prime = nu_bar;
      
      PairKernelsM1Test(eos_pars, &kernel_pars->pair_kernel_params, &pair_1, &pair_2);
      
      out->m1_mat_em[id_nue][n+i][n+j] = pair_1.em_e;
      out->m1_mat_em[id_nue][n+j][n+i] = pair_2.em_e;
      
      out->m1_mat_em[id_anue][n+i][n+j] = pair_1.em_e;
      out->m1_mat_em[id_anue][n+j][n+i] = pair_2.em_e;

      out->m1_mat_em[id_nux][n+i][n+j] = pair_1.em_x;
      out->m1_mat_em[id_nux][n+j][n+i] = pair_2.em_x;

      out->m1_mat_em[id_anux][n+i][n+j] = pair_1.em_x;
      out->m1_mat_em[id_anux][n+j][n+i] = pair_2.em_x;    

      out->m1_mat_ab[id_nue][n+i][n+j] = pair_1.abs_e;
      out->m1_mat_ab[id_nue][n+j][n+i] = pair_2.abs_e;

      out->m1_mat_ab[id_anue][n+i][n+j] = pair_1.abs_e;
      out->m1_mat_ab[id_anue][n+j][n+i] = pair_2.abs_e;

      out->m1_mat_ab[id_nux][n+i][n+j] = pair_1.abs_x;
      out->m1_mat_ab[id_nux][n+j][n+i] = pair_2.abs_x;

      out->m1_mat_ab[id_anux][n+i][n+j] = pair_1.abs_x;
      out->m1_mat_ab[id_anux][n+j][n+i] = pair_2.abs_x;
    }

    for (int j = 0; j < quad->nx; j++) {

      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the pair kernels
      kernel_pars->pair_kernel_params.omega = nu;
      kernel_pars->pair_kernel_params.omega_prime = nu_bar;
      
      PairKernelsM1Test(eos_pars, &kernel_pars->pair_kernel_params, &pair_1, &pair_2);
      
      out->m1_mat_em[id_nue][i][n+j] = pair_1.em_e;
      out->m1_mat_em[id_nue][n+j][i] = pair_2.em_e;
      
      out->m1_mat_em[id_anue][i][n+j] = pair_1.em_e;
      out->m1_mat_em[id_anue][n+j][i] = pair_2.em_e;

      out->m1_mat_em[id_nux][i][n+j] = pair_1.em_x;
      out->m1_mat_em[id_nux][n+j][i] = pair_2.em_x;

      out->m1_mat_em[id_anux][i][n+j] = pair_1.em_x;
      out->m1_mat_em[id_anux][n+j][i] = pair_2.em_x;    

      out->m1_mat_ab[id_nue][i][n+j] = pair_1.abs_e;
      out->m1_mat_ab[id_nue][n+j][i] = pair_2.abs_e;

      out->m1_mat_ab[id_anue][i][n+j] = pair_1.abs_e;
      out->m1_mat_ab[id_anue][n+j][i] = pair_2.abs_e;

      out->m1_mat_ab[id_nux][i][n+j] = pair_1.abs_x;
      out->m1_mat_ab[id_nux][n+j][i] = pair_2.abs_x;

      out->m1_mat_ab[id_anux][i][n+j] = pair_1.abs_x;
      out->m1_mat_ab[id_anux][n+j][i] = pair_2.abs_x;
    }
  }

  return;
}
