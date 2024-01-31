//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  nu_abs_em_brem.c
//  \brief computes opacities for the bremsstrahlung reaction

#include <math.h>
#include "bns_nurates.h"
#include "constants.h"
#include "kernels.h"
#include "functions.h"
#include "integration.h"

MyKernelQuantity BremEmissivityAbsorptivityIntegrandFermi(double *omega_prime, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {
  (void)omega_prime;
  (void)my_eos_params;
  (void)my_kernel_params;

  MyKernelQuantity result = {.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
  return result;
}

MyKernelQuantity BremOpacitiesFermi(MyQuadrature *quad, MyEOSParams *my_eos_params, MyKernelParams *my_kernel_params) {

  double t = 1.5 * my_eos_params->temp;
  MyFunctionSpecial func;
  func.function = &BremEmissivityAbsorptivityIntegrandFermi;
  func.dim = 1;
  func.eos_params = my_eos_params;
  func.kernel_params = my_kernel_params;

  MyKernelQuantity opacities = GaussLegendreIntegrateZeroInfSpecial(quad, &func, t);

  return opacities;
}

void BremOpacitiesTable(MyQuadrature *quad, GreyOpacityParams *grey_pars, double t, M1Matrix *out) {
  double nu, nu_bar;

  const int n = quad->nx;

  MyKernelOutput brem_ker;

  grey_pars->kernel_pars.brem_kernel_params.l = 0;

  for (int i = 0; i < quad->nx; i++) {

    for (int j = i; j < quad->nx; j++) {
      
      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t * quad->points[j];
     
      // compute the brem kernels
      grey_pars->kernel_pars.brem_kernel_params.omega = nu;
      grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_bar;

      brem_ker = BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
      
      out->m1_mat_em[0][i][j] = brem_ker.em[0];
      out->m1_mat_em[0][j][i] = brem_ker.em[0];

      out->m1_mat_ab[0][i][j] = brem_ker.abs[0];
      out->m1_mat_ab[0][j][i] = brem_ker.abs[0];


      // energies and parameters
      nu = t / quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the brem kernels
      grey_pars->kernel_pars.brem_kernel_params.omega = nu;
      grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_bar;
      
      brem_ker = BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
      
      out->m1_mat_em[0][n+i][n+j] = brem_ker.em[0];
      out->m1_mat_em[0][n+j][n+i] = brem_ker.em[0];

      out->m1_mat_ab[0][n+i][n+j] = brem_ker.abs[0];
      out->m1_mat_ab[0][n+j][n+i] = brem_ker.abs[0];
      
    }

    for (int j = 0; j < quad->nx; j++) {

      // energies and parameters
      nu = t * quad->points[i];
      nu_bar = t / quad->points[j];

      // compute the brem kernels
      grey_pars->kernel_pars.brem_kernel_params.omega = nu;
      grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_bar;
      
      brem_ker = BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
      
      out->m1_mat_em[0][i][n+j] = brem_ker.em[0];
      out->m1_mat_em[0][n+j][i] = brem_ker.em[0];

      out->m1_mat_ab[0][i][n+j] = brem_ker.abs[0];
      out->m1_mat_ab[0][n+j][i] = brem_ker.abs[0];

    }
  }

  return;
}

