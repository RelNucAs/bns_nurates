// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_integrate_kernels.c
//  \brief special integration routines for integrating over opacity kernels
//         multiple integrals are performed inside a single routine

#include <math.h>
#include "integration.h"

/* Compute 4 different integration results at once (emission & absorption for e & x neutrinos)
 * Use this function for integrating over kernel quantities
 *
 * Inputs:
 *      quad: a properly generated quadrature structure
 *      func: a special type of function structure which holds a function (returns 4 values), eos parameters and all kernel parameters
 *      t:    the value of x at which the integration is to be broken into two
 *
 * Outputs:
 *      the emissivities and absoptivities for e and x neutrinos (take care of constant multiplications separately outside)
 */
MyOpacityQuantity GaussLegendreIntegrateZeroInfSpecial(MyQuadrature *quad, MyFunctionSpecial *func, double t) {

  double f1_em_e[quad->n], f2_em_e[quad->n];
  double f1_abs_e[quad->n], f2_abs_e[quad->n];
  double f1_em_x[quad->n], f2_em_x[quad->n];
  double f1_abs_x[quad->n], f2_abs_x[quad->n];

  for (int i = 0; i < quad->n; i++) {
    func->kernel_params->brem_kernel_params.omega_prime = quad->x[i];
    func->kernel_params->pair_kernel_params.omega_prime = quad->x[i];

    MyOpacityQuantity f1_vals = func->function(t * quad->x[i], func->eos_params, func->kernel_params);
    MyOpacityQuantity f2_vals = func->function(t / quad->x[i], func->eos_params, func->kernel_params);

    f1_em_e[i] = f1_vals.em_e;
    f1_abs_e[i] = f1_vals.abs_e;
    f1_em_x[i] = f1_vals.em_x;
    f1_abs_x[i] = f1_vals.abs_x;

    f2_em_e[i] = f2_vals.em_e / (quad->x[i] * quad->x[i]);
    f2_abs_e[i] = f2_vals.abs_e / (quad->x[i] * quad->x[i]);
    f2_em_x[i] = f2_vals.em_x / (quad->x[i] * quad->x[i]);
    f2_abs_x[i] = f2_vals.abs_x / (quad->x[i] * quad->x[i]);
  }

  double result_em_e = t * (DoIntegration(quad->n, quad->w, f1_em_e) + DoIntegration(quad->n, quad->w, f2_em_e));
  double result_abs_e = t * (DoIntegration(quad->n, quad->w, f1_abs_e) + DoIntegration(quad->n, quad->w, f2_abs_e));
  double result_em_x = t * (DoIntegration(quad->n, quad->w, f1_em_x) + DoIntegration(quad->n, quad->w, f2_em_x));
  double result_abs_x = t * (DoIntegration(quad->n, quad->w, f1_abs_x) + DoIntegration(quad->n, quad->w, f2_abs_x));

  MyOpacityQuantity result = {.em_e = result_em_e, .abs_e = result_abs_e, .em_x = result_em_x, .abs_x = result_abs_x};

  return result;

}