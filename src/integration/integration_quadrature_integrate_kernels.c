// ================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_integrate_kernels.c
//  \brief special integration routines for integrating over opacity kernels
//         multiple integrals are performed inside a single routine

#include <math.h>
#include "integration.h"

/* Compute 4 different integration results at once (emission & absorption for e & points neutrinos)
 * Use this function for integrating over kernel quantities
 *
 * Inputs:
 *      quad: a properly generated quadrature structure
 *      func: a special type of function structure which holds a function (returns 4 values), eos parameters and all kernel parameters
 *      t:    the value of points at which the integration is to be broken into two
 *
 * Outputs:
 *      the emissivities and absoptivities for e and points neutrinos (take care of constant multiplications separately outside)
 */
MyOpacityQuantity GaussLegendreIntegrateZeroInfSpecial(MyQuadrature *quad, MyFunctionSpecial *func, double t) {

  double f1_em_e[quad->nx], f2_em_e[quad->nx];    // emissivity e neutrino
  double f1_abs_e[quad->nx], f2_abs_e[quad->nx];  // absoptivity e neutrino
  double f1_em_x[quad->nx], f2_em_x[quad->nx];    // emissivity mu/tau neutrino
  double f1_abs_x[quad->nx], f2_abs_x[quad->nx];  // absorptivity mu/tau neutrino

  double f_em_e_y[quad->ny];   // emissivity e neutrino
  double f_abs_e_y[quad->ny];  // absoptivity e neutrino
  double f_em_x_y[quad->ny];    // emissivity mu/tau neutrino
  double f_abs_x_y[quad->ny];  // absorptivity mu/tau neutrino

  double f_em_e_z[quad->ny];    // emissivity e neutrino
  double f_abs_e_z[quad->ny];  // absoptivity e neutrino
  double f_em_x_z[quad->ny];    // emissivity mu/tau neutrino
  double f_abs_x_z[quad->ny];  // absorptivity mu/tau neutrino

  double w_y[quad->ny], w_z[quad->nz];
  double var[3];

  for (int k = 0; k < quad->nz; k++) {
    for (int j = 0; j < quad->ny; j++) {
      for (int i = 0; i < quad->nx; i++) {
        var[0] = t * quad->points[i];
        var[1] = quad->points[quad->nx + j];
        var[2] = quad->points[quad->nx + quad->ny + k];

        MyOpacityQuantity f1_vals = func->function(var, func->eos_params, func->kernel_params);
        f1_em_e[i] = f1_vals.em_e;
        f1_abs_e[i] = f1_vals.abs_e;
        f1_em_x[i] = f1_vals.em_x;
        f1_abs_x[i] = f1_vals.abs_x;

        var[0] = t / quad->points[i];
        MyOpacityQuantity f2_vals = func->function(var, func->eos_params, func->kernel_params);
        f2_em_e[i] = f2_vals.em_e / (quad->points[i] * quad->points[i]);
        f2_abs_e[i] = f2_vals.abs_e / (quad->points[i] * quad->points[i]);
        f2_em_x[i] = f2_vals.em_x / (quad->points[i] * quad->points[i]);
        f2_abs_x[i] = f2_vals.abs_x / (quad->points[i] * quad->points[i]);
      }
      f_em_e_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_em_e) + DoIntegration(quad->nx, quad->w, f2_em_e));
      f_abs_e_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_abs_e) + DoIntegration(quad->nx, quad->w, f2_abs_e));
      f_em_x_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_em_x) + DoIntegration(quad->nx, quad->w, f2_em_x));
      f_abs_x_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_abs_x) + DoIntegration(quad->nx, quad->w, f2_abs_x));
      w_y[j] = quad->w[quad->nx + j];
    }
    f_em_e_z[k] = DoIntegration(quad->ny, w_y, f_em_e_y);
    f_abs_e_z[k] = DoIntegration(quad->ny, w_y, f_abs_e_y);
    f_em_x_z[k] = DoIntegration(quad->ny, w_y, f_em_x_y);
    f_abs_x_z[k] = DoIntegration(quad->ny, w_y, f_abs_x_y);
    w_z[k] = quad->w[quad->nx + quad->ny + k];
  }

  MyOpacityQuantity result = {.em_e = DoIntegration(quad->nz, f_em_e_z, w_z), .abs_e = DoIntegration(quad->nz, f_abs_e_z, w_z),
      .em_x = DoIntegration(quad->nz, f_em_x_z, w_z), .abs_x = DoIntegration(quad->nz, f_abs_x_z, w_z)};

  return result;

}


