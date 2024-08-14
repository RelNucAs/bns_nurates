#ifndef BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
#define BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration.h
//  \brief header file for all integration routines

#include "bns_nurates.hpp"
#include "functions.hpp"

// routines for generating quadratures
void GaussLegendre(MyQuadrature* quad);
void GaussLegendreMultiD(MyQuadrature* quad);
void GaussLaguerre(MyQuadrature* quad);

// routines for saving and loading quadratures
void SaveQuadrature(char* filedir, MyQuadrature* quad);
void LoadQuadrature(char* filedir, MyQuadrature* quad);

/* Make a convolution between sampled function and weights
 *
 * Inputs:
 *    nx:        number of points in the quadrature scheme
 *    wtarray:  array of quadrature weights
 *    fnarray:  array of function values at corresponding quadrature positions
 */
KOKKOS_INLINE_FUNCTION
double DoIntegration(const int n, const double* wtarray, const double* fnarray)
{
    double integral = 0.;

    for (int i = 0; i < n; i++)
    {
        integral += wtarray[i] * fnarray[i];
    }

    return integral;
}

/* Integrate a function from 0 to inf using a Gauss-Legendre quadrature by
 * breaking the integral into two parts, from [0,t] and [t,inf]. The other
 * dimensions (if applicable) have finite limits.
 *
 * Note: Uses (4.4.2) of Numerical Recipes 3rd edition (Press et. al.) to
 * compute the integral from [t,inf] The quad struct should contain a
 * Gauss-Legendre quadrature in [0,1] only!
 *
 * Inputs:
 *    quad: A properly populated Gauss-Legendre quadrature struct from x1 = 0 to
 * x2 = 1. func: The function struct to be integrated t:    The value of points
 * at which to break the integral into two
 */
inline double GaussLegendreIntegrateZeroInf(MyQuadrature* quad, MyFunction* func,
                                     double t)
{

    double f1_x[quad->nx], f2_x[quad->nx], f_y[quad->ny], f_z[quad->nz];
    double w_y[quad->ny], w_z[quad->nz];
    // var3d var = var3d_default;
    double var[3];

    for (int k = 0; k < quad->nz; k++)
    {
        for (int j = 0; j < quad->ny; j++)
        {
            for (int i = 0; i < quad->nx; i++)
            {
                var[0] = t * quad->points[i];
                var[1] = quad->points[quad->nx + j];
                var[2] = quad->points[quad->nx + quad->ny + k];

                f1_x[i] = func->function(var, func->params);

                var[0] = t / quad->points[i];

                f2_x[i] = func->function(var, func->params) /
                          (quad->points[i] * quad->points[i]);
            }
            f_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_x) +
                          DoIntegration(quad->nx, quad->w, f2_x));
            w_y[j] = quad->w[quad->nx + j];
        }
        f_z[k] = DoIntegration(quad->ny, w_y, f_y);
        w_z[k] = quad->w[quad->nx + quad->ny + k];
    }

    return DoIntegration(quad->nz, f_z, w_z);
}

/* Compute 4 different integration results at once (emission & absorption for e
 * & points neutrinos) Use this function for integrating over kernel quantities
 *
 * Inputs:
 *      quad: a properly generated quadrature structure
 *      func: a special type of function structure which holds a function
 * (returns 4 values), eos parameters and all kernel parameters t:    the value
 * of points at which the integration is to be broken into two
 *
 * Outputs:
 *      the emissivities and absoptivities for e and points neutrinos (take care
 * of constant multiplications separately outside)
 */
inline MyKernelQuantity GaussLegendreIntegrateZeroInfSpecial(MyQuadrature* quad,
                                                      MyFunctionSpecial* func,
                                                      double t)
{

    double f1_em_e[quad->nx], f2_em_e[quad->nx];   // emissivity e neutrino
    double f1_abs_e[quad->nx], f2_abs_e[quad->nx]; // absoptivity e neutrino
    double f1_em_x[quad->nx], f2_em_x[quad->nx];   // emissivity mu/tau neutrino
    double f1_abs_x[quad->nx],
        f2_abs_x[quad->nx]; // absorptivity mu/tau neutrino

    double f_em_e_y[quad->ny];  // emissivity e neutrino
    double f_abs_e_y[quad->ny]; // absoptivity e neutrino
    double f_em_x_y[quad->ny];  // emissivity mu/tau neutrino
    double f_abs_x_y[quad->ny]; // absorptivity mu/tau neutrino

    double f_em_e_z[quad->ny];  // emissivity e neutrino
    double f_abs_e_z[quad->ny]; // absoptivity e neutrino
    double f_em_x_z[quad->ny];  // emissivity mu/tau neutrino
    double f_abs_x_z[quad->ny]; // absorptivity mu/tau neutrino

    double w_y[quad->ny], w_z[quad->nz];
    double var[3];

    for (int k = 0; k < quad->nz; k++)
    {
        for (int j = 0; j < quad->ny; j++)
        {
            for (int i = 0; i < quad->nx; i++)
            {
                var[0] = t * quad->points[i];
                var[1] = quad->points[quad->nx + j];
                var[2] = quad->points[quad->nx + quad->ny + k];

                MyKernelQuantity f1_vals =
                    func->function(var, func->eos_params, func->kernel_params);
                f1_em_e[i]  = f1_vals.em_e;
                f1_abs_e[i] = f1_vals.abs_e;
                f1_em_x[i]  = f1_vals.em_x;
                f1_abs_x[i] = f1_vals.abs_x;

                var[0] = t / quad->points[i];
                MyKernelQuantity f2_vals =
                    func->function(var, func->eos_params, func->kernel_params);
                f2_em_e[i] = f2_vals.em_e / (quad->points[i] * quad->points[i]);
                f2_abs_e[i] =
                    f2_vals.abs_e / (quad->points[i] * quad->points[i]);
                f2_em_x[i] = f2_vals.em_x / (quad->points[i] * quad->points[i]);
                f2_abs_x[i] =
                    f2_vals.abs_x / (quad->points[i] * quad->points[i]);
            }
            f_em_e_y[j]  = t * (DoIntegration(quad->nx, quad->w, f1_em_e) +
                               DoIntegration(quad->nx, quad->w, f2_em_e));
            f_abs_e_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_abs_e) +
                                DoIntegration(quad->nx, quad->w, f2_abs_e));
            f_em_x_y[j]  = t * (DoIntegration(quad->nx, quad->w, f1_em_x) +
                               DoIntegration(quad->nx, quad->w, f2_em_x));
            f_abs_x_y[j] = t * (DoIntegration(quad->nx, quad->w, f1_abs_x) +
                                DoIntegration(quad->nx, quad->w, f2_abs_x));
            w_y[j]       = quad->w[quad->nx + j];
        }
        f_em_e_z[k]  = DoIntegration(quad->ny, w_y, f_em_e_y);
        f_abs_e_z[k] = DoIntegration(quad->ny, w_y, f_abs_e_y);
        f_em_x_z[k]  = DoIntegration(quad->ny, w_y, f_em_x_y);
        f_abs_x_z[k] = DoIntegration(quad->ny, w_y, f_abs_x_y);
        w_z[k]       = quad->w[quad->nx + quad->ny + k];
    }

    MyKernelQuantity result = {.em_e  = DoIntegration(quad->nz, f_em_e_z, w_z),
                               .abs_e = DoIntegration(quad->nz, f_abs_e_z, w_z),
                               .em_x  = DoIntegration(quad->nz, f_em_x_z, w_z),
                               .abs_x =
                                   DoIntegration(quad->nz, f_abs_x_z, w_z)};

    return result;
}

/* Integrate a function from 0 to inf using a Gauss-Laguerre quadrature
 *
 * Inputs:
 *    quad: A properly populated Gauss-Laguerre quadrature struct
 *    func:    The function struct to be integrated
 */
inline double GaussLaguerreIntegrateZeroInf(MyQuadrature* quad, MyFunction* func)
{

    double f[quad->nx];

    if (quad->alpha == 0.)
    {
        for (int i = 0; i < quad->nx; i++)
        {
            f[i] = func->function(&quad->points[i],
                                  func->params); // * exp(quad->points[i]);
        }
    }
    else
    {
        for (int i = 0; i < quad->nx; i++)
            f[i] = func->function(
                &quad->points[i],
                func->params); // * exp(quad->points[i]) / pow(quad->points[i],
                               // quad->alpha);
    }

    return DoIntegration(quad->nx, quad->w, f);
}

/* Perform 2d integration of multiple functions using a Gauss-Legendre
 * quadrature
 *
 * quad:    must be a properly populated quadrature (upto 2d)
 * func:    the function(s) to be integrated
 * t:       the value at which to break the integral into two
 */
// @TODO: rewrite loops in row-wise order
inline
MyQuadratureIntegrand GaussLegendreIntegrateFixedSplit2D(MyQuadrature* quad,
                                                         MyFunctionMultiD* func,
                                                         double t)
{

    int num_integrands = func->my_quadrature_integrand.n;

    double f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    double f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    double w_y[quad->ny];
    double var[2];

    MyQuadratureIntegrand f1_vals, f2_vals;
    MyQuadratureIntegrand result;

    for (int j = 0; j < quad->ny; j++)
    {

        var[1] = t * quad->points[quad->nx + j];

        for (int i = 0; i < quad->nx; i++)
        {

            var[0]  = t * quad->points[i];
            f1_vals = func->function(var, func->params);

            var[0]  = t / quad->points[i];
            f2_vals = func->function(var, func->params);

            for (int k = 0; k < num_integrands; k++)
            {
                f1_x[k][i] = f1_vals.integrand[k];
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }
        }

        for (int k = 0; k < num_integrands; k++)
        {
            f1_y[k][j] = t * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                              DoIntegration(quad->nx, quad->w, f2_x[k]));
        }

        var[1] = t / quad->points[quad->nx + j];

        for (int i = 0; i < quad->nx; i++)
        {

            var[0]  = t * quad->points[i];
            f1_vals = func->function(var, func->params);

            var[0]  = t / quad->points[i];
            f2_vals = func->function(var, func->params);

            for (int k = 0; k < num_integrands; k++)
            {
                f1_x[k][i] = f1_vals.integrand[k];
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }
        }


        for (int k = 0; k < num_integrands; k++)
        {
            f2_y[k][j] = t * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                              DoIntegration(quad->nx, quad->w, f2_x[k]));
        }


        w_y[j] = quad->w[quad->nx + j];
    }

    for (int k = 0; k < num_integrands; k++)
    {
        result.integrand[k] = t * (DoIntegration(quad->ny, w_y, f1_y[k]) +
                                   DoIntegration(quad->ny, w_y, f2_y[k]));
    }

    return result;
}


/* Perform 1d integration of multiple functions using a Gauss-Legendre
 * quadrature
 *
 * quad:    must be a properly populated quadrature (upto 2d)
 * func:    the function(s) to be integrated
 * t:       the value at which to break the integral into two
 */
// @TODO: rewrite loops in row-wise order
inline
MyQuadratureIntegrand GaussLegendreIntegrateFixedSplit1D(MyQuadrature* quad,
                                                         MyFunctionMultiD* func,
                                                         double t)
{

    int num_integrands = func->my_quadrature_integrand.n;
    double f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    double var[2];

    MyQuadratureIntegrand f1_vals, f2_vals;

    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int i = 0; i < quad->nx; i++)
    {

        var[0]  = t * quad->points[i];
        f1_vals = func->function(var, func->params);

        var[0]  = t / quad->points[i];
        f2_vals = func->function(var, func->params);

        for (int k = 0; k < num_integrands; k++)
        {
            f1_x[k][i] = f1_vals.integrand[k];
            f2_x[k][i] =
                f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
        }
    }

    for (int k = 0; k < num_integrands; k++)
    {
        result.integrand[k] = t * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                                   DoIntegration(quad->nx, quad->w, f2_x[k]));
    }

    return result;
}

/* Perform 2d integration of multiple functions using a Gauss-Legendre
 * quadrature
 *
 * quad:    must be a properly populated quadrature (upto 2d)
 * func:    the function(s) to be integrated
 * t:       the value at which to break the integral into two
 */
inline
MyQuadratureIntegrand GaussLegendreIntegrate2D(MyQuadrature* quad,
                                               MyFunctionMultiD* func,
                                               double* tx, double* ty)
{

    int num_integrands = func->my_quadrature_integrand.n;

    double f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    double f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    double w_y[quad->ny];
    double var[2];

    MyQuadratureIntegrand result;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int j = 0; j < quad->ny; j++)
        {

            var[1] = ty[k] * quad->points[quad->nx + j];

            for (int i = 0; i < quad->nx; i++)
            {

                var[0] = tx[k] * quad->points[i];
                MyQuadratureIntegrand f1_vals =
                    func->function(var, func->params);
                f1_x[k][i] = f1_vals.integrand[k];

                var[0] = tx[k] / quad->points[i];
                MyQuadratureIntegrand f2_vals =
                    func->function(var, func->params);
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }

            f1_y[k][j] = tx[k] * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                                  DoIntegration(quad->nx, quad->w, f2_x[k]));

            var[1] = ty[k] / quad->points[quad->nx + j];

            for (int i = 0; i < quad->nx; i++)
            {

                var[0] = tx[k] * quad->points[i];
                MyQuadratureIntegrand f1_vals =
                    func->function(var, func->params);
                f1_x[k][i] = f1_vals.integrand[k];

                var[0] = tx[k] / quad->points[i];
                MyQuadratureIntegrand f2_vals =
                    func->function(var, func->params);
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }

            f2_y[k][j] = tx[k] * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                                  DoIntegration(quad->nx, quad->w, f2_x[k]));

            w_y[j] = quad->w[quad->nx + j];
        }

        result.integrand[k] = ty[k] * (DoIntegration(quad->ny, w_y, f1_y[k]) +
                                       DoIntegration(quad->ny, w_y, f2_y[k]));
    }

    return result;
}

/* Perform 1d integration of multiple functions using a Gauss-Legendre
 * quadrature
 *
 * quad:    must be a properly populated quadrature (upto 2d)
 * func:    the function(s) to be integrated
 * t:       the value at which to break the integral into two
 */
inline
MyQuadratureIntegrand
GaussLegendreIntegrate1D(MyQuadrature* quad, MyFunctionMultiD* func, double* t)
{

    int num_integrands = func->my_quadrature_integrand.n;
    double f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    double var[2];
    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int i = 0; i < quad->nx; i++)
        {

            var[0]                        = t[k] * quad->points[i];
            MyQuadratureIntegrand f1_vals = func->function(var, func->params);
            f1_x[k][i]                    = f1_vals.integrand[k];

            var[0]                        = t[k] / quad->points[i];
            MyQuadratureIntegrand f2_vals = func->function(var, func->params);
            f2_x[k][i] =
                f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
        }

        result.integrand[k] =
            t[k] * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                    DoIntegration(quad->nx, quad->w, f2_x[k]));
    }

    return result;
}

KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand GaussLegendreIntegrate2DMatrix(const MyQuadrature* quad,
                                                     const M1MatrixKokkos2D* mat, double t)
{
    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const double t_sqr = t * t;

    double w_i, w_j, w_ij;
    double x_i, x_j;
    double x2_i, x2_j, x2_ij;

    MyQuadratureIntegrand result = {0};

    for (int idx = 0; idx < total_num_species; idx++)
    {

        for (int i = 0; i < n; i++)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            for (int j = 0; j < n; j++)
            {

                x_j = quad->points[j];
                w_j = quad->w[j];

                w_ij  = w_i * w_j;
                x2_j  = x_j * x_j;
                x2_ij = x2_i * x2_j;

                result.integrand[0 + idx] +=
                    w_ij * (mat->m1_mat_em[idx][i][j] +
                            mat->m1_mat_em[idx][n + i][j] / x2_i +
                            mat->m1_mat_em[idx][i][n + j] / x2_j +
                            mat->m1_mat_em[idx][n + i][n + j] / x2_ij);

                result.integrand[total_num_species + idx] +=
                    w_ij * (mat->m1_mat_ab[idx][i][j] +
                            mat->m1_mat_ab[idx][n + i][j] / x2_i +
                            mat->m1_mat_ab[idx][i][n + j] / x2_j +
                            mat->m1_mat_ab[idx][n + i][n + j] / x2_ij);
            }
        }
    }

    for (int idx = 0; idx < num_integrands; idx++)
    {
        result.integrand[idx] = t_sqr * result.integrand[idx];
    }

    return result;
}


KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand GaussLegendreIntegrate2DMatrixFlatten(const MyQuadrature* quad,
                                                     const M1MatrixKokkos2DFlatten* mat, double t)
{
    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const double t_sqr = t * t;

    double w_i, w_j, w_ij;
    double x_i, x_j;
    double x2_i, x2_j, x2_ij;

    MyQuadratureIntegrand result = {0};

    for (int idx = 0; idx < total_num_species; idx++)
    {

        for (int i = 0; i < n; i++)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            for (int j = 0; j < n; j++)
            {

                x_j = quad->points[j];
                w_j = quad->w[j];

                w_ij  = w_i * w_j;
                x2_j  = x_j * x_j;
                x2_ij = x2_i * x2_j;

                result.integrand[0 + idx] +=
                    w_ij * (mat->m1_mat_em[idx][n_max * i + j] +
                            mat->m1_mat_em[idx][n_max * (n + i) + j] / x2_i +
                            mat->m1_mat_em[idx][n_max * i + n + j] / x2_j +
                            mat->m1_mat_em[idx][n_max * (n + i) + n + j] / x2_ij);

                result.integrand[total_num_species + idx] +=
                    w_ij * (mat->m1_mat_ab[idx][n_max * i + j] +
                            mat->m1_mat_ab[idx][n_max * (n + i) + j] / x2_i +
                            mat->m1_mat_ab[idx][n_max * i + n + j] / x2_j +
                            mat->m1_mat_ab[idx][n_max * (n + i) + n + j] / x2_ij);
            }
        }
    }

    for (int idx = 0; idx < num_integrands; idx++)
    {
        result.integrand[idx] = t_sqr * result.integrand[idx];
    }

    return result;
}


KOKKOS_INLINE_FUNCTION
void GaussLegendreIntegrate2DMatrixOnlyNumber(const MyQuadrature* quad,
                                                     const M1MatrixKokkos2D* mat, double t, MyQuadratureIntegrand* result_1, MyQuadratureIntegrand* result_2)
{
    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const double t_sqr = t * t;

    double w_i, w_j, w_ij;
    double x_i, x_j;
    double x2_i, x2_j, x2_ij;

    for (int idx = 0; idx < total_num_species; idx++)
    {

        for (int i = 0; i < n; i++)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            for (int j = 0; j < n; j++)
            {

                x_j = quad->points[j];
                w_j = quad->w[j];

                w_ij  = w_i * w_j;
                x2_j  = x_j * x_j;
                x2_ij = x2_i * x2_j;

                result_1->integrand[0 + idx] +=
                    w_ij * (mat->m1_mat_em[idx][i][j] +
                            mat->m1_mat_em[idx][n + i][j] / x2_i +
                            mat->m1_mat_em[idx][i][n + j] / x2_j +
                            mat->m1_mat_em[idx][n + i][n + j] / x2_ij);

                result_1->integrand[total_num_species + idx] +=
                    w_ij * (mat->m1_mat_ab[idx][i][j] +
                            mat->m1_mat_ab[idx][n + i][j] / x2_i +
                            mat->m1_mat_ab[idx][i][n + j] / x2_j +
                            mat->m1_mat_ab[idx][n + i][n + j] / x2_ij);
                
                result_2->integrand[0 + idx] +=
                    w_ij * (mat->m1_mat_em[idx][i][j] * t * x_i +
                            mat->m1_mat_em[idx][n + i][j] * (t / x_i) / x2_i +
                            mat->m1_mat_em[idx][i][n + j] * t * x_i / x2_j +
                            mat->m1_mat_em[idx][n + i][n + j] * (t / x_i) / x2_ij);

                result_2->integrand[total_num_species + idx] +=
                    w_ij * (mat->m1_mat_ab[idx][i][j] * t * x_i +
                            mat->m1_mat_ab[idx][n + i][j] * (t / x_i) / x2_i +
                            mat->m1_mat_ab[idx][i][n + j] * t * x_i / x2_j +
                            mat->m1_mat_ab[idx][n + i][n + j] * (t / x_i) / x2_ij);
            }
        }
    }

    for (int idx = 0; idx < num_integrands; idx++)
    {
        result_1->integrand[idx] = t_sqr * result_1->integrand[idx];
        result_2->integrand[idx] = t_sqr * result_2->integrand[idx];
    }

    return; // result;
}
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand GaussLegendreIntegrate1DMatrix(const MyQuadrature* quad,
                                                     const M1MatrixKokkos1D *mat, double t)
{
    const int n              = quad->nx;
    const int num_integrands = 12;

    double x_i, w_i, x2_i;
    double tmp_0, tmp_1;

    MyQuadratureIntegrand result = {0};

    for (int idx = 0; idx < num_integrands; idx++)
    {

        for (int i = 0; i < n; i++)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            //tmp_0 = mat[idx * n_max + i];
            //tmp_1 = mat[idx * n_max + n + i];

            result.integrand[idx] +=
                quad->w[i] * (mat->m1_mat[idx][i] + mat->m1_mat[idx][n + i] / x2_i);
        }
    }

    for (int idx = 0; idx < num_integrands; idx++)
    {
        result.integrand[idx] = t * result.integrand[idx];
    }

    return result;
}


/* Perform 1d integration (finite interval) of multiple functions using a
 * Gauss-Legendre quadrature
 *
 * quad:    must be a properly populated 1d quadrature generated from between
 * interval of integration func:    the function(s) to be integrated
 */
inline
MyQuadratureIntegrand
GaussLegendreIntegrate1DFiniteInterval(MyQuadrature* quad,
                                       MyFunctionMultiD* func, double* t)
{
    (void)t;

    int num_integrands = func->my_quadrature_integrand.n;
    double f_x[num_integrands][quad->nx];
    double var[2];
    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int k = 0; k < num_integrands; ++k)
    {
        for (int i = 0; i < quad->nx; i++)
        {

            var[0]                        = quad->points[i];
            MyQuadratureIntegrand f1_vals = func->function(var, func->params);
            f_x[k][i]                     = f1_vals.integrand[k];
        }
        result.integrand[k] = DoIntegration(quad->nx, quad->w, f_x[k]);
    }

    return result;
}

/* Perform 2d integration (finite interval) of multiple functions using a
 * Gauss-Legendre quadrature
 *
 * The first variable is integrated within a finite interval, the second
 * variable is integrated from 0 to inf
 *
 * quad:    must be a properly populated 2d quadrature generated from between
 * interval of integration func:    the function(s) to be integrated
 */
inline
MyQuadratureIntegrand GaussLegendreIntegrate2DFiniteInterval(
    MyQuadrature* quad, MyFunctionMultiD* func, double* tx, double* ty)
{
    (void)tx;

    int num_integrands = func->my_quadrature_integrand.n;

    double f1_x[num_integrands][quad->nx];
    double f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    double w_y[quad->ny];
    double var[2];

    MyQuadratureIntegrand result;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int j = 0; j < quad->ny; j++)
        {

            var[1] = ty[k] * quad->points[quad->nx + j];
            for (int i = 0; i < quad->nx; i++)
            {
                var[0] = quad->points[i];
                MyQuadratureIntegrand f1_vals =
                    func->function(var, func->params);
                f1_x[k][i] = f1_vals.integrand[k];
            }

            f1_y[k][j] = DoIntegration(quad->nx, quad->w, f1_x[k]);

            var[1] = ty[k] / quad->points[quad->nx + j];
            for (int i = 0; i < quad->nx; i++)
            {
                var[0] = quad->points[i];
                MyQuadratureIntegrand f1_vals =
                    func->function(var, func->params);
                f1_x[k][i] = f1_vals.integrand[k];
            }

            f2_y[k][j] = DoIntegration(quad->nx, quad->w, f1_x[k]);

            w_y[j] = quad->w[quad->nx + j];
        }

        result.integrand[k] = ty[k] * (DoIntegration(quad->ny, w_y, f1_y[k]) +
                                       DoIntegration(quad->ny, w_y, f2_y[k]));
    }

    return result;
}

#endif // BNS_NURATES_SRC_INTEGRATION_INTEGRATION_H_
