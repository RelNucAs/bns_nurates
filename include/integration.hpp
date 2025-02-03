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

#define num_max_integrands 10

/* Generate Gauss-Legendre quadratures in [x1,x2].
 *
 * Note: This routine can only generate quadrature in 1d.
 *
 * For this routine to generate data successfully, quad struct
 * must have dim, type, nx, x1, x2 metadata populated.
 *
 * Inputs:
 *      quad: A MyQuadrature structure to hold quadratures.
 *            This already contains metadata for the quadrature, the routine
 *            only populates the quadrature points and weights
 */
inline void GaussLegendre(MyQuadrature* quad)
{
    constexpr BS_REAL half = 0.5;

    const double kEps = 1.0e-10; // 1.0e-14;

    BS_ASSERT(quad->dim == 1);
    BS_ASSERT(quad->type == kGauleg);

    double z1, z, xm, xl, pp, p3, p2, p1;

    int n = quad->nx;
    int m = (n + 1) / 2;
    xm    = half * (quad->x2 + quad->x1);
    xl    = half * (quad->x2 - quad->x1);

    for (int i = 0; i < m; ++i)
    {
        z = cos(kBS_Pi * (i + 0.75) / (n + 0.5));
        do
        {
            p1 = 1.0;
            p2 = 0.0;

            for (int j = 0; j < n; ++j)
            {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }

            pp = n * (z * p1 - p2) / (POW2(z) - 1.0);
            z1 = z;
            z  = z1 - p1 / pp;

        } while (fabs(z - z1) > kEps);

        quad->points[i]         = xm - xl * z;
        quad->points[n - 1 - i] = xm + xl * z;
        quad->w[i]              = 2.0 * xl / ((1.0 - POW2(z)) * POW2(pp));
        quad->w[n - 1 - i]      = quad->w[i];

        BS_ASSERT(isfinite(quad->points[i]), "quad->points[%d]=%e", i,
                  quad->points[i]);
        BS_ASSERT(quad->points[i] > -1., "quad->points[%d]=%e", i,
                  quad->points[i]);
        BS_ASSERT(quad->points[i] < +1., "quad->points[%d]=%e", i,
                  quad->points[i]);

        BS_ASSERT(isfinite(quad->points[n - 1 - i]), "quad->points[%d]=%e",
                  n - 1 - i, quad->points[n - 1 - i]);
        BS_ASSERT(quad->points[n - 1 - i] > -1., "quad->points[%d]=%e",
                  n - 1 - i, quad->points[n - 1 - i]);
        BS_ASSERT(quad->points[n - 1 - i] < +1., "quad->points[%d]=%e",
                  n - 1 - i, quad->points[n - 1 - i]);


        BS_ASSERT(isfinite(quad->w[i]) && quad->w[i] > 0. && quad->w[i] < 1.);
        BS_ASSERT(isfinite(quad->w[n - 1 - i]) && quad->w[n - 1 - i] > 0. &&
                  quad->w[n - 1 - i] < 1.);
    }
}

/* Make a convolution between sampled function and weights
 *
 * Inputs:
 *    nx:        number of points in the quadrature scheme
 *    wtarray:  array of quadrature weights
 *    fnarray:  array of function values at corresponding quadrature positions
 */
KOKKOS_INLINE_FUNCTION
BS_REAL DoIntegration(const int n, const BS_REAL* wtarray,
                      const BS_REAL* fnarray)
{
    BS_REAL integral = 0.;

    for (int i = 0; i < n; ++i)
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
inline BS_REAL GaussLegendreIntegrateZeroInf(MyQuadrature* quad,
                                             MyFunction* func, BS_REAL t)
{

    BS_REAL f1_x[quad->nx], f2_x[quad->nx], f_y[quad->ny], f_z[quad->nz];
    BS_REAL w_y[quad->ny], w_z[quad->nz];
    // var3d var = var3d_default;
    BS_REAL var[3];

    for (int k = 0; k < quad->nz; ++k)
    {
        for (int j = 0; j < quad->ny; ++j)
        {
            for (int i = 0; i < quad->nx; ++i)
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
inline MyKernelQuantity
GaussLegendreIntegrateZeroInfSpecial(MyQuadrature* quad,
                                     MyFunctionSpecial* func, BS_REAL t)
{

    BS_REAL f1_em_e[quad->nx], f2_em_e[quad->nx];   // emissivity e neutrino
    BS_REAL f1_abs_e[quad->nx], f2_abs_e[quad->nx]; // absoptivity e neutrino
    BS_REAL f1_em_x[quad->nx], f2_em_x[quad->nx]; // emissivity mu/tau neutrino
    BS_REAL f1_abs_x[quad->nx],
        f2_abs_x[quad->nx]; // absorptivity mu/tau neutrino

    BS_REAL f_em_e_y[quad->ny];  // emissivity e neutrino
    BS_REAL f_abs_e_y[quad->ny]; // absoptivity e neutrino
    BS_REAL f_em_x_y[quad->ny];  // emissivity mu/tau neutrino
    BS_REAL f_abs_x_y[quad->ny]; // absorptivity mu/tau neutrino

    BS_REAL f_em_e_z[quad->ny];  // emissivity e neutrino
    BS_REAL f_abs_e_z[quad->ny]; // absoptivity e neutrino
    BS_REAL f_em_x_z[quad->ny];  // emissivity mu/tau neutrino
    BS_REAL f_abs_x_z[quad->ny]; // absorptivity mu/tau neutrino

    BS_REAL w_y[quad->ny], w_z[quad->nz];
    BS_REAL var[3];

    for (int k = 0; k < quad->nz; ++k)
    {
        for (int j = 0; j < quad->ny; ++j)
        {
            for (int i = 0; i < quad->nx; ++i)
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
inline BS_REAL GaussLaguerreIntegrateZeroInf(MyQuadrature* quad,
                                             MyFunction* func)
{

    constexpr BS_REAL zero = 0;

    BS_REAL f[quad->nx];

    if (quad->alpha == zero)
    {
        for (int i = 0; i < quad->nx; ++i)
        {
            f[i] = func->function(&quad->points[i],
                                  func->params); // * exp(quad->points[i]);
        }
    }
    else
    {
        for (int i = 0; i < quad->nx; ++i)
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
inline MyQuadratureIntegrand
GaussLegendreIntegrateFixedSplit2D(MyQuadrature* quad, MyFunctionMultiD* func,
                                   BS_REAL t)
{

    int num_integrands = func->my_quadrature_integrand.n;

    BS_REAL f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    BS_REAL f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    BS_REAL w_y[quad->ny];
    BS_REAL var[2];

    MyQuadratureIntegrand f1_vals, f2_vals;
    MyQuadratureIntegrand result;

    for (int j = 0; j < quad->ny; ++j)
    {

        var[1] = t * quad->points[quad->nx + j];

        for (int i = 0; i < quad->nx; ++i)
        {

            var[0]  = t * quad->points[i];
            f1_vals = func->function(var, func->params);

            var[0]  = t / quad->points[i];
            f2_vals = func->function(var, func->params);

            for (int k = 0; k < num_integrands; ++k)
            {
                f1_x[k][i] = f1_vals.integrand[k];
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }
        }

        for (int k = 0; k < num_integrands; ++k)
        {
            f1_y[k][j] = t * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                              DoIntegration(quad->nx, quad->w, f2_x[k]));
        }

        var[1] = t / quad->points[quad->nx + j];

        for (int i = 0; i < quad->nx; ++i)
        {

            var[0]  = t * quad->points[i];
            f1_vals = func->function(var, func->params);

            var[0]  = t / quad->points[i];
            f2_vals = func->function(var, func->params);

            for (int k = 0; k < num_integrands; ++k)
            {
                f1_x[k][i] = f1_vals.integrand[k];
                f2_x[k][i] =
                    f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
            }
        }


        for (int k = 0; k < num_integrands; ++k)
        {
            f2_y[k][j] = t * (DoIntegration(quad->nx, quad->w, f1_x[k]) +
                              DoIntegration(quad->nx, quad->w, f2_x[k]));
        }


        w_y[j] = quad->w[quad->nx + j];
    }

    for (int k = 0; k < num_integrands; ++k)
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
inline MyQuadratureIntegrand
GaussLegendreIntegrateFixedSplit1D(MyQuadrature* quad, MyFunctionMultiD* func,
                                   BS_REAL t)
{

    int num_integrands = func->my_quadrature_integrand.n;
    BS_REAL f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    BS_REAL var[2];

    MyQuadratureIntegrand f1_vals, f2_vals;

    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int i = 0; i < quad->nx; ++i)
    {

        var[0]  = t * quad->points[i];
        f1_vals = func->function(var, func->params);

        var[0]  = t / quad->points[i];
        f2_vals = func->function(var, func->params);

        for (int k = 0; k < num_integrands; ++k)
        {
            f1_x[k][i] = f1_vals.integrand[k];
            f2_x[k][i] =
                f2_vals.integrand[k] / (quad->points[i] * quad->points[i]);
        }
    }

    for (int k = 0; k < num_integrands; ++k)
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
inline MyQuadratureIntegrand GaussLegendreIntegrate2D(MyQuadrature* quad,
                                                      MyFunctionMultiD* func,
                                                      BS_REAL* tx, BS_REAL* ty)
{

    int num_integrands = func->my_quadrature_integrand.n;

    BS_REAL f1_x[num_integrands][quad->nx], f2_x[num_integrands][quad->nx];
    BS_REAL f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    BS_REAL w_y[quad->ny];
    BS_REAL var[2];

    MyQuadratureIntegrand result;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int j = 0; j < quad->ny; ++j)
        {

            var[1] = ty[k] * quad->points[quad->nx + j];

            for (int i = 0; i < quad->nx; ++i)
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

            for (int i = 0; i < quad->nx; ++i)
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
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand
GaussLegendreIntegrate1D(MyQuadrature* quad, MyFunctionMultiD* func, BS_REAL* t)
{

    int num_integrands = func->my_quadrature_integrand.n;
    BS_ASSERT(num_integrands <= num_max_integrands);
    BS_REAL f1_x[num_max_integrands][n_max], f2_x[num_max_integrands][n_max];
    BS_REAL var[2];
    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int i = 0; i < quad->nx; ++i)
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
MyQuadratureIntegrand
GaussLegendreIntegrate2DMatrix(const MyQuadrature* quad,
                               const M1MatrixKokkos2D* mat, BS_REAL t)
{
    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const BS_REAL t_sqr = t * t;

    BS_REAL w_i, w_j, w_ij;
    BS_REAL x_i, x_j;
    BS_REAL x2_i, x2_j, x2_ij;

    MyQuadratureIntegrand result = {0};

    for (int idx = 0; idx < total_num_species; ++idx)
    {

        for (int i = 0; i < n; ++i)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            for (int j = 0; j < n; ++j)
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

    for (int idx = 0; idx < num_integrands; ++idx)
    {
        result.integrand[idx] = t_sqr * result.integrand[idx];
    }

    return result;
}


KOKKOS_INLINE_FUNCTION
void GaussLegendreIntegrate2DMatrixForM1Coeffs(const MyQuadrature* quad,
                                               const M1MatrixKokkos2D* mat,
                                               BS_REAL t,
                                               MyQuadratureIntegrand* result_1,
                                               MyQuadratureIntegrand* result_2)
{
    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const BS_REAL t_sqr = t * t;

    BS_REAL w_i, w_j, w_ij;
    BS_REAL x_i, x_j;
    BS_REAL x2_i, x2_j, x2_ij;
    BS_REAL aux_1, aux_2, aux_3, aux_4;

    for (int idx = 0; idx < total_num_species; ++idx)
    {

        for (int i = 0; i < n; ++i)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            aux_1 = t * x_i;
            aux_2 = t / (x_i * x2_i);

            for (int j = 0; j < n; ++j)
            {

                x_j = quad->points[j];
                w_j = quad->w[j];

                w_ij  = w_i * w_j;
                x2_j  = x_j * x_j;
                x2_ij = x2_i * x2_j;

                aux_3 = aux_1 / x2_j;
                aux_4 = aux_2 / x2_j;

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
                    w_ij * (mat->m1_mat_em[idx][i][j] * aux_1 +
                            mat->m1_mat_em[idx][n + i][j] * aux_2 +
                            mat->m1_mat_em[idx][i][n + j] * aux_3 +
                            mat->m1_mat_em[idx][n + i][n + j] * aux_4);

                result_2->integrand[total_num_species + idx] +=
                    w_ij * (mat->m1_mat_ab[idx][i][j] * aux_1 +
                            mat->m1_mat_ab[idx][n + i][j] * aux_2 +
                            mat->m1_mat_ab[idx][i][n + j] * aux_3 +
                            mat->m1_mat_ab[idx][n + i][n + j] * aux_4);
            }
        }
    }

    for (int idx = 0; idx < num_integrands; ++idx)
    {
        result_1->integrand[idx] *= t_sqr;
        result_2->integrand[idx] *= t_sqr;
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void GaussLegendreIntegrate2DMatrixForNEPS(const MyQuadrature* quad,
                                           const M1MatrixKokkos2D* mat,
                                           BS_REAL t,
                                           MyQuadratureIntegrand* result_1,
                                           MyQuadratureIntegrand* result_2)
{
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL one  = 1;

    const int n              = quad->nx;
    const int num_integrands = 2 * total_num_species;

    const BS_REAL t_sqr = POW2(t);

    BS_REAL w_i, w_j, w_ij;
    BS_REAL x_i, x_j;
    BS_REAL x3_i;
    BS_REAL aux_1, aux_2, aux_3, aux_4;

    for (int idx = 0; idx < total_num_species; ++idx)
    {

        for (int i = 0; i < n; ++i)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x3_i = POW3(x_i);

            for (int j = 0; j < n; ++j)
            {

                x_j = quad->points[j];
                w_j = quad->w[j];

                w_ij = w_i * w_j;

                aux_1 = half * t * x_i * (one - x_j);
                aux_2 = half * t * x_i * (one + x_j);
                aux_3 = half * t * (one - x_j) / x_i;
                aux_4 = half * t * (one + x_j) / x_i;

                result_1->integrand[0 + idx] +=
                    w_ij * (x_i * (mat->m1_mat_em[idx][i][j] +
                                   mat->m1_mat_em[idx][i][n + j]) +
                            (mat->m1_mat_em[idx][n + i][j] +
                             mat->m1_mat_em[idx][n + i][n + j]) /
                                x3_i);

                result_1->integrand[total_num_species + idx] +=
                    w_ij * (x_i * (mat->m1_mat_ab[idx][i][j] +
                                   mat->m1_mat_ab[idx][i][n + j]) +
                            (mat->m1_mat_ab[idx][n + i][j] +
                             mat->m1_mat_ab[idx][n + i][n + j]) /
                                x3_i);

                result_2->integrand[0 + idx] +=
                    w_ij * (x_i * (aux_1 * mat->m1_mat_em[idx][i][j] +
                                   aux_2 * mat->m1_mat_em[idx][i][n + j]) +
                            (aux_3 * mat->m1_mat_em[idx][n + i][j] +
                             aux_4 * mat->m1_mat_em[idx][n + i][n + j]) /
                                x3_i);

                result_2->integrand[total_num_species + idx] +=
                    w_ij * (x_i * (aux_1 * mat->m1_mat_ab[idx][i][j] +
                                   aux_2 * mat->m1_mat_ab[idx][i][n + j]) +
                            (aux_3 * mat->m1_mat_ab[idx][n + i][j] +
                             aux_4 * mat->m1_mat_ab[idx][n + i][n + j]) /
                                x3_i);
            }
        }
    }

    for (int idx = 0; idx < num_integrands; ++idx)
    {
        result_1->integrand[idx] *= half * t_sqr;
        result_2->integrand[idx] *= half * t_sqr;
    }

    return;
}

KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand GaussLegendreIntegrate1DMatrix(const MyQuadrature* quad,
                                                     const int num_integrands,
                                                     const BS_REAL mat[][n_max],
                                                     BS_REAL* t)
{
    const int n = quad->nx;
    BS_REAL x_i, w_i, x2_i;

    MyQuadratureIntegrand result = {0};

    for (int idx = 0; idx < num_integrands; ++idx)
    {

        for (int i = 0; i < n; ++i)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = x_i * x_i;

            result.integrand[idx] +=
                w_i * (mat[idx][i] + mat[idx][n + i] / x2_i);
        }
    }

    for (int idx = 0; idx < num_integrands; ++idx)
    {
        result.integrand[idx] *= t[idx];
    }

    return result;
}


KOKKOS_INLINE_FUNCTION
void GaussLegendreIntegrate1DMatrixOnlyNumber(const MyQuadrature* quad,
                                              const int num_integrands,
                                              const BS_REAL mat[][n_max],
                                              BS_REAL* t,
                                              MyQuadratureIntegrand* out_n,
                                              MyQuadratureIntegrand* out_j)
{
    const int n = quad->nx;
    BS_REAL x_i, w_i, x2_i;

    for (int idx = 0; idx < num_integrands; ++idx)
    {
        for (int i = 0; i < n; ++i)
        {

            x_i  = quad->points[i];
            w_i  = quad->w[i];
            x2_i = POW2(x_i);

            out_n->integrand[idx] +=
                w_i * (mat[idx][i] + mat[idx][n + i] / x2_i);

            BS_ASSERT(isfinite(out_n->integrand[idx]) &&
                      out_n->integrand[idx] >= 0.);

            out_j->integrand[idx] +=
                w_i * (mat[idx][i] * t[idx] * x_i +
                       mat[idx][n + i] * (t[idx] / x_i) / x2_i);

            BS_ASSERT(isfinite(out_j->integrand[idx]) &&
                      out_j->integrand[idx] >= 0.);
        }

        out_n->integrand[idx] *= t[idx];
        out_j->integrand[idx] *= t[idx];
    }

    return;
}

/* Perform 1d integration (finite interval) of multiple functions using a
 * Gauss-Legendre quadrature
 *
 * quad:    must be a properly populated 1d quadrature generated from between
 * interval of integration func:    the function(s) to be integrated
 */
inline MyQuadratureIntegrand
GaussLegendreIntegrate1DFiniteInterval(MyQuadrature* quad,
                                       MyFunctionMultiD* func, BS_REAL* t)
{
    (void)t;

    int num_integrands = func->my_quadrature_integrand.n;
    BS_REAL f_x[num_integrands][quad->nx];
    BS_REAL var[2];
    MyQuadratureIntegrand result;

    result.n = num_integrands;

    for (int k = 0; k < num_integrands; ++k)
    {
        for (int i = 0; i < quad->nx; ++i)
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
inline MyQuadratureIntegrand GaussLegendreIntegrate2DFiniteInterval(
    MyQuadrature* quad, MyFunctionMultiD* func, BS_REAL* tx, BS_REAL* ty)
{
    (void)tx;

    int num_integrands = func->my_quadrature_integrand.n;

    BS_REAL f1_x[num_integrands][quad->nx];
    BS_REAL f1_y[num_integrands][quad->ny], f2_y[num_integrands][quad->ny];

    BS_REAL w_y[quad->ny];
    BS_REAL var[2];

    MyQuadratureIntegrand result;

    for (int k = 0; k < num_integrands; ++k)
    {

        for (int j = 0; j < quad->ny; ++j)
        {

            var[1] = ty[k] * quad->points[quad->nx + j];
            for (int i = 0; i < quad->nx; ++i)
            {
                var[0] = quad->points[i];
                MyQuadratureIntegrand f1_vals =
                    func->function(var, func->params);
                f1_x[k][i] = f1_vals.integrand[k];
            }

            f1_y[k][j] = DoIntegration(quad->nx, quad->w, f1_x[k]);

            var[1] = ty[k] / quad->points[quad->nx + j];
            for (int i = 0; i < quad->nx; ++i)
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
