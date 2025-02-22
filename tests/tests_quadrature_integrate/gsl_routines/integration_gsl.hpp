//
// Created by maitraya on 7/20/23.
//

#ifndef BNS_NURATES_TESTS_TESTS_QUADRATURE_INTEGRATE_GSL_CODE_INTEGRATION_GSL_H_
#define BNS_NURATES_TESTS_TESTS_QUADRATURE_INTEGRATE_GSL_CODE_INTEGRATION_GSL_H_

#ifdef GSL_INCLUDES_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

/*===========================================================================*/

// gsl_routines/integration_GSL.c

// Gauss-Laguerre integration
double GSLLagQuadrature(const int n, const double a, const double b,
                        const double alpha, const gsl_function* f);


// Gauss-Legendre integration
/* Using general fixed-quadrature framework */
double GSLLegQuadrature(const int n, const double a, const double b,
                        const gsl_function* f);

/* Specific implementation */
double GSLLegIntegrate(const int n, const double a, const double b,
                       const gsl_function* f);

// Computation of (0,+inf) integral using two-piece Gauss-Legendre integration
double GslLegInfSplit(const int n, const gsl_function* f, const double s);

/*===========================================================================*/

// gsl_routines/mnewt_GSL.c

/* One-dimensional root finding functions */

// 1D Newton-Raphson with analytic derivative
double MNewt1dGSL(const double guess, gsl_function_fdf* fdf);


/* Two-dimensional root finding functions */

// 2D Newton-Raphson with analytic Jacobian
int MNewt2dGSL(double* xi, double* xf, gsl_multiroot_function_fdf* fdf);

// 2D Newton-Raphson with finite difference Jacobian
int MNewt2dGSL_fd(double* xi, double* xf, gsl_multiroot_function* f);

// ===================================================================

#endif // GSL_INCLUDES_H_

#endif // BNS_NURATES_TESTS_TESTS_QUADRATURE_INTEGRATE_GSL_CODE_INTEGRATION_GSL_H_
