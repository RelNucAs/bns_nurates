//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_TESTS_TESTS_H_
#define BNS_NURATES_TESTS_TESTS_H_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

/*===========================================================================*/

// GSL_code/integration_GSL.c

// Gauss-Laguerre integration
double GSL_lag_quadr(const int n,
                     const double a, const double b,
                     const double alpha, const double beta,
                     const gsl_function *f);


// Gauss-Legendre integration
/* Using general fixed-quadrature framework */
double GSL_leg_quadr(const int n,
                     const double a, const double b,
                     const double alpha, const double beta,
                     const gsl_function *f);


/* Specific implementation */
double GSL_leg_integ(const int n,
                     const double a, const double b,
                     const gsl_function *f);


// Computation of (0,+inf) integral using two-piece Gauss-Legendre integration
double GSL_leg_split(const int n, const gsl_function *f, const double s);

/*===========================================================================*/

// GSL_code/mnewt_GSL.c

/* One-dimensional root finding functions */

// 1D Newton-Raphson with analytic derivative 
double MNewt1d_GSL(const double guess, gsl_function_fdf *fdf);


/* Two-dimensional root finding functions */

// 2D Newton-Raphson with analytic Jacobian
int MNewt2d_GSL(double *xi, double *xf, gsl_multiroot_function_fdf *fdf);

// 2D Newton-Raphson with finite difference Jacobian
int MNewt2d_fd_GSL(double *xi, double *xf, gsl_multiroot_function * f);

#endif //BNS_NURATES_TESTS_TESTS_H_


