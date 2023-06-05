//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_TESTS_TESTS_H_
#define BNS_NURATES_TESTS_TESTS_H_

#include <gsl/gsl_integration.h>

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

#endif //BNS_NURATES_TESTS_TESTS_H_


