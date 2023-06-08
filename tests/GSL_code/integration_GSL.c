#include <gsl/gsl_integration.h>

#include "../tests.h"

// GSL implementation of numerical fixed-point quadrature integration routines
// Based on GSL 2.7.1 (https://www.gnu.org/software/gsl/doc/html/integration.html)


// Fixed-point quadrature integration with generic weight w(x)
// For integrals written as \int_{a}^{b} f(x)*w(x) dx

/*
 * n : Number of quadrature points
 * PairF : Function to be integrated (GSL variable: struct with func + void pointer to params)
 * gauss_type: type of weighting function
*/

double GSL_integration(const int n,
                       const double a, const double b,
		       const double alpha, const double beta,
		       const gsl_function *f,
		       const gsl_integration_fixed_type * gauss_type) {
  gsl_integration_fixed_workspace *w;
  double result;
  w = gsl_integration_fixed_alloc(gauss_type, n, a, b, alpha, beta);
  gsl_integration_fixed(f, &result, w);
  gsl_integration_fixed_free(w);
  return result;
}


// Gauss-Laguerre integration (a,+inf): gauss_type -> gsl_integration_fixed_laguerre
// w(x) = (x-a)^alpha * exp(-b*(x-a))
// Constraints: alpha > -1, b > 0 (beta is ignored)
double GSLLagQuadrature(const int n,
                        const double a, const double b,
                        const double alpha,
                        const gsl_function *f) {
  const gsl_integration_fixed_type * lag_type = gsl_integration_fixed_laguerre;	
  return GSL_integration(n, a, b, alpha, 0., f, lag_type);
}


// Gauss-Legendre integration (a,b): gauss_type -> gsl_integration_fixed_legendre
// w(x) = 1
// Constraints: b > a (alpha and beta are ignored)
double GSLLegQuadrature(const int n,
                        const double a, const double b,
                        const gsl_function *f) {
  const gsl_integration_fixed_type * leg_type = gsl_integration_fixed_legendre;	
  return GSL_integration(n, a, b, 0., 0., f, leg_type);
}

// Specific GSL implementation of Gauss-Legendre integration
double GSLLegIntegrate(const int n,
                       const double a, const double b,
                       const gsl_function *f) {
  gsl_integration_glfixed_table * t;
  t = gsl_integration_glfixed_table_alloc(n);
  double result = gsl_integration_glfixed(f, a, b, t);
  gsl_integration_glfixed_table_free(t);
  return result;
}


// Computation of (0,+inf) integral using two-piece Gauss-Legendre integration
// Split at s in two integrals and map them into (0,1) integration thorugh a change of variables
double GslLegInfSplit(const int n, const gsl_function *f, const double s) {
  double xi, wi;
  double result = 0.;
  gsl_integration_glfixed_table * t;
  t = gsl_integration_glfixed_table_alloc(n);
  for (int i=0;i<n;i++) {
    gsl_integration_glfixed_point(0, 1, i, &xi, &wi, t);
    result += wi * f->function(s*xi, f->params);
    result += wi * f->function(s/xi, f->params) / (xi*xi);
  }
  return s * result;
} 

