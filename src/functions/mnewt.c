#include <math.h>

#include "functions.h"


// Ansatz for neutrino distribution function in optically thick and thin regimes, respectively
// @TODO: Eventually move this somewhere else (f_nu reconstruction section)
// 	 and understand how to implement computation of Fermi integrals in the code
//
//void thick_fdf(double *x, double *C, double *fvec, double *fjac) {
//	const double n_nu = C[0];
//	const double J_nu = C[1];
//	
//	fvec[0] = n_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],3.);
//	fvec[1] = J_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p3(x[1])/pow(x[0],4.);
//
//	fjac[0] =  3. * 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],4.);
//      fjac[1] = -2. * 4.*pi*Fermi_integral_p1(x[1])/pow(x[0],3.);
//      fjac[2] =  4. * 4.*pi*Fermi_integral_p3(x[1])/pow(x[0],5.);
//      fjac[3] = -3. * 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],4.);
//      return;
//}

//void thin_fdf(double *x, double *C, double *fvec, double *fjac) {
//  const double n_nu = C[0];
//  const double J_nu = C[1];
//  fvec[0] = n_nu*pow(h*c,3.) - 4.*pi*exp(Gammln(x[0]+3.))/pow(x[1],x[0]+3.);
//  fvec[1] = J_nu*pow(h*c,3.) - 4.*pi*exp(Gammln(x[0]+4.))/pow(x[1],x[0]+4.);
  
//  fjac[0] = 4.*pi*exp(Gammln(x[0]+3.))/pow(x[1],x[0]+3.)*(log(x[1])-SFPsi(x[0]+3.));
//  fjac[1] = 4.*pi*(x[0]+3.)*exp(Gammln(x[0]+3.))/pow(x[1],x[0]+4.);
//  fjac[2] = 4.*pi*exp(Gammln(x[0]+4.))/pow(x[1],x[0]+4.)*(log(x[1])-SFPsi(x[0]+4.));
//  fjac[3] = 4.*pi*(x[0]+4.)*exp(Gammln(x[0]+4.))/pow(x[1],x[0]+5.);
//  //fjac[0] = 4.*pi*exp(gammln(x[0]+3.)/pow(x[1],x[0]+3.)*(log(x[1])-gsl_sf_psi(x[0]+3.));
//  //fjac[1] = 4.*pi*(x[0]+3.)*exp(gammln(x[0]+3.)/pow(x[1],x[0]+4.);
//  //fjac[2] = 4.*pi*exp(gammln(x[0]+4.)/pow(x[1],x[0]+4.)*(log(x[1])-gsl_sf_psi(x[0]+4.));
//  //fjac[3] = 4.*pi*(x[0]+4.)*exp(gammln(x[0]+4.)/pow(x[1],x[0]+5.);
//  return;
//}

/*===========================================================================*/

// One-dimensional root finding functions

// 1D root-finding parameters
// @TODO: find suitable values for the following parameters
const int  ntrial_1d = 150;      // Maximum allowed number of iterations.
const double xacc_1d = 1.0E-07;  // Set the accuracy for Newton Raphson

// 1D Newton-Raphson with analytic derivative 
double MNewt1d(double guess, double x1, double x2, double f0,
               MyFunction *func, MyFunction *dfunc) {
// Using a combination of Newton-Raphson and bisection, return the solution
// of f(x) = f0 bracketed between x1 and x2. The root will be refined
// until its accuracy is known within xacc. f is a user-supplied struct
// that returns the function value at the point x. df is a user-supplied
// struct that returns the value of the function first derivative at the
// point x.
  double xh, xl;
  double fl = func->function(&x1, func->params) - f0;
  double fh = func->function(&x2, func->params) - f0;
  
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    printf("xl = %.3e, fl = %.3e\n", x1, fl);
    printf("xh = %.3e, fh = %.3e\n", x2, fh);
    throw("Root must be bracketed in rtsafe");
  }

  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  
  if (fl < 0.0) { //Orient the search so that f(xl) < 0.
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }

  double rts=guess; //0.5*(x1+x2);  // Initialize the guess for root,
  double dxold=fabs(x2-x1);         // the “stepsize before last,”
  double dx=dxold;                  // and the last step.

  double  f =  func->function(&rts,  func->params) - f0;
  double df = dfunc->function(&rts, dfunc->params);
  
  for (int j=0;j<ntrial_1d;j++) { //Loop over allowed iterations.
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) { //Bisect if Newton out of range, or not decreasing fast enough.
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else { //Change in root is negligible. Newton step acceptable. Take it.
      dxold=dx;
      dx=f/df;
      double temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }

    if (fabs(dx) < xacc_1d) return rts; //Convergence criterion.

    f =  func->function(&rts,  func->params) - f0; // The one new function evaluation per iteration.
    df = dfunc->function(&rts, dfunc->params);

    if (f < 0.0) { // Maintain the bracket on the root.
       xl=rts;
    } else {
      xh=rts;
    }
  }
  throw("Maximum number of iterations exceeded in rtsafe");
}

/*===========================================================================*/

// Two-dimensional root finding functions

// Invert two-dimensional matrix
/*
 * in  : input 2x2 matrix
 * out : output 2x2 matrix
 *
 * Entries: 
 * 	0 --> [0,0]
 * 	1 --> [0,1]
 *	2 --> [1,0]
 *	3 --> [1,1]
*/
void Invert2DMat(double *in, double *out) {
  const double den = in[0]*in[3] - in[1]*in[2];
  out[0] =  in[3] / den;
  out[1] = -in[1] / den;
  out[2] = -in[2] / den;
  out[3] =  in[0] / den;
  return;
}

// Implementation of 2D Newton-Raphson root finding

/*
 * Finding solution of A = C - F[x] = 0 (system of two non-linear coupled equations) 
 *
 * Input:
 * 	- x (dim = 2) --> initial guess
 * 	- C (dim = 2) --> inhomogeneous term
 *	- PairF (dim = 2) --> calculated using pointer to fdf(x,C,fvec,fjac) function
 * Output:
 * 	- x (1x2 matrix) --> Newton-Raphson solution after convergence
 *
 * - fvec: A (dim = 2)
 * - fjac: Jacobian matrix of A (dim = 2x2)
*/

// 2D root-finding parameters
// @TODO: do some tests and optimize the following parameters
const int  ntrial_2d = 1000;    // Max number of NR iterations
const double tolx_2d = 1.e-5;  // Accuracy level on the variable
const double tolf_2d = 1.e-7;  // Accuracy level on the functions

void MNewt2d(double *x,
             double C[2],
             void (*fdf)(double*, double*, double*, double*)) {
  const int n = 2; // Matrix dimension
  int i;
  double p[n], fvec[n];
  double fjac[n*n]; // Jacobian matrix
  double finv[n*n]; // Inverse of Jacobian matrix
  for (int k=0;k<ntrial_2d;k++) {
    fdf(x,C,fvec,fjac);
    double errf=0.0;
    for (i=0;i<n;i++) errf += fabs(fvec[i]);
    if (errf <= tolf_2d) return;
    Invert2DMat(fjac,finv);
    for (i=0;i<n;i++) p[i] = -(finv[i*n]*fvec[0] + finv[i*n+1]*fvec[1]); //-fvec[i];
    double errx=0.0;
    for (i=0;i<n;i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
      //printf("x[%d] = %.3e\n", i, x[i]);
    }
    if (errx <= tolx_2d) return;
  }
  return;
}

/*===========================================================================*/