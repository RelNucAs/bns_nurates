#include <math.h>

#include "functions.h"

// TODO: change names of variables and functions following Google C style

// Ansatz for neutrino distribution function in optically thick and thin regimes, respectively
// TODO: Eventually move this somewhere else (f_nu reconstruction section)
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

// TODO: do some tests and optimize the following parameters
const int ntrial = 1000;    // Max number of NR iterations
const double tolx = 1.e-5;  // Accuracy level on the variable
const double tolf = 1.e-7;  // Accuracy level on the functions

void MNewt2d(double *x,
             double C[2],
             void (*fdf)(double*, double*, double*, double*)) {
  const int n = 2; // Matrix dimension
  int i;
  double p[n], fvec[n];
  double fjac[n*n]; // Jacobian matrix
  double finv[n*n]; // Inverse of Jacobian matrix
  for (int k=0;k<ntrial;k++) {
    fdf(x,C,fvec,fjac);
    double errf=0.0;
    for (i=0;i<n;i++) errf += fabs(fvec[i]);
    if (errf <= tolf) return;
    Invert2DMat(fjac,finv);
    for (i=0;i<n;i++) p[i] = -(finv[i*n]*fvec[0] + finv[i*n+1]*fvec[1]); //-fvec[i];
    double errx=0.0;
    for (i=0;i<n;i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
      //printf("x[%d] = %.3e\n", i, x[i]);
    }
    if (errx <= tolx) return;
  }
  return;
}
