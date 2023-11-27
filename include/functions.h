//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
#define BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>

#include "bns_nurates.h"


// Exception handling from Numerical Recipes
// @TODO: decide how to handle errors in the code

#ifndef _USENRERRORCLASS_
#define throw(message) { \
  printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); \
  exit(1); \
}
#else
struct NRerror {
  char *message;
  char *file;
  int line;
  NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
  printf("ERROR: %s\n     in file %s at line %d\n",
          err.message, err.file, err.line);
  exit(1);
}
#endif


/*===========================================================================*/

// digamma.c

// Evaluation of PairPsi (Digamma) function
double SFPsi(const double x);

/*===========================================================================*/

// fermi_integrals.c

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -9/2 */
double FDI_m92(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -7/2 */
double FDI_m72(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k=-5/2 */
double FDI_m52(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -3/2 */
double FDI_m32(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -1/2 */
double FDI_m12(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 0 */
double FDI_0(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1/2 */
double FDI_p12(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1 */
double FDI_p1(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3/2 */
double FDI_p32(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 2 */
double FDI_p2(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5/2 */
double FDI_p52(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3 */
double FDI_p3(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7/2 */
double FDI_p72(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 4 */
double FDI_p4(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9/2 */
double FDI_p92(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5 */
double FDI_p5(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 11/2 */
double FDI_p112(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 6 */
double FDI_p6(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 13/2 */
double FDI_p132(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7 */
double FDI_p7(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 15/2 */
double FDI_p152(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 8 */
double FDI_p8(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 17/2 */
double FDI_p172(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9 */
double FDI_p9(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 19/2 */
double FDI_p192(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 10 */
double FDI_p10(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 21/2 */
double FDI_p212(const double x);

/*===========================================================================*/

// fermi_distr.c

// Computation of Fermi-Dirac distribution function
double FermiDistr(const double e, const double temp, const double mu);

/*===========================================================================*/

// gamma.c

// Computation of gamma function
double Gammln(const double xx);

// Computation of gamma function using Stirling's approximation
double GammaStirling(const double x);

/*===========================================================================*/

// lambert.c

// Recursive formula for W0 real branch of Lambert function
double W0(double x);

/*===========================================================================*/

// mnewt.c

// Implementation of 1D Newton-Raphson root finding
double MNewt1d(double guess, double x1, double x2, double f0,
               MyFunction *func, MyFunction *dfunc);

// Implementation of 2D Newton-Raphson root finding
void MNewt2d(double *x,
             double C[2],
             void (*fdf)(double *, double *, double *, double *));

/*===========================================================================*/

// safe_exp.c

// Safe exp function to avoid underflow/overflow
double SafeExp(const double x);

/*===========================================================================*/

// theta.c

// Step function implementation with if statement
double HeavisidePiecewise(const double x);

// Step function approximation with tanh - (Eq.5)
double HeavisideTanhApprox(const double x);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
