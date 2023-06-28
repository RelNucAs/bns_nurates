//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
#define BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>

#include "../bns_nurates.h"


// Exception handling from Numerical Recipes
// @TODO: decide how to handling errors in the code

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

// fermi.c

// Computation of Fermi-Dirac distribution function
double FermiDistr(const double e, const double temp, const double mu);

/*===========================================================================*/

// gamma.c

// Computation of gamma function
double Gammln(const double xx);

/*===========================================================================*/

// mnewt.c

// Implementation of 1D Newton-Raphson root finding
double MNewt1d(double x, double guess,
             double x1, double x2, double f0,
             MyFunction *func, MyFunction *dfunc);

// Implementation of 2D Newton-Raphson root finding
void MNewt2d(double *x,
             double C[2],
             void (*fdf)(double*, double*, double*, double*));

/*===========================================================================*/

// safe_exp.c

// Safe exp function to avoid underflow/overflow
double SafeExp(const double x);

/*===========================================================================*/

#endif //BNS_NURATES_SRC_FUNCTIONS_FUNCTIONS_H_
