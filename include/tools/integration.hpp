#ifndef INTEGRATION_H
#define INTEGRATION_H

// \file  integration.hpp
// \brief Quadrature methods (Gauss-Legendre and Gauss-Laguerre) for intergration

#include "tools/nr3.h"
#include "spectral_function.hpp"


/* Save abscissas and weights to txt file */
// Gauss-Legendre
void save_gauleg(int n, const Doub x1, const Doub x2);
// Gauss-Laguerre
void save_gaulag(int n, const Doub alf);

/* Read abscissas and weights from txt file */
// Gauss-Legendre
void read_gleg_wgts(int n, double *x, double *w);
// Gauss-Laguerre
void read_glag_wgts(int n, double *x, double *w, double alpha);

/* Make convolution between sampled function and weights */
double do_integration(int n, double *wgt, double *spectrArray);

/* Perform Gauss-Legendere integration */
double glag_integ(int n, my_function F, double alpha);

/* Perform Gauss-Laguerre integration */
double gleg_integ(int n, my_function F, double t);

/* Estimate maximum of the function to be integrated */
double find_max(int nslice, double x0, double x1, double (*func)(double, struct my_f_params), struct my_f_params p);


double find_guess(int nslice, my_function F);

//double split_value(double T, double eta_e);

#endif
