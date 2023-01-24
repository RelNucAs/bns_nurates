#pragma once
#include "tools/nr3.h"
#include "spectral_function.hpp"

void save_gauleg(int n, const Doub x1, const Doub x2);

void save_gaulag(int n, const Doub alf);

void read_gleg_wgts(int n, double *x, double *w);

void read_glag_wgts(int n, double *x, double *w, double alpha);

double do_integration(int n, double *wgt, double *spectrArray);

double glag_integ(int n, my_function F, double alpha);

double split_value(double T, double eta_e);

double gleg_integ(int n, my_function F, double t);

double find_max(int nslice, double x0, double x1, double (*func)(double, struct my_f_params), struct my_f_params p);

double find_guess(int nslice, my_function F);

