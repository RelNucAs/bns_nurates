#pragma once

void usrfun(VecDoub_I &x, VecDoub_O &fvec, MatDoub_O &fjac);

void mnewt(const Int ntrial, VecDoub_IO &x, const Doub tolx, const Doub tolf);

void thick_f(double *x, double *C, double *f);

void thick_df(double *x, double *C, double *J);

void thin_f(double *x, double *C, double *f);

void thin_df(double *x, double *C, double *f);

