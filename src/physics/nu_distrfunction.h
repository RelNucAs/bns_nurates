#pragma once

#include "constants.hpp"
#include "tools/fermi_integrals.hpp"
#include "tools/gamma.hpp"
#include "tools/digamma.hpp"

void mnewt(const Int ntrial, VecDoub_IO &x, const Doub tolx, const Doub tolf) {
        Int i,n=x.size();
        VecDoub p(n),fvec(n);
        MatDoub fjac(n,n);
        for (Int k=0;k<ntrial;k++) {
                usrfun(x,fvec,fjac);
                Doub errf=0.0;
                for (i=0;i<n;i++) errf += abs(fvec[i]);
                if (errf <= tolf) return;
                for (i=0;i<n;i++) p[i] = -fvec[i];
                LUdcmp alu(fjac);
                alu.solve(p,p);
                Doub errx=0.0;
                for (i=0;i<n;i++) {
                        errx += abs(p[i]);
                        x[i] += p[i];
                }
                if (errx <= tolx) return;
        }
        return;
}


void thick_f(double *x, double *C, double *f) {
	double A = x[0];
	double B = x[1];
	f[0] = C[0]*pow(h*c,3.) - 4.*pi*Fermi_integral_p2(B)/pow(A,3.);
	f[1] = C[1]*pow(h*c,3.) - 4.*pi*Fermi_integral_p3(B)/pow(A,4.);
        return;
}


void thick_df(double *x, double *C, double *J) {
	double A = x[0];
	double B = x[1];
        J[0] =  3. * 4.*pi*Fermi_integral_p2(B)/pow(A,4.);
        J[1] = -2. * 4.*pi*Fermi_integral_p1(B)/pow(A,3.);
        J[2] =  4. * 4.*pi*Fermi_integral_p3(B)/pow(A,5.);
        J[3] = -3. * 4.*pi*Fermi_integral_p2(B)/pow(A,4.);

        return;
}

void thin_f(double *x, double *C, double *f) {
	double A = x[0];
	double B = x[1];
        f[0] = C[0]*pow(h*c,3.) - 4.*pi*gammln(A+3.)/pow(B,A+3.);
        f[0] = C[0]*pow(h*c,3.) - 4.*pi*gammln(A+3.)/pow(B,A+3.);
        f[1] = C[1]*pow(h*c,3.) - 4.*pi*gammln(A+4.)/pow(B,A+4.);
        return;
}

void thin_df(double *x, double *C, double *f) {
	double A = x[0];
	double B = x[1];
        J[0] = 4.*pi*gammln(A+3.)/pow(B,A+3.)*(log(B)-sf_psi(A+3.)); //4.*pi*gsl_sf_gamma(A+3.)/pow(B,A+3.)*(log(B)-gsl_sf_psi(A+3.));
        J[1] = 4.*pi*(A+3.)*gammln(A+3.)/pow(B,A+4.); //4.*pi*(A+3.)*gsl_sf_gamma(A+3.)/pow(B,A+4.);
        J[2] = 4.*pi*gammln(A+4.)/pow(B,A+4.)*(log(B)-sf_psi(A+4.)); //4.*pi*gsl_sf_gamma(A+4.)/pow(B,A+4.)*(log(B)-gsl_sf_psi(A+4.));
        J[3] = 4.*pi*(A+4.)*gammln(A+4.)/pow(B,A+5.); //4.*pi*(A+4.)*gsl_sf_gamma(A+4.)/pow(B,A+5.);
        return;
}

