#pragma once
#include <functional>
#include "tools/gamma.h"
#include "tools/digamma.h"

void thick_fdf(const double x[2], const double C[2], double *fvec, double *fjac) {
	const double n_nu = C[0];
	const double J_nu = C[1];
	
	fvec[0] = n_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],3.);
	fvec[1] = J_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p3(x[1])/pow(x[0],4.);

	fjac[0] =  3. * 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],4.);
        fjac[1] = -2. * 4.*pi*Fermi_integral_p1(x[1])/pow(x[0],3.);
        fjac[2] =  4. * 4.*pi*Fermi_integral_p3(x[1])/pow(x[0],5.);
        fjac[3] = -3. * 4.*pi*Fermi_integral_p2(x[1])/pow(x[0],4.);
        return;
}

void thin_fdf(const double x[2], const double C[2], double *fvec, double *fjac) {
	const double n_nu = C[0];
	const double J_nu = C[1];

        fvec[0] = n_nu*pow(h*c,3.) - 4.*pi*exp(gammln(x[0]+3.))/pow(x[1],x[0]+3.);
        fvec[1] = J_nu*pow(h*c,3.) - 4.*pi*exp(gammln(x[0]+4.))/pow(x[1],x[0]+4.);

        fjac[0] = 4.*pi*exp(gammln(x[0]+3.))/pow(x[1],x[0]+3.)*(log(x[1])-sf_psi(x[0]+3.));
        fjac[1] = 4.*pi*(x[0]+3.)*exp(gammln(x[0]+3.))/pow(x[1],x[0]+4.);
        fjac[2] = 4.*pi*exp(gammln(x[0]+4.))/pow(x[1],x[0]+4.)*(log(x[1])-sf_psi(x[0]+4.));
        fjac[3] = 4.*pi*(x[0]+4.)*exp(gammln(x[0]+4.))/pow(x[1],x[0]+5.);
        //fjac[0] = 4.*pi*exp(gammln(x[0]+3.)/pow(x[1],x[0]+3.)*(log(x[1])-gsl_sf_psi(x[0]+3.));
        //fjac[1] = 4.*pi*(x[0]+3.)*exp(gammln(x[0]+3.)/pow(x[1],x[0]+4.);
        //fjac[2] = 4.*pi*exp(gammln(x[0]+4.)/pow(x[1],x[0]+4.)*(log(x[1])-gsl_sf_psi(x[0]+4.));
        //fjac[3] = 4.*pi*(x[0]+4.)*exp(gammln(x[0]+4.)/pow(x[1],x[0]+5.);
        return;
}

void invert_2Dmat(const double in[2*2], double *out) {
	const double den = in[0]*in[3] - in[1]*in[2];
	out[0] =  in[3] / den;
	out[1] = -in[1] / den;
	out[2] = -in[2] / den;
	out[3] =  in[0] / den;
	return;
}

void mnewt2d(double *x, const double C[2], std::function<void(double*&, const double*&, double [2], double [4])> f_df) {
	//size_t x_size = sizeof(x) / sizeof(double); //equal to one since pointer points to first element of the array
	//if (x_size!=2) {
	//	printf("Wrong dimension (!=2) of input array c in 'mnewt2d'\n");
	//	return;
	//}
	const int ntrial = 1000, n = 2;
	int i;
	const double tolx = 1.e-5;
	const double tolf = 1.e-7;
	double p[n], fvec[n];
	double fjac[n*n], finv[n*n];
	for (int k=0;k<ntrial;k++) {
		f_df(x,C,fvec,fjac);
		double errf=0.0;
		for (i=0;i<n;i++) errf += abs(fvec[i]);
		if (errf <= tolf) return;
		invert_2Dmat(fjac,finv);
		for (i=0;i<n;i++) p[i] = -(finv[i*n]*fvec[0] + finv[i*n+1]*fvec[1]); //-fvec[i];
		double errx=0.0;
		for (i=0;i<n;i++) {
			errx += abs(p[i]);
			x[i] += p[i];
			//printf("x[%d] = %.3lf\n", i, x[i]);
		}
		if (errx <= tolx) return;
	}
	return;
}
