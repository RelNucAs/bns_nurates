#pragma once
#include "../constants_bis.h"
//#include "../tools/nr3.h"
#include "complete_FG.h"

using namespace constants;

struct n_net_FD {double n_net; double dn_net; double f12; double f32; double a_f12; double a_f32; };

//.......Given the density rho, the non relativistic degeneracy eta, the temperature T
//.......and the particle species, the subroutine computes the net fraction ynet and the first eta derivative.
//.......UNITS: rho in g/cm**3; T in MeV
double n_net_f(const double eta, const double T, const double mLep) {
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        const double theta = T/mLep;
        double f1, f2, f3, f4;

        f1 = compute_res( eta,          theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 0.5);
        f3 = compute_res( eta,          theta, 1.5);
        f4 = compute_res(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}


double n_net_df(const double eta, const double T, const double mLep) {
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        const double theta = T/mLep;
        double f1, f2, f3, f4;

        f1 = compute_res_ed( eta,          theta, 0.5);
        f2 = compute_res_ed(-eta-2./theta, theta, 0.5);
        f3 = compute_res_ed( eta,          theta, 1.5);
        f4 = compute_res_ed(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}

void n_net_fdf(const double eta, const double T, const double mLep, struct n_net_FD &FD) {
	const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.);
	const double theta = T/mLep;
	double f1, f2, f3, f4;

	FD.f12   = compute_res( eta,          theta, 0.5);
	FD.f32   = compute_res( eta,          theta, 1.5);
	FD.a_f12 = compute_res(-eta-2./theta, theta, 0.5);
	FD.a_f32 = compute_res(-eta-2./theta, theta, 1.5);

        f1 = compute_res_ed( eta,          theta, 0.5);
        f2 = compute_res_ed(-eta-2./theta, theta, 0.5);
        f3 = compute_res_ed( eta,          theta, 1.5);
        f4 = compute_res_ed(-eta-2./theta, theta, 1.5);

	FD.n_net  = K*pow(theta,1.5)*(FD.f12-FD.a_f12+theta*(FD.f32-FD.a_f32));
	FD.dn_net = K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
	return;
}
	

template<int species>
double find_guess_eta(double nLep, double T) {
	double nq;

	if constexpr(species == 1) {
		const double t = T/kB; //temp in K
		const double nfm = ne*1.e-39; //number density in fm^{-3}
		const double pF = h*c*pow(3.*ne/(8.*pi),1./3.); //[MeV]

	//.........Find the guess
		if ((t<=1.e10) && (nfm <= 1.e-17)) { //classical NR
			nq = pow(2.*pi*me*T/(h*h*c*c),1.5); //cm^{-3}
			return -log(fabs(2.*nq/ne));
		} else if ((t<=1.e10) && (nfm>1.e-17) && (nfm <= 1.e-10)) { //degenerate NR
			return (pF*pF/(2.*me*T)) - (me/T);       
		} else if ((t<=1.e10) && (nfm>1.e-10)) { //degenerate UR
			return (pF/T) - (me/T);
		} else if (t>1.e10) {
			nq = 8.*pi*pow(T/(h*c),3.); //classical UR
			return -log(fabs(2.*nq/ne)) - (me/T);
		}
	} else if constexpr(species == 2) {
	//.........Find the guess
		nq = pow(2.*pi*mmu*T/(h*h*c*c),1.5); //cm^{-3}
		return -log(fabs(2.*nq/nLep));
	}

}


template<int species>
double rtsafe(const double nLep, const double T, const double guess) {
//Using a combination of Newton-Raphson and bisection, return the root of a function bracketed
//between x1 and x2. The root will be reﬁned until its accuracy is known within ˙xacc. funcd
//is a user-supplied struct that returns the function value as a functor and the ﬁrst derivative of
//the function at the point x as the function df (see text).
	const int MAXIT=150; //Maximum allowed number of iterations.
	const double xacc = 1.e-7; //set the accuracy for Newton Rapsondouble xh,xl;
	double x1, x2, mLep;

	if constexpr(species == 1) {
		mLep = me;
		x1 = e_tab_lim.eta_min;
		x2 = e_tab_lim.eta_max;
	} else if constexpr(species == 2) {
		mLep = mmu;
		x1 = mu_tab_lim.eta_min;
		x2 = mu_tab_lim.eta_max;
	}

	double xh, xl;
	double fl=n_net_f(x1, T, mLep)-nLep;
	double fh=n_net_f(x2, T, mLep)-nLep;

	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		printf("nLep = %.3e cm-3, T = %.3e MeV, fl = %.3e cm-3, fh = %.3e cm-3\n", nLep, T, fl, fh);
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
	double rts=guess; //0.5*(x1+x2);  //Initialize the guess for root,
	double dxold=fabs(x2-x1); //the “stepsize before last,”
	double dx=dxold;         //and the last step.
	double f=n_net_f(rts, T, mLep)-nLep;
	double df=n_net_df(rts, T, mLep);
	for (int j=0;j<MAXIT;j++) { //Loop over allowed iterations.
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
		if (fabs(dx) < xacc) return rts; //Convergence criterion.

		f=n_net_f(rts, T, mLep)-nLep;
		df=n_net_df(rts, T, mLep); //The one new function evaluation per iteration.

		if (f < 0.0) { //Maintain the bracket on the root.
			xl=rts;
		} else {
			xh=rts;
		}
	}
	throw("Maximum number of iterations exceeded in rtsafe");
}



template<int species>
double rtsafe_mod(const double nLep, const double T, const double guess, struct n_net_FD &FD) {
//Using a combination of Newton-Raphson and bisection, return the root of a function bracketed
//between x1 and x2. The root will be reﬁned until its accuracy is known within ˙xacc. funcd
//is a user-supplied struct that returns the function value as a functor and the ﬁrst derivative of
//the function at the point x as the function df (see text).
	const int MAXIT=150; //Maximum allowed number of iterations.
	const double xacc = 1.e-7; //set the accuracy for Newton Rapsondouble xh,xl;
	double x1, x2, mLep;
	
	if constexpr(species == 1) {
		mLep = me;
		x1 = e_tab_lim.eta_min;
		x2 = e_tab_lim.eta_max;
	} else if constexpr(species == 2) {
		mLep = mmu;
		x1 = mu_tab_lim.eta_min;
		x2 = mu_tab_lim.eta_max;
	}

	double xh, xl;
	double fl=n_net_f(x1, T, mLep)-nLep;
	double fh=n_net_f(x2, T, mLep)-nLep;

	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		printf("nLep = %.3e cm-3, T = %.3e MeV, fl = %.3e cm-3, fh = %.3e cm-3\n", nLep, T, fl, fh);
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
	double rts=guess; //0.5*(x1+x2);  //Initialize the guess for root,
	double dxold=fabs(x2-x1); //the “stepsize before last,”
	double dx=dxold;         //and the last step.
	n_net_fdf(rts, T, mLep, FD);
	double f  = FD.n_net - nLep; //n_net_f(rts, T, mLep)-nLep;
	double df = FD.dn_net; //n_net_df(rts, T, mLep);
	for (int j=0;j<MAXIT;j++) { //Loop over allowed iterations.
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) { //Bisect if Newton out of range, or not decreasing fast enough.
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			//printf("If 1\n");
			if (xl == rts) {
				//printf("%d", j);
				return rts;
			}
		} else { //Change in root is negligible. Newton step acceptable. Take it.
			dxold=dx;
			dx=f/df;
			double temp=rts;
			rts -= dx;
			//printf("If 2\n");
			if (temp == rts) {
				//printf("%d", j);
				return rts;
			}
		}
		if (fabs(dx) < xacc) { //Convergence criterion.
			//printf("%d", j);
			return rts;
		}
		//printf("dx = %.3e\n", dxold);
		//printf("eta = %.3e\n", rts);
		n_net_fdf(rts, T, mLep, FD);
		f  = FD.n_net - nLep; //n_net_f(rts, T, mLep)-nLep;
		df = FD.dn_net; //n_net_df(rts, T, mLep); //The one new function evaluation per iteration.


		//printf("rts = %.3e\n", rts);
		//printf("f = %.3e\n", f);
		//printf("df = %.3e\n\n", df);


		if (f < 0.0) { //Maintain the bracket on the root.
			xl=rts;
		} else {
			xh=rts;
		}
	}
	throw("Maximum number of iterations exceeded in rtsafe");
}



