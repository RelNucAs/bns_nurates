// \file NN_brem.cpp
// \brief Compute the nucleon-nucleon bremsstrahlung using the analytic fitting formula
//        in Hannestad & Raffelt 1998, Apj, 507, 339 (https://iopscience.iop.org/article/10.1086/306303/pdf)

#include <iostream>

#include "constants.hpp"       // Header file for physical constants
#include "physics/NN_brem.hpp" // Header file for NN bremsstrahlung

using namespace std;
using namespace constants;
using namespace brems;

// Input arguments:
//   w, wp   : neutrino (antineutrino) energy [MeV]
//   rho     : matter density [g/cm3]
//   tmev    : temperature [MeV]
//   xn      : free neutron mass fraction
//   xp      : free proton mass fraction
//
//
// Output arguments:
//   s1, s2 : differential neutrino pair production kernel (prod/abs kernel [MeV-1])
 
/* todo: rewrite code to use nb instead of rho + understand origin of rho_14*(xn+xp) in s1/s2 */


/* Compute NN bremsstrahlung scattering kernel */
std::tuple<double,double> brem_kernel(const double w, const double wp, const double rho, const double tmev, const double xn, const double xp) {
	/* Initialize */
	const double rho_14 = rho/1.0e+14; //density in units of 1.0e+14 g/cm3
	const double t_10   = tmev/10.;    //temperature in units of 10 MeV

	/* Compute eta_star, the neutron effective degeneracy parameter */
	const double eta_star = 3.04 * pow(rho_14,2./3.) / t_10; //Eq.(36)

	/* Compute gamma, the spin fluctuation rate */
	const double gamma = 8.6 * rho_14 / sqrt(t_10); //Eq.(37)

	/* Compute y, the pion mass parameter */
	double y = 1.94 / t_10; //Eq.(38)

	/* Compute x, the dimensionless neutrino energy sum */
	double x = (w+wp) / tmev; 

	//cout << "Before: ";
	//cout << "x = " << x << ", y = " << y << endl;
		
	/* Compute the dimensionless fitting parameter, sb */
      	const double sb = s_brem(x,y,eta_star);  //N.B.: x and y values are not changed by the function

	/* Compute the dimensionless fitting parameter, gb */
	const double gb = g_brem(y,eta_star);    //N.B.: x and y values are not changed by the function
	
	//cout << "After: ";
	//cout << "x = " << x << ", y = " << y << endl;
	
	/* Compute the differential absorption kernel, s_a */
	double s2 = gamma / (x*x + pow(0.5*gamma*gb,2.)) * sb/tmev * rho_14 * (xn+xp); //differential? production? kernel (Eq.(35))
	// !Understand rho_14 * (xn+xp) origin!
	double s1 = s2 * fexp(-x); //differential? absorption? kernel (from detailed balance)
	
	/* Multiply absorption kernel by antineutrino distribution function before integrating over wp */
	//s2 = s2 * fermi(wp/tmev);

	return std::make_tuple(s1,s2);
}


// Compute s-component of the kernel for neutrino bremsstrahlung
// and inelastic neutrino scattering in a nucleon field
//
// Input arguments:
//  y           : m_{pi}^{2}/M_{N}c^{2}kT (pion mass parameter)
//  eta_star    : neutron degeneracy parameter
//
// Input arguments:
// g           : dimensionless quantity related to S_{sigma}

double s_brem(double x, double y, double eta_star) {

	/* Prevent singular behavior of terms */
	if (x >= 0.) {
		x = max(x,x_min);
	} else {
		x = min(x,-x_min);
	}
	y        = max(y,y_min);
	eta_star = max(eta_star,eta_min);

	/* Compute S_ND, nondegenerate approximation */
	const double s_ND_num   = 2. * pi1_2 * pow(x + 2. - exp(-y/12.),1.5) *
	   (x*x + 2.*x*y + 5./3.*y*y + 1.);
	const double s_ND_denom = pi1_2 + pow(pi1_8 + x + y,4.);
	double s_ND = s_ND_num/s_ND_denom; //Eq.(45)

	if (x < 0.) s_ND = s_ND * fexp(-x); //Eq.(41)

	/* Compute S_D, degenerate approximation */
	const double u     = sqrt(y / (2.*eta_star)) + 1.e-10;
	const double u2    = u*u;
	const double u_arg = u2 / (2.*sqrt(2.*u2 + 4.));
	const double f_u   = 1. - 5./6.*u* atan(2./u) + u2/(3.*(u2+4.)) +
	   atan(1./u_arg)*u_arg/3.; //Eq.(47)
	const double s_D   = 3. * pow(0.5*pi,2.5) * pow(eta_star,-2.5) *
	   (x*x + 4.*pi2) * x * f_u / (4. * pi2 * (1. - fexp(-x))); //Eq.(46)

	//if (s_D < 0.) {
	//	cout <<                "s_D = " <<         s_D << endl;
	//	cout << "one minus fexp(-x) = " << 1.-fexp(-x) << endl;
	//	cout <<                  "x = " <<           x << endl;
	//	cout <<                  "y = " <<           y << endl;
	//	cout <<                  "u = " <<           u << endl;
        //	cout <<           "eta_star = " <<    eta_star << endl;
	//	cout <<                "f_u = " <<         f_u << endl;
	//}

//      write(6,*)S_D,eta_star,f_u,x,fexp(-x)
//      write(6,*)''

	/* Compute F */
	const double F_denom = (3. +  pow(x-1.2,2.) + pow(x,-4.)) *
	   (1. + eta_star*eta_star) * (1. + pow(y,4.));
	const double F       = 1. + 1./F_denom; //Eq.(50)

	/* Compute G */
	const double G = 1. - 0.0044*pow(x,1.1) * y / (0.8 + 0.06*pow(y,1.05)) *
	   sqrt(eta_star) / (eta_star + 0.2); //Eq.(50)

	/* Compute C */
	const double H = 0.1 * eta_star / (2.39 + 0.1*pow(eta_star,1.1)); //Eq.(50)
	const double C = 1.1 * pow(x,1.1) * H /
	   (2.3 + H*pow(x,0.93) + 0.0001*pow(x,1.2)) * 30. / (30. + 0.005*pow(x,2.8)); //Eq.(50)

	/* Compute p */
	const double p = 0.67 + 0.18*pow(y,0.4); //Eq.(50)
//      write(6,*)s_ND,S_D,p,C,G,F

	if (s_ND < 0.) cout << "s_ND = " << s_ND << endl;
	if (s_D < 0.)  cout << "s_D  = " << s_D  << endl;

	/* Compute s */
	return pow( pow(s_ND,-p) + pow(s_D,-p) , -1./p) * F * (1. + C*G); //Eq.(49)
}



// Compute g-component of the kernel for neutrino bremsstrahlung
// and inelastic neutrino scattering in a nucleon field
//
// Input arguments:
//   y           : m_{pi}^{2}/M_{N}c^{2}kT (pion mass parameter)
//   eta_star    : neutron degeneracy parameter
//
// Input arguments:
//   g           : dimensionless quantity related to S_{sigma}

double g_brem(double y, double eta_star) {

	/* Prevent singular behavior of terms */
	y        = max(y,y_min);
        eta_star = max(eta_star,eta_min);

	/* Compute alpha1 */
	const double y2 = y*y;
	const double eta_star_1 = 1./eta_star;
	const double alpha1_denom = 25.*y2 + 1.;
	
	const double alpha1 = (0.5 + eta_star_1) / (1. + eta_star_1) *
	   1./alpha1_denom + (0.5 + eta_star/15.6) * 25.*y2 / alpha1_denom; //Eq.(53)

	/* Compute alpha2 */
	const double alpha2 = (0.63 + 0.04*pow(eta_star,1.45)) / (1. + 0.02*pow(eta_star,2.5)); //Eq.(53)

	/* Compute alpha3 */
	const double alpha3 = 1.2 * fexp(0.6*eta_star - 0.4*pow(eta_star,1.5)); //Eq.(53)

	/* Compute p1 */
	const double p1 = (1.8 + 0.45*eta_star) / (1. + 0.15*pow(eta_star,1.5)); //Eq.(53)

	/* Compute p2 */
	const double p2 = 2.3 - 0.05*eta_star / (1. + 0.025*eta_star); //Eq.(53)

	/* Compute g */
	return (alpha1 + alpha2*pow(y,p1)) / (1. + alpha3*pow(y,p2) + alpha2*pow(y,p1+2.)/13.75); //Eq.(52)

}

/* Safe exp function to avoid overflow */
double fexp(const double x) {
	const double expmax = 0.1*std::numeric_limits<double>::max_exponent; //DBL_MAX_EXP
	return exp(min(expmax,max(-expmax,x)));
}
