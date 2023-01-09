#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <gsl/gsl_integration.h>
#include "../source/parameters.h"
//#include "../source/constants_bis.h"
#include "../source/tools/nr3.h"
#include "../source/eos/interp.h"
#include "../source/eos/eos_fermions.h"
#include "../source/eos/find_eta.h"

using namespace constants;

int main (){
        std::vector<double> n_arr, t_arr, eta_arr;
        double ne, T, theta;
	double eta, eta_NR, guess;
        double mu1, mu2, mu3;
        double p1, p2, p3;
        double a_p1, a_p2, a_p3;
        double e1, e2, e3;
        double a_e1, a_e2, a_e3;
        double s1, s2, s3;
        double a_s1, a_s2, a_s3;

	const int sp = 1; 

	/*Read EOS tables*/
        struct EOSeta eta_table = read_eta_table(sp);
        struct EOScomplete EOS_table = read_total_table(sp);

        //EOS arrays
        n_arr = eta_table.d;
        t_arr = eta_table.t;
        eta_arr = eta_table.eta;

	ne = 1.e0; //e-7;
	//eta = 3.84251563e+02;
	T = 6.e-2; //MeV
	theta = T/me;
	
	//printf("%.12e\n", compute_res(10., theta, 0.5));
	printf("Input quantities: ne = %.3e 1/fm^3, T = %.3e MeV (theta = %.10e)\n\n", ne, T, theta);

	//double h1 = compute_res(eta, theta, 0.5);
	//double h2 = compute_res(eta, theta, 1.5);
	//double h3 = compute_res(eta, theta, 2.5);
	//double h2 = compute_res_ed(eta, theta, k);
	//printf("f = %.15e, df = %.15e\n", h1, h2);
	//printf("h1 = %.6e, h2 = %.6e, h3 = %.6e\n", h1, h2, h3);
        
	//double test_1 = exp(500.);
	//long double test_2 = exp((long double) 500.);
	//long double test_3 = expl((long double) 500.);

	//long double tmp = -25.;
	//long double tt = 1. - exp(tmp)/(1.+exp(tmp));
	//long double tt1 = pow(1.+exp(tmp),-1.);
	//long double tt2 = 1. - exp(tmp);
	//printf("t1 = %.10e, t2 = %.10Le, t3 = %.10Le\n", test_1, test_2, test_3);
	//printf("tt = %.30Le, tt1 = %.30Le, tt2 = %.30Le\n", tt, tt1, tt2);

	//First method: completely on-the-fly (both eta and e.g. energy)
	printf("First method: completely on-the-fly\n");
	guess = find_guess_e(1.e39*ne, T);
	eta_NR = rtsafe(1.e39*ne, T, me, guess, eta1_e, eta2_e);
	printf("eta_NR = %.20e\n", eta_NR);
	double antieta = -(eta_NR+2./theta);
	printf("antieta_NR = %.20e\n", antieta);
	mu1 = T*eta_NR + me;
	//p1  = P_part(T, eta_NR, me);
	//e1  = e_part(T, eta_NR, me);
	//s1  = s_part(T, eta_NR, me);
	a_p1 = P_antipart(T, eta_NR, me);
	//a_e1 = e_antipart(T, eta_NR, me);
	//a_s1 = s_antipart(T, eta_NR, me);
	printf("Positron pressure = %.20e erg/cm^3\n", a_p1);
	//double F1 = compute_res(antieta, theta, 0.5);
	//double F2 = compute_res(antieta, theta, 1.5);
	//double F3 = compute_res(antieta, theta, 2.5);
	//printf("FD integral (k = 0.5, eta = %.3e, theta = %.3e): %.12e\n", antieta, theta, F1);
	//printf("FD integral (k = 1.5, eta = %.3e, theta = %.3e): %.12e\n", antieta, theta, F2);
	//printf("FD integral (k = 2.5, eta = %.3e, theta = %.3e): %.12e\n", antieta, theta, F3);

	exit(EXIT_SUCCESS);
	printf("Electrons: P = %.6e, e = %.6e, s = %.6e\n", p1, e1, s1);
	printf("Positrons: P = %.6e, e = %.6e, s = %.6e\n\n", a_p1, a_e1, a_s1);


	//Second method: eta interpolation + e.g. energy on-the-fly
	printf("Second method: eta interpolation, then on-the-fly\n");
	eta = eos_tintep(ne, T, n_arr, t_arr, eta_arr);
	printf("eta = %.6e\n", eta);
	mu2 = T*eta + me;
	p2  = P_part(T, eta, me);
	e2  = e_part(T, eta, me);
	s2  = s_part(T, eta, me);
	a_p2 = P_antipart(T, eta, me);
	a_e2 = e_antipart(T, eta, me);
	a_s2 = s_antipart(T, eta, me);
	printf("Electrons: P = %.6e, e = %.6e, s = %.6e\n", p2, e2, s2);
	printf("Positrons: P = %.6e, e = %.6e, s = %.6e\n\n", a_p2, a_e2, a_s2);

	//Third method: direct interpolation of e.g. energy
	printf("Third method: direct interpolation\n");
	mu3 = eos_tintep(ne, T, n_arr, t_arr, EOS_table.mu_part);
	p3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.P_part);
	e3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.e_part);
	s3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.s_part);
	a_p3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.P_antipart);
	a_e3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.e_antipart);
	a_s3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.s_antipart);
	printf("Electrons: P = %.6e, e = %.6e, s = %.6e\n", p3, e3, s3);
	printf("Positrons: P = %.6e, e = %.6e, s = %.6e\n\n", a_p3, a_e3, a_s3);

	//double K = 8*sqrt(2)*pi*pow(me/(h*c),3.);
	//printf("k = %.15e\n", K*MeV/(3.*pow(me,3.)));
	//printf("h = %.15e\n", K*MeV/(8.*sqrt(2.)*pi*pow(me,3.)));
	//printf("d = %.15e\n", K*MeV*kB/(pow(me,3.)));
	//printf("g = %.15e\n", K*mu/(pow(me,3.)));
	//std::exit(EXIT_SUCCESS);

	return 0;
}

