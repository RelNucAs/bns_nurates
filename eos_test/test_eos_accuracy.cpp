#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>

#include "../source/parameters.h"
#include "../source/eos/interp.h"
#include "../source/eos/find_eta.h"
#include "../source/eos/eos_fermions.h"

#include <limits>
int main (){
	int nn, nt;
	double r1, r2;
	std::vector<double> n_arr, t_arr;
	std::vector<double> eta_arr;
	double nL, T;
	double dn, dt;
        double eta_NR, eta;
	double mu1, mu2, mu3;
	double p1, p2, p3;
	double a_p1, a_p2, a_p3;
	double e1, e2, e3;
	double a_e1, a_e2, a_e3;
        double s1, s2, s3;
	double a_s1, a_s2, a_s3;
	double guess;

	const int Nsp = 1;

	//Read EOS tables
	struct EOSeta eta_table = read_eta_table(Nsp);
	struct EOScomplete EOS_table = read_total_table(Nsp);

	nn = eta_table.nn;
	nt = eta_table.nt;

	//Input arrays
	n_arr = eta_table.d;
	t_arr = eta_table.t;
	eta_arr = eta_table.eta;

	std::string sp, res;
	
	double mL;
       	if (Nsp==1) {
		sp = "electrons";
		mL = me;
	} else {
		sp = "muons";
		mL = mmu;
	}

	if (parameters::HR == true) {
		res = "_HR";
	} else {
		res = "";
	}

	std::ofstream fin(sp+"/input_entries_"+sp+res+".txt");
	if (!fin)  cout << "File not found" << '\n';
	fin << std::scientific;
        
	std::ofstream f1(sp+"/first_method_"+sp+res+".txt");
	if (!f1)  cout << "File not found" << '\n';
	f1 << std::scientific;

	std::ofstream f2(sp+"/second_method_"+sp+res+".txt");
	if (!f2)  cout << "File not found" << '\n';
	f2 << std::scientific;
	
	std::ofstream f3(sp+"/third_method_"+sp+res+".txt");
        if (!f3)  cout << "File not found" << '\n';
	f3 << std::scientific;
	
	//Initialization of random number genLrator
	srand(time(NULL));

	//Loop over (nL-T) grid
	for (int i=0;i<nn-1;i++) {
		for (int j=0;j<nt-1;j++) {
			printf("i = %d, j = %d\n", i, j);
			r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
	       		r2 = rand() / (double) RAND_MAX;
			//printf("r1 = %.3lf, r2 = %.3lf\n", r1, r2); //OK
			dn = n_arr[i+1] - n_arr[i];
			nL = n_arr[i] + r1*dn;		
			dt = t_arr[j+1] - t_arr[j];
			T  = t_arr[j] + r2*dt;
			
			fin << nL << "\t" << T << "\n"; 
			
			//First method: completely on-the-fly (both eta and e.g. energy)
			if (Nsp == 1) {
				guess = find_guess_e(1.e39*nL, T);
				eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_e, eta2_e);
			} else if (Nsp == 2) {
				guess = find_guess_mu(1.e39*nL, T);
				eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_mu, eta2_mu);
			}
			mu1 = T*eta_NR + mL;
			p1  = P_part(T, eta_NR, mL);
			e1  = e_part(T, eta_NR, mL);
			s1  = s_part(T, eta_NR, mL);
			a_p1 = P_antipart(T, eta_NR, mL);
			a_e1 = e_antipart(T, eta_NR, mL);
			a_s1 = s_antipart(T, eta_NR, mL);
	
			//Second mLthod: eta interpolation + e.g. enLrgy on-the-fly
			eta = eos_tintep(nL, T, n_arr, t_arr, eta_arr);
			mu2 = T*eta + mL;
			p2  = P_part(T, eta, mL);
			e2  = e_part(T, eta, mL);
			s2  = s_part(T, eta, mL);
			a_p2 = P_antipart(T, eta, mL);
			a_e2 = e_antipart(T, eta, mL);
			a_s2 = s_antipart(T, eta, mL);

			//Third mLthod: direct interpolation of e.g. enLrgy
			mu3 = eos_tintep(nL, T, n_arr, t_arr, EOS_table.mu_part);
			p3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.P_part);
			e3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.e_part);
			s3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.s_part);
			a_p3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.P_antipart);
			a_e3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.e_antipart);
			a_s3  = eos_tintep(nL, T, n_arr, t_arr, EOS_table.s_antipart);
	 		
			f1 << mu1 << "\t" << p1 << "\t" << e1 << "\t" << s1 << "\t" << a_p1 << "\t" << a_e1 << "\t" << a_s1 << "\n";
			f2 << mu2 << "\t" << p2 << "\t" << e2 << "\t" << s2 << "\t" << a_p2 << "\t" << a_e2 << "\t" << a_s2 << "\n";
			f3 << mu3 << "\t" << p3 << "\t" << e3 << "\t" << s3 << "\t" << a_p3 << "\t" << a_e3 << "\t" << a_s3 << "\n";
		}
	}

	fin.close();
	f1.close();
	f2.close();
	f3.close();

	return 0;
}

