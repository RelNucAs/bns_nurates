#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>

#include "source/parameters.h"
#include "source/eos/interp.h"
#include "source/eos/find_eta.h"
#include "source/eos/eos_fermions.h"

int main (){
	double r1, r2;
	double n_arr[nne], t_arr[nt], eta_arr[nne*nt]; //, e_arr[nne*nt];
	double ne, T;
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

	//Read EOS tables
	struct EOSeta eta_table = read_eta_table();
	struct EOScomplete EOS_table = read_total_table();

	//Input arrays
	for (int i=0;i<nne;i++) n_arr[i] = eta_table.d[i];
	for (int i=0;i<nt;i++)  t_arr[i] = eta_table.t[i];
	for (int i=0;i<nne*nt;i++) eta_arr[i] = eta_table.eta[i];
	//for (int i=0;i<nne*nt;i++) {
		//e_arr[i]   = EOS_table.e_part[i];
		//printf("%.3e\n", EOS_table.e_part[i]);
	//}
	
	

        std::ofstream fin("input_entries.txt");
        if (!fin)  cout << "File not found" << '\n';
	fin << std::scientific;
        
	std::ofstream f1("first_method.txt");
        if (!f1)  cout << "File not found" << '\n';
	f1 << std::scientific;

	std::ofstream f2("second_method.txt");
        if (!f2)  cout << "File not found" << '\n';
	f2 << std::scientific;
	
	std::ofstream f3("third_method.txt");
        if (!f3)  cout << "File not found" << '\n';
	f3 << std::scientific;
	
	//Initialization of random number generator
	srand(time(NULL));

	//Loop over (ne-T) grid
	for (int i=0;i<nne-1;i++) {
		for (int j=0;j<nt-1;j++) {
			printf("i = %d, j = %d\n", i, j);
			r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
	       		r2 = rand() / (double) RAND_MAX;
			//printf("r1 = %.3lf, r2 = %.3lf\n", r1, r2); //OK
			dn = n_arr[i+1] - n_arr[i];
			ne = n_arr[i] + r1*dn;		
			dt = t_arr[j+1] - t_arr[j];
			T  = t_arr[j] + r2*dt;
			
			fin << ne << "\t" << T << "\n"; 
			
			//First method: completely on-the-fly (both eta and e.g. energy)
			guess = find_guess_e(1.e39*ne, T);
			eta_NR = rtsafe(1.e39*ne, T, me, guess, eta1, eta2);
			mu1 = T*eta_NR + me;
			p1  = P_part(T, eta_NR, me);
			e1  = e_part(T, eta_NR, me);
			s1  = s_part(T, eta_NR, me);
			a_p1 = P_antipart(T, eta_NR, me);
			a_e1 = e_antipart(T, eta_NR, me);
			a_s1 = s_antipart(T, eta_NR, me);
	
			//Second method: eta interpolation + e.g. energy on-the-fly
			eta = eos_tintep(ne, T, n_arr, t_arr, eta_arr);
			mu2 = T*eta + me;
			p2  = P_part(T, eta, me);
			e2  = e_part(T, eta, me);
			s2  = s_part(T, eta, me);
			a_p2 = P_antipart(T, eta, me);
			a_e2 = e_antipart(T, eta, me);
			a_s2 = s_antipart(T, eta, me);

			//Third method: direct interpolation of e.g. energy
			mu2 = T*eta + me;
			p3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.P_part);
			e3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.e_part);
			s3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.s_part);
			a_p3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.P_antipart);
			a_e3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.e_antipart);
			a_s3  = eos_tintep(ne, T, n_arr, t_arr, EOS_table.s_antipart);
	 		
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

