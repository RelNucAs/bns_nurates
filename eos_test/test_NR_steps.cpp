#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
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

	//std::ofstream fout(sp+"/test_NR_steps_"+sp+res+".txt");
	//if (!fout)  cout << "File not found" << '\n';
	//fout << std::scientific;
        
	//Initialization of random number genLrator
	srand(time(NULL));

	//Loop over (nL-T) grid
	for (int i=0;i<nn-1;i++) {
		for (int j=0;j<nt-1;j++) {
			//printf("i = %d, j = %d\n", i, j);
			r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
	       		r2 = rand() / (double) RAND_MAX;
			//printf("r1 = %.3lf, r2 = %.3lf\n", r1, r2); //OK
			dn = n_arr[i+1] - n_arr[i];
			nL = n_arr[i] + r1*dn;		
			dt = t_arr[j+1] - t_arr[j];
			T  = t_arr[j] + r2*dt;
			
			//fout << nL << "\t" << T << "\n"; 
			
			//First method: completely on-the-fly (both eta and e.g. energy)
			if (Nsp == 1) {
				guess = find_guess_e(1.e39*nL, T);
				eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_e, eta2_e);
			} else if (Nsp == 2) {
				guess = find_guess_mu(1.e39*nL, T);
				eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_mu, eta2_mu);
			}
			printf("\t");
			guess = eos_tintep(nL, T, n_arr, t_arr, eta_arr);
			
			if (Nsp == 1) {
                                eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_e, eta2_e);
                        } else if (Nsp == 2) {
                                eta_NR = rtsafe(1.e39*nL, T, mL, guess, eta1_mu, eta2_mu);
                        }
			printf("\n");
			
			//fout  a_e3 << "\t" << a_s3 << "\n";
		}
	}

	//fout.close();

	return 0;
}

