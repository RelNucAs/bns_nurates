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

template<int species>
void test_accuracy() {
	double r1, r2;
	double nLep, T;
	double dn, dt;

	//Read EOS tables
	struct EOSeta eta_table = read_eta_table(species);
	struct EOScomplete EOS_table = read_total_table(species);

	const int nn = eta_table.nn;
	const int nt = eta_table.nt;

	std::vector<double> n_arr = eta_table.nL;
	std::vector<double> t_arr = eta_table.t;
	std::array<double,9> eos_array;

	std::string sp, res;
	
       	if constexpr(species==1) {
		sp = "electrons";
	} else if constexpr(species==2){
		sp = "muons";
	}

	if (parameters::HR == true) {
		res = "_HR";
	} else {
		res = "";
	}

	std::ofstream fout(sp+"/test_accuracy_"+sp+res+".txt");
	if (!fout)  cout << "File not found" << endl;
	fout << std::scientific;
        
        cout << "Testing EOS accuracy for " << sp << "...." << endl;
	
	//Initialization of random number genLrator
	srand(time(NULL));

	//Loop over (nLep-T) grid
        cout << "Loop over number density array (from 0 to " << (nn-2) << ")" << endl;
	for (int i=0;i<nn-1;i++) {
		cout << "i = " << i << endl;
		for (int j=0;j<nt-1;j++) {
			//printf("i = %d, j = %d\n", i, j);
			r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
	       		r2 = rand() / (double) RAND_MAX;
			dn = n_arr[i+1] - n_arr[i];
			dt = t_arr[j+1] - t_arr[j];
			nLep = pow(10.,n_arr[i] + r1*dn);		
			T    = pow(10.,t_arr[j] + r2*dt);
			
			fout << nLep << " " << T; 
			
			//First method: completely on-the-fly
			eos_array = eos_ferm_array<1,species>(nLep, T, eta_table, EOS_table);
			for (int i=0;i<9;i++) fout << eos_array[i] << " ";
	
			//Second method: eta interpolation + on-the-fly computation
			eos_array = eos_ferm_array<2,species>(nLep, T, eta_table, EOS_table);
			for (int i=0;i<9;i++) fout << eos_array[i] << " ";

			//Third method: direct interpolation
			eos_array = eos_ferm_array<3,species>(nLep, T, eta_table, EOS_table);
			for (int i=0;i<9;i++) fout << eos_array[i] << " ";
			fout << endl;
		}
	}

	cout << endl;
	fout.close();

	return;
}

int main() {
	test_accuracy<1>(); //1: electrons
	test_accuracy<2>(); //2: muons
	return 0;
}
