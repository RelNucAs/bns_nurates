#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <chrono>
#include <bits/stdc++.h>
#include <limits>

#include "../source/parameters.h"
#include "../source/eos/interp.h"
#include "../source/eos/find_eta.h"
#include "../source/eos/eos_fermions.h"

using namespace std::chrono;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> chrono_t; //Chrono type.
typedef std::chrono::duration<long int, std::ratio<1, 1000> > chrono_int;

int main (){
	int niter = 10000;

	const int sp = 1; //1:electrons, 2:muons

	chrono_t begin, stop;
	chrono_int duration;
	
	double r1, r2;
        double ne, t;

        double yemin = 5.e-3;
        double yemax = 5.5e-1;

        double denmin = 5.e+3;
        double denmax = 2.e+16;

        double ne_min = yemin*denmin/(1.e39*mu);
        double ne_max = yemax*denmax/(1.e39*mu);

        double t_min = 5.e-2;
        double t_max = 1.05e+2;

	double dn = ne_max-ne_min;
	double dt = t_max-t_min;

	std::vector<double> n_arr, t_arr;
	std::array<double,9> lep_EOS;

	//Read EOS tables
        struct EOSeta eta_table = read_eta_table(sp);
        struct EOScomplete EOS_table = read_total_table(sp);
	


	//Initialization of random number generator
	printf("Measuring numerical performance of EOS module, number of iterations: %d\n\n", niter);

	//Initialization of random number generator
	srand(time(NULL));

	//Loop over (ne-t) grid
	for (int i=0; i<niter; i++) {
		r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
		r2 = rand() / (double) RAND_MAX;

		n_arr.push_back(ne_min + r1*dn);
		t_arr.push_back(t_min  + r2*dt);
	}






	printf("Measuring performance: first method\n");
	begin = high_resolution_clock::now();
	//First method: completely on-the-fly (both eta and e.g. energy)
	for (int i=0; i<niter; i++) {
		ne = n_arr[i];
		t  = t_arr[i];
		lep_EOS = eos_ferm_array<1,sp>(ne, t, eta_table, EOS_table);

	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, first method: " << duration.count()*1.e-3 << " s" << endl << endl;




	printf("Measuring performance: second method\n");
	begin = high_resolution_clock::now();
	//Second method: eta interpolation + e.g. energy on-the-fly
	for (int i=0; i<niter; i++) {
		ne = n_arr[i];
		t  = t_arr[i];
		lep_EOS = eos_ferm_array<2,sp>(ne, t, eta_table, EOS_table);
	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, second method: " << duration.count()*1.e-3 << " s" << endl << endl;




	printf("Measuring performance: third method\n");
	begin = high_resolution_clock::now();
	//Third method: direct interpolation of e.g. energy
	for (int i=0; i<niter; i++) {
		ne = n_arr[i];
		t  = t_arr[i];
		lep_EOS = eos_ferm_array<3,sp>(ne, t, eta_table, EOS_table);
	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, third method: " << duration.count()*1.e-3 << " s" << endl << endl;

	return 0;
}

