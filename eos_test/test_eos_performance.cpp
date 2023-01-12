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
	double r1, r2;
	std::vector<double> n_arr, t_arr;
	std::vector<double> n_tab, t_tab, eta_tab;

	chrono_t begin, stop;
	chrono_int duration;
	
	double s1, s2, s3;
        double rho, ne, t;
	double ye = mu*1.e39;
        double guess;

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

	const int sp = 1; //1:electrons, 2:muons
	//Read EOS tables
        struct EOSeta eta_table = read_eta_table(sp);
        struct EOScomplete EOS_table = read_total_table(sp);

	std::array<double,9> lep_EOS;

        //Input table
        n_tab = eta_table.d;
        t_tab = eta_table.t;
        eta_tab = eta_table.eta;

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

		//guess = find_guess_e(1.e39*ne, t);
		//eta_NR = rtsafe(1.e39*ne, t, me, guess, eta1_e, eta2_e);
		//s1  = s_part(t, eta_NR, me);
		
		lep_EOS = eos_ferm_allmethods(ne, ye, t, sp, 1, &eta_table, &EOS_table);
	}
	stop = high_resolution_clock::now();
	for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, first method: " << duration.count()*1.e-3 << " s" << endl << endl;

	printf("Measuring performance: second method\n");
	begin = high_resolution_clock::now();
	//Second method: eta interpolation + e.g. energy on-the-fly
	for (int i=0; i<niter; i++) {
		ne = n_arr[i];
		t  = t_arr[i];
		
		//eta = eos_tintep(ne, t, n_tab, t_tab, eta_tab);
		//s2  = s_part(t, eta, me);
		
		lep_EOS = eos_ferm_allmethods(ne, ye, t, sp, 2, &eta_table, &EOS_table);
	}
	stop = high_resolution_clock::now();
	for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, second method: " << duration.count()*1.e-3 << " s" << endl << endl;
	
	printf("Measuring performance: third method\n");
	begin = high_resolution_clock::now();
	//Third method: direct interpolation of e.g. energy
	for (int i=0; i<niter; i++) {
		ne = n_arr[i];
		t  = t_arr[i];
		
		//s3  = eos_tintep(ne, t, n_tab, t_tab, EOS_table.s_part);
		
		lep_EOS = eos_ferm_allmethods(ne, ye, t, sp, 3, &eta_table, &EOS_table);
	}
	stop = high_resolution_clock::now();
	for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, third method: " << duration.count()*1.e-3 << " s" << endl << endl;

	return 0;
}

