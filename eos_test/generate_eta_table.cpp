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
#include "../source/eos/find_eta.h"
#include "../source/eos/eos_fermions.h"

int main (){
	int n1 = 1400; //700; //nne
	int n2 = 300; //150; //nt

	double ne, T, eta;
	double guess;

	double yemin = 5.e-3;
	double yemax = 5.5e-1;
	
	double denmin = 5.e+3;
	double denmax = 2.e+16;

	double ne_min = yemin*denmin/(1.e39*mu);
	double ne_max = yemax*denmax/(1.e39*mu);
	
	double t_min = 5.e-2;
	double t_max = 1.05e+2;
	
        double lne_edges[n1+1];
	double lt_edges[n2+1];
	double ne_array[n1];
	double t_array[n2];

        string table_filename = "../eos_table/eos_electrons_leo_HR.txt";
        ofstream Iout(table_filename);
        Iout << n1 << "\n";
	Iout << n2 << "\n";
	Iout << std::scientific << setprecision(8);

	for (int i=0;i<n1+1;i++) lne_edges[i] = log10(ne_min) + ((double) i)*(log10(ne_max)-log10(ne_min))/((double) n1);
	for (int i=0;i<n1;i++) {
		ne_array[i]  = pow(10.,0.5*(lne_edges[i+1]+lne_edges[i]));
		Iout << ne_array[i] << " ";
	}
	Iout << "\n";

	for (int i=0;i<n2+1;i++) lt_edges[i] = log10(t_min) + ((double) i)*(log10(t_max)-log10(t_min))/((double) n2);
	for (int i=0;i<n2;i++) {
		t_array[i]  = pow(10.,0.5*(lt_edges[i+1]+lt_edges[i]));
		Iout << t_array[i] << " ";
	}
	Iout << "\n";
	
	for (int i=0;i<n1;i++) {
		for (int j=0;j<n2;j++) {
			printf("i = %d, j = %d\n", i, j);
			ne = ne_array[i];
			T = t_array[j];

			guess = find_guess_e(1.e39*ne, T);
                        eta = rtsafe(1.e39*ne, T, me, guess, eta1, eta2);
         
	       		Iout << eta << " ";
		}
		Iout << "\n";
	}		
	
	Iout.close();

	return 0;
}

