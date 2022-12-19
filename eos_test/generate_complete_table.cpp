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
#include "../source/eos/interp.h"

int main (){
	int n1 = 700; //1400; //nne
	int n2 = 150; //300; //nt

	double ne, T, eta;
	double guess;
	
	double mue, p, e, s, ap, ae, as;

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
	double eta_array[n1*n2];
	std::vector<double> eta_arr(n1*n2);

        string table_filename = "../eos_table/eos_electrons_complete_with_eta_ele.txt"; //eos_electrons_complete_leo_HR.txt";
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

	for (int i=0;i<n2+1;i++) lt_edges[i] = log10(t_min)+ ((double) i)*(log10(t_max)-log10(t_min))/((double) n2);
	for (int i=0;i<n2;i++) {
		t_array[i]  = pow(10.,0.5*(lt_edges[i+1]+lt_edges[i]));
		Iout << t_array[i] << " ";
	}
	Iout << "\n";
	
        struct EOSeta eta_table = read_eta_table();
        eta_arr = eta_table.eta;

	for (int i=0;i<n1;i++) {
		for (int j=0;j<n2;j++) {
			printf("i = %d, j = %d\n", i, j);
	//		ne = ne_array[i];
	//		T = t_array[j];

	//		guess = find_guess_e(1.e39*ne, T);
        //                eta = rtsafe(1.e39*ne, T, me, guess, eta1, eta2);
                        
	         	eta_array[n2*i+j] = eta_arr[n2*i+j];
                        //s = s_part(t_array[j], eta_array[n2*i+j], me);
		}
	//	printf("\n");
	}

        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			mue = t_array[j]*eta_array[n2*i+j] + me;
			Iout << mue << " ";
		}
		Iout << "\n";
	}

        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			p = P_part(t_array[j], eta_array[n2*i+j], me);
			Iout << p << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			ap = P_antipart(t_array[j], eta_array[n2*i+j], me);
			Iout << ap << " ";
		}
		Iout << "\n";
	}

	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			e = e_part(t_array[j], eta_array[n2*i+j], me);
			Iout << e << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			ae = e_antipart(t_array[j], eta_array[n2*i+j], me);
			Iout << ae << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			s = s_part(t_array[j], eta_array[n2*i+j], me);
			Iout << s << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			as = s_antipart(t_array[j], eta_array[n2*i+j], me);
			Iout << as << " ";
		}
		Iout << "\n";
	}

	Iout.close();

	return 0;
}

