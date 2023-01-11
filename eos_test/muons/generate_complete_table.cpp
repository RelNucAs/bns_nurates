#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>

#include "../../source/parameters.h"
#include "../../source/eos/find_eta.h"
#include "../../source/eos/eos_fermions.h"
#include "../../source/eos/interp.h"

int main (){
	int n1, n2;
	if (parameters::HR == true){
		n1 = parameters::HR_tab.nnm;
		n2 = parameters::HR_tab.nt;
	} else {
		n1 = parameters::SR_tab.nnm;
		n2 = parameters::SR_tab.nt;
	}

	//int n1 = 750; //1500; //nnm
	//int n2 = 150; //300; //nt

	double nmu, T, eta;
	double guess;
	
	double mu_mu, p, n, e, s, an, ap, ae, as;

        double ymumin = 5.e-8;
        double ymumax = 5.e-1;
        
	double denmin = 1.e+08;
        double denmax = 2.e+16;

        double nmu_min = (ymumin*denmin)/(1.e39*mu);
        double nmu_max = (ymumax*denmax)/(1.e39*mu);

        double t_min = 5.e-02;
        double t_max = 1.05e+2;

        double lnm_edges[n1+1];
	double lt_edges[n2+1];
	double nm_array[n1];
	double t_array[n2];
	double eta_array[n1*n2];
	//std::vector<double> eta_arr(n1*n2);

        double mL = mmu;
	
       	string res;
	if (parameters::HR == true) {
		res = "_HR";
	} else {
		res = "";
	}
	string table_filename = abs_path + "eos_table/muons/eos_muons_complete_leo"+res+"_withn.txt";
        ofstream Iout(table_filename);
        Iout << n1 << "\n";
	Iout << n2 << "\n";
	Iout << std::scientific << setprecision(8);

	for (int i=0;i<n1+1;i++) lnm_edges[i] = log10(nmu_min) + ((double) i)*(log10(nmu_max)-log10(nmu_min))/((double) n1);
	for (int i=0;i<n1;i++) {
		nm_array[i]  = pow(10.,0.5*(lnm_edges[i+1]+lnm_edges[i]));
		Iout << nm_array[i] << " ";
	}
	Iout << "\n";

	for (int i=0;i<n2+1;i++) lt_edges[i] = log10(t_min)+ ((double) i)*(log10(t_max)-log10(t_min))/((double) n2);
	for (int i=0;i<n2;i++) {
		t_array[i]  = pow(10.,0.5*(lt_edges[i+1]+lt_edges[i]));
		Iout << t_array[i] << " ";
	}
	Iout << "\n";
	
        //struct EOSeta eta_table = read_eta_table();
        //eta_arr = eta_table.eta;

	for (int i=0;i<n1;i++) {
		for (int j=0;j<n2;j++) {
			printf("i = %d, j = %d\n", i, j);
			nmu = nm_array[i];
			T = t_array[j];

			guess = find_guess_mu(1.e39*nmu, T);
                        eta = rtsafe(1.e39*nmu, T, mL, guess, eta1_mu, eta2_mu);

                        eta_array[n2*i+j] = eta;
		}
	//	printf("\n");
	}

        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			mu_mu = t_array[j]*eta_array[n2*i+j] + mL;
			Iout << mu_mu << " ";
		}
		Iout << "\n";
	}

	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        n = n_part(t_array[j], eta_array[n2*i+j], mL);
                        Iout << n << " ";
                }
                Iout << "\n";
        }

	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        an = n_antipart(t_array[j], eta_array[n2*i+j], mL);
                        Iout << an << " ";
                }
                Iout << "\n";
        }

        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			p = P_part(t_array[j], eta_array[n2*i+j], mL);
			Iout << p << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			ap = P_antipart(t_array[j], eta_array[n2*i+j], mL);
			Iout << ap << " ";
		}
		Iout << "\n";
	}

	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			e = e_part(t_array[j], eta_array[n2*i+j], mL);
			Iout << e << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			ae = e_antipart(t_array[j], eta_array[n2*i+j], mL);
			Iout << ae << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			s = s_part(t_array[j], eta_array[n2*i+j], mL);
			Iout << s << " ";
		}
		Iout << "\n";
	}
	
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			as = s_antipart(t_array[j], eta_array[n2*i+j], mL);
			Iout << as << " ";
		}
		Iout << "\n";
	}

	Iout.close();

	return 0;
}

