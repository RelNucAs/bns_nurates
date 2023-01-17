#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <iomanip>

#include "../source/parameters.h"
#include "../source/eos/find_eta.h"

using namespace std;
using namespace parameters;

template<int species>
void make_eta_table() {
	int n1, n2;
	double   dmin,   dmax;
	double   ymin,   ymax;
	double   tmin,   tmax;
	double etamin, etamax;
	string sp, res;
	if constexpr(species == 1) {
		dmin = e_tab_lim.den_min;
		dmax = e_tab_lim.den_max;
		ymin = e_tab_lim.y_min;
		ymax = e_tab_lim.y_max;
		tmin = e_tab_lim.t_min;
		tmax = e_tab_lim.t_max;
		etamin = e_tab_lim.eta_min;
		etamax = e_tab_lim.eta_max;
		sp = "electrons";
		if (parameters::HR == true){
	                n1 = HR_tab.nne;
        	        n2 = HR_tab.nt;
			res = "_HR";
        	} else {
                	n1 = SR_tab.nne;
                	n2 = SR_tab.nt;
			res = "";
        	}

	} else if constexpr(species == 2) {
		dmin = mu_tab_lim.den_min;
		dmax = mu_tab_lim.den_max;
		ymin = mu_tab_lim.y_min;
		ymax = mu_tab_lim.y_max;
		tmin = mu_tab_lim.t_min;
		tmax = mu_tab_lim.t_max;
		etamin = mu_tab_lim.eta_min;
		etamax = mu_tab_lim.eta_max;
		sp = "muons";
                if (parameters::HR == true){
                        n1 = HR_tab.nnm;
                        n2 = HR_tab.nt;
			res = "_HR";
                } else {
                        n1 = SR_tab.nnm;
                        n2 = SR_tab.nt;
			res = "";
                }

	} else {
		cout << "Make_eta_table: Error in species identifier" << endl;
		exit(EXIT_FAILURE);
	}

	const double nmin = (ymin*dmin)/(1.e39*mu);
        const double nmax = (ymax*dmax)/(1.e39*mu);

	const double log_nmin = log10(nmin);
	const double log_nmax = log10(nmax);

	const double log_tmin = log10(tmin);
	const double log_tmax = log10(tmax);

	double nLep, temp, eta;
	double guess;

	//std::array<double,n1+1>  ln_edges;
	//std::array<double,n2+1>  lt_edges;
	//std::array<double,n1>     n_array;
	//std::array<double,n2>     t_array;

	std::vector<double> ln_edges;
	std::vector<double> lt_edges;
	std::vector<double>  n_array;
	std::vector<double>  t_array;

        string table_filename = abs_path + "eos_table/" + sp + "/eos_" + sp + "_leo"+res+".txt";
 
 	cout << "Generating eta table for " << sp << " with :" << endl;
	cout << "# n points = " << n1 << ", # t points = " << n2 << " (HR = " << HR << ")" << endl;
	cout << "n_min = "   <<   nmin << " fm-3, n_max = " <<   nmax << " fm-3" << endl;
	cout << "t_min = "   <<   tmin <<  " MeV, t_max = " <<   tmax <<  " MeV" << endl;
	cout << "eta_min = " << etamin <<    ", eta_max = " << etamax <<            endl;

	ofstream Iout(table_filename);
        Iout << n1 << "\n";
	Iout << n2 << "\n";
	Iout << std::scientific << setprecision(8);

	for (int i=0;i<n1+1;i++) ln_edges.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n1));
	for (int i=0;i<n2+1;i++) lt_edges.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n2));

	for (int i=0;i<n1;i++) {
		n_array.push_back(0.5*(ln_edges[i+1]+ln_edges[i]));
		Iout << n_array[i] << " ";
	}
	Iout << "\n";


	for (int i=0;i<n2;i++) {
		t_array.push_back(0.5*(lt_edges[i+1]+lt_edges[i]));
		Iout << t_array[i] << " ";
	}
	Iout << "\n";
	
	cout << "Loop over number density array (from 0 to " << (n1-1) << ")" << endl;
	for (int i=0;i<n1;i++) {
		cout << "i = " << i << endl;
		for (int j=0;j<n2;j++) {
			nLep = pow(10., n_array[i]);
			temp = pow(10., t_array[j]);

			guess = find_guess_eta<species>(1.e39*nLep, temp);
                        eta = rtsafe<species>(1.e39*nLep, temp, guess);
         
	       		Iout << eta << " ";
		}
		Iout << "\n";
	}		
	
	Iout.close();
	cout << "Done!" << endl << endl;
	return;
}



int main (){

	make_eta_table<1>();
	make_eta_table<2>();

	return 0;

}

