#pragma once

#include <fstream>

#include "./nu_abs_em.h"
#include "./tools/nr3.h"
#include "./tools/gamma.h"
#include "./tools/gauss_wgts.h"

using namespace nuabsem;

void save_gaulag(int n, const Doub alf) {
	VecDoub_O x(n);
	VecDoub_O w(n);
	
	for (int i=0;i<n;i++) {
                x[i] = 0.;
                w[i] = 0.;
                //printf("%.3lf\n", x[i]);
                //printf("%.3lf\n", w[i]);
        }
	
        gaulag(x, w, alf);
	
        ofstream fout("gaulag_weights.txt");
        fout << "# Abscissas and weights for Gauss-Laguerre integration\n";
	
	for (int i=0;i<n;i++) {
		fout << std::scientific << x[i] << "\t"  << w[i] << "\n";
	}
        fout.close();
}


double em_integ(int n, double nb, double T, double lep_mass, double yp, double yn, double mu_l, double mu_hat, double deltaU, int g_type) {
       	//g_type = 0: Gauss-Laguerre 

	double integral = 0.;
	double absc, wgt;
	double spectrEm, spectrAb;
	string filename, dummyLine;

	if (g_type == 0) {
		filename = "gaulag_weights.txt";
	} else {
		printf("Index of Gaussian Quadrature method not recognized");
		exit (EXIT_FAILURE);
	}

        // open input file
        std::fstream fin(filename);

	if (!fin) {
		cout << "File regarding Gaussian quadrature not found" << '\n';
     	        exit (EXIT_FAILURE);
	}
        
	getline(fin, dummyLine); //skip header line

	for (int i=0;i++;i<n) {
                //std::stringstream ss(dummyLine);
                fin >> absc >> wgt;
		
		printf("%.3lf, %.3lf\n", absc, wgt);
                std::tie(spectrEm,spectrAb) = nu_n_abs(absc, nb, T, lep_mass, yp, yn, mu_l, mu_hat, deltaU);
                
		integral += wgt * spectrEm * pow(absc,3.);
		//printf("%.6e\n", integral);
        }

	return 2*pi*integral;
}

