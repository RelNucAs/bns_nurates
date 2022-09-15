#pragma once

#include <fstream>

#include "./tools/fermi_integrals.h"
#include "./tools/nr3.h"
#include "./tools/gamma.h"
#include "./tools/gauss_wgts.h"


void save_gauleg(int n, const Doub x1, const Doub x2) {
        VecDoub_O x(n);
        VecDoub_O w(n);

	char fileHeader[100];

        for (int i=0;i<n;i++) {
                x[i] = 0.;
                w[i] = 0.;
        }

	gauleg(x1, x2, x, w);

        ofstream fout("./input/gauleg_weights.txt");
        sprintf(fileHeader, "# Abscissas and weights for Gauss-Legendre integration from x1 = %.3lf to x2 = %.3lf\n", x1, x2);
	
	fout << fileHeader;

        for (int i=0;i<n;i++) {
                fout << std::scientific << x[i] << "\t"  << w[i] << "\n";
        }
        fout.close();
}


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
	
        ofstream fout("./input/gaulag_weights.txt");
        fout << "# Abscissas and weights for Gauss-Laguerre integration\n";
	
	for (int i=0;i<n;i++) {
		fout << std::scientific << x[i] << "\t"  << w[i] << "\n";
	}
        fout.close();
}

void read_wgts(int n, double *x, double *w, int g_type) {
	//g_type = 0: Gauss-Legendre
	//g_type = 1: Gauss-Laguerre
	
	string filename, dummyLine;

	if (g_type == 0) {
                filename = "./input/gauleg_weights.txt";
        } else if (g_type == 1) {
                filename = "./input/gaulag_weights.txt";
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

        for (int i=0;i<n;i++) {
                std::stringstream ss(dummyLine);
                fin >> x[i] >> w[i];
	}

	fin.close();
}


double do_integration(int n, double *wgt, double *spectrArray) {
	
	double integral = 0.;

	for (int i=0;i<n;i++) {
		//printf("%.3lf, %.3lf\n",wgt[i], spectrArray[i]);
		integral += wgt[i] * spectrArray[i];
        }

	return integral;

}



double compute_integral(int n, double *spectrArray, int g_type, double tmp) {
	double x[n];
	double w[n];

	printf("%.3e, %.3e\n", tmp, Fermi_integral_p5(tmp));
	read_wgts(n, x, w, g_type);
	
	//for (int i=0;i<n;i++) {
	//	printf("%.3lf, %.3lf\n",x[i], w[i]);
	//}
	return do_integration(n, w, spectrArray);

}


double split_value(


