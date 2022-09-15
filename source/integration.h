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

double compute_integral(int n, double *f, int g_type) {
	double x[n];
	double w[n];

	read_wgts(n, x, w, g_type);
	
	return do_integration(n, w, f);

}

double compute_split_gauleg(int n, double *f1, double *f2, double a) {
	double x[n], w[n];
	double I1, I2;
	double f3[n];

	read_wgts(n, x, w, 0);

	for (int i=0;i<n;i++) f3[i] = f2[i] / (x[i] * x[i]);

	I1 = a * do_integration(n, w, f1);
	I2 = a * do_integration(n, w, f3);	
	
	return I1 + I2;
}


double split_value(double T, double eta_e) {
	return T * Fermi_integral_p5(eta_e) / Fermi_integral_p4(eta_e);
}

