//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  integration_quadrature_io.c
//  \brief write and read quadratures for integration

#include <stdio.h>
#include "../bns_nurates.h"
#include "integration.h"
#include "../functions/functions.h"

void save_gauss_legendre(int n, const double x1, const double x2) {

  char fileHeader[100];

  MyQuadrature quad;
  quad.type = kGauleg;
  quad.n = n;
  quad.dim = 1;

  gauss_legendre(&quad);

  char outname[50];
  sprintf(outname, "../quadratures/gauleg_weights_n_%d.txt", n);

  ofstream fout(abs_path+outname);
  sprintf(fileHeader, "# Abscissas and weights for Gauss-Legendre integration from x1 = %.3lf to x2 = %.3lf\n", x1, x2);

  fout << fileHeader;

  fout << std::scientific;
  fout << std::setprecision(16);
        for (int i=0;i<n;i++) {
                fout << x[i] << "\t"  << w[i] << "\n";
        }
        fout.close();
}

// Gauss-Laguerre
void save_gaulag(int n, const Doub alf) {
	VecDoub_O x(n);
	VecDoub_O w(n);

	for (int i=0;i<n;i++) {
                x[i] = 0.;
                w[i] = 0.;
        }

        gaulag(x, w, alf);

        char outname[50];
        sprintf(outname, "input/gaulag_weights_%.0lf_n_%d.txt", alf, n);

        ofstream fout(abs_path+outname);
        fout << "# Abscissas and weights for Gauss-Laguerre integration\n";

	fout << std::scientific;
	fout << std::setprecision(16);
	for (int i=0;i<n;i++) {
		fout << x[i] << "\t"  << w[i] << "\n";
	}
        fout.close();
}

void read_gleg_wgts(int n, double *x, double *w) {
  	char gaussname[50];
	string dummyLine;
        sprintf(gaussname, "input/gauleg_weights_n_%d.txt", n);

        // open input file
        std::fstream fin(abs_path+gaussname);

        if (!fin) {
                cout << "File regarding Gauleg quadrature not found: n = " << n << "\n";
                exit (EXIT_FAILURE);
        }

        getline(fin, dummyLine); //skip header line

        for (int i=0;i<n;i++) {
                std::stringstream ss(dummyLine);
                fin >> x[i] >> w[i];
	}

	fin.close();
}

// Gauss-Laguerre
void read_glag_wgts(int n, double *x, double *w, double alpha) {
        char gaussname[50];
        string dummyLine;

        sprintf(gaussname, "input/gaulag_weights_%d_n_%d.txt", int(alpha), n);

	// open input file
        std::fstream fin(abs_path+gaussname);

        if (!fin) {
                cout << "File regarding Gaulag quadrature not found: n = " << n << ", alpha = " << alpha << "\n";
                exit (EXIT_FAILURE);
        }

        getline(fin, dummyLine); //skip header line

        for (int i=0;i<n;i++) {
                std::stringstream ss(dummyLine);
                fin >> x[i] >> w[i];
        }

        fin.close();
}