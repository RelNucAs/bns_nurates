#pragma once

#include <fstream>

#include "./tools/fermi_integrals.h"
#include "./tools/nr3.h"
#include "./tools/gamma.h"
#include "./tools/gauss_wgts.h"
#include "./spectral_function.h"

void save_gauleg(int n, const Doub x1, const Doub x2) {
        VecDoub_O x(n);
        VecDoub_O w(n);

	char fileHeader[100];

        for (int i=0;i<n;i++) {
                x[i] = 0.;
                w[i] = 0.;
        }

	gauleg(x1, x2, x, w);

        char outname[50];
        sprintf(outname, "input/gauleg_weights_n_%d.txt", n);

        ofstream fout(abs_path+outname);
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

        char outname[50];
        sprintf(outname, "input/gaulag_weights_%.0lf_n_%d.txt", alf, n);

        ofstream fout(abs_path+outname);
        fout << "# Abscissas and weights for Gauss-Laguerre integration\n";
	
	for (int i=0;i<n;i++) {
		fout << std::scientific << x[i] << "\t"  << w[i] << "\n";
	}
        fout.close();
}

void read_gleg_wgts(int n, double *x, double *w) {
  	char gaussname[50];
	string dummyLine;
        sprintf(gaussname, "input/gauleg_weights_n_%d.txt", n);

	//if (g_type == 0) {
        //        sprintf(gaussname, "../input/gauleg_weights_n_%d.txt", n);
        //} else if (g_type == 1) {
        //        sprintf(gaussname, "../input/gaulag_weights_0_n_%d.txt", n);
        //} else if (g_type == 2) {
        //        sprintf(gaussname, "../input/gaulag_weights_5_n_%d.txt", n);
        //} else {
        //        printf("Index of Gaussian Quadrature method not recognized");
        //        exit (EXIT_FAILURE);
        //}

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


double do_integration(int n, double *wgt, double *spectrArray) {
	
	double integral = 0.;

	for (int i=0;i<n;i++) {
		//printf("%.3lf, %.3lf\n",wgt[i], spectrArray[i]);
		integral += wgt[i] * spectrArray[i];
        }

	return integral;

}

double glag_integ(int n, my_function F, double alpha) {
        double x[n], w[n];
        double f[n];

        read_glag_wgts(n, x, w, alpha);

	if (alpha == 0.) {
		for (int i=0;i<n;i++) f[i] = F.function(x[i], F.params) * exp(x[i]);
	} else {
		for (int i=0;i<n;i++) f[i] = F.function(x[i], F.params) * exp(x[i]) / pow(x[i],alpha);
	}

        return do_integration(n, w, f);
}


double split_value(double T, double eta_e) {
	return T * Fermi_integral_p5(eta_e) / Fermi_integral_p4(eta_e);
}


double gleg_integ(int n, my_function F, double t) {
        double x[n], w[n];
        double f1[n], f2[n];

        read_gleg_wgts(n, x, w);

        for (int i=0;i<n;i++) {
                f1[i] = F.function(t*x[i], F.params);
                f2[i] = F.function(t/x[i], F.params) / (x[i]*x[i]);
        }
        return t * (do_integration(n, w, f1) + do_integration(n, w, f2));
}


double GSL_integration(int n, double a, double b, double alpha, double beta, gsl_function F, const gsl_integration_fixed_type * gauss_type) {
        gsl_integration_fixed_workspace *w;
        double result;
        w = gsl_integration_fixed_alloc(gauss_type, n, a, b, alpha, beta);
        gsl_integration_fixed(&F, &result, w);
        gsl_integration_fixed_free (w);
        return result;
}

double find_max(int nslice, double x0, double x1, double (*func)(double, struct my_f_params), struct my_f_params p) {
        double val, max = 0;
        double x_step, guess;
        double dx = (log10(x1)-log10(x0)) / (double) nslice;
        for (double x=x0;x<x1;x+=dx) {
                x_step = pow(10.,x);
                val = func(x_step, p);
                //printf("%.3e  ",val);
                if (val > max) {
                        max = val;
                        guess = x_step;
                }
        }
        //printf("\n");
        return guess;
}

double find_guess(int nslice, my_function F) {
        double val, max = 0;
        double guess;
        void *p = F.params;
        struct FD_params * params = (struct FD_params *)p;
        double eta = (params->eta);
        double x0, x1;

        if (abs(eta) < 10.) {
                x0 = 0.1*abs(eta),
                x1 = 100*abs(eta);
        } else {
                x0 = 0.5*abs(eta);
                x1 = 2.0*abs(eta);
        }

        double dx = (x1-x0) / (double) nslice;
        for (double x=x0;x<x1;x+=dx) {
                val = F.function(x,p);
                //printf("x = %.3e, Val = %.3e\n", x, val);
                if (val > max) {
                        max = val;
                        guess = x;
                }
        }
        return guess;
}

