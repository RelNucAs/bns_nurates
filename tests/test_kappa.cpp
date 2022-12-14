#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <string> 
#include <vector>
#include <sstream>
#include <cmath> //math library for basic operations, pi
#include <cstdio> //libabry for printf 
#include <cstdlib> //libabry for abs()

#include "../source/constants.h" //Header file containing all relevant constants
#include "../source/parameters.h" //Header file containing all parameters
#include "../source/weak_magnetism.h" //Header file containing functions to compute WM
#include "../source/nu_abs_em.h" //Header file containing all rates
#include "../source/nu_elastic_scatt.h" //Header file containing elastic scattering rates
#include "../source/spectral_function.h"
#include "../source/tools/fermi_integrals.h"
#include "../source/tools/NewtonRaphson.h"
#include "../source/tools/nr3.h"
#include "../source/mnewt.h"

using namespace constants;
using namespace parameters;
using namespace nuabsem;
using namespace elastic_scatt;
using namespace weakmag;
using namespace std;

int main (){
	int zone;
        double r, m, d, T, ye, s, e, ynu, ynubar, ynum, znu, znubar, znum;
	double nb;
	double A, B, n_nu, j_nu;
	double mu_e, mu_hat, mu_nu;
	double x0[2], x[2];
	string finLine, dinLine;
	gsl_multiroot_function_fdf FDF;
	gsl_multiroot_function F;
	double tmp[2];
	int opt = 0; //0:thick, 1:thin

	std::fstream din("../ccsn_data/abun1d00380.d");
        if (!din)
        {
                cout << "File not found" << '\n';
        }
        for (int i = 0; i < 4; i++) { //skip first two lines and read first one with data
                getline(din, dinLine);
        }


        // open input profile
        std::fstream fin("../input/nurates_1.008E+01.txt"); //read file
        if (!fin)
        {
                cout << "File not found" << '\n';
        }
        for (int i = 0; i < 3; i++) { //skip first two lines and read first one with data
                getline(fin, finLine);
        }
	
	int i = 0;

        while (!din.fail()){
                // read data from input table
                std::stringstream sstr(dinLine);
                std::stringstream ss(finLine);

                sstr >> zone >> r >> m >> d >> T >> ye >> s >> e;
	        sstr >> ynu >> ynubar >> ynum >> znu >> znubar >> znum;
		
		//printf("zone = %d, r = %.3e, m = %.3e, d = %.3e, T = %.3e, Ye = %.3lf, s = %.3e, e = %.3e, ", zone,r,m,d,T,ye,s,e);
		//printf("ynu = %.3e, ynubar = %.3e, ynum = %.3e, znu = %.3e, znubar = %.3e, znum = %.3e\n", ynu,ynubar,ynum,znu,znubar,znum);

    		ss >> zone >> r >> d >> T >> ye >> mu_e >> mu_hat;
		
		//printf("zone = %d, r = %.3e, d = %.3e, T = %.3e, Ye = %.3lf, mu_e = %.3e, mu_hat = %.3e\n", zone,r,d,T,ye,mu_e,mu_hat);
                //printf("%s\n", dinLine.c_str());
		nb = d/mb;

		mu_nu = mu_e - (Q+mu_hat);
		
		//printf("First = %.3e, second = %.3e\n", 4*pi*Fermi_integral_p2(mu_nu/T)*pow(T,3.)/pow(h*c,3.), 4*pi*Fermi_integral_p3(mu_nu/T)*pow(T,4.)/pow(h*c,3.));
		//printf("%.3e, %.3e\n", ynu*nb, znu*nb);

		A = 1./T;
		B = mu_nu/T;

		//// Check NR correctness
		//printf("Check Newton-Raphson correctness\n");

		if (opt == 0) { //optically thick
			n_nu = 4.*pi*Fermi_integral_p2(B)/pow(A,3.);
			j_nu = 4.*pi*Fermi_integral_p3(B)/pow(A,4.);
		} else if (opt == 1) { //optically thin
			n_nu = 4.*pi*gsl_sf_gamma(A+3.)/pow(B,A+3.);
			j_nu = 4.*pi*gsl_sf_gamma(A+4.)/pow(B,A+4.);
		}
		
		double C[2] = {n_nu/pow(h*c,3.),j_nu/pow(h*c,3.)};
	        x[0] = A * 1.02;
        	x[1] = B * 0.98;

		x0[0] = A*1.02;
		x0[1] = B*0.98;
		
		//printf("Starting point: x0 = %.5lf, x1 = %.5lf\n", x[0], x[1]);
		
		struct M1_values p = {C[0],C[1]};
		
		if (opt == 0) {
			FDF = {&thick_f, &thick_df, &thick_fdf_GSL, 2, &p};
			F = {&thick_f, 2, &p};
                 	mnewt2d(x,C,&thick_fdf);
	 	} else if (opt == 1) {
			FDF = {&thin_f, &thin_df, &thin_fdf_GSL, 2, &p};
			F = {&thin_f, 2, &p};
                 	mnewt2d(x,C,&thin_fdf);
		}

		//printf("NR function:    A = %.5lf, B = %.5lf\n", x[0], x[1]);

	        NRaphson_2D(x0, x, FDF);
		
		//printf("GSL function:   A = %.5lf, B = %.5lf\n\n", x[0], x[1]);
	        //NRaphson_2D_findiff(x0, x, F);
		
		//if (abs(A-x[0])/min(A,x[0]) > 1.e-3) {
		//	printf("Exact A = %.3lf,	Recovered A = %.3lf,	diff = %.3e\n",   A, x[0], abs(A-x[0])/min(A,x[0]));
		//} else {
		//	printf("Exact A = %.3lf,	Recovered A = %.3lf\n",   A, x[0]);
		//}

		//if (abs(B-x[1])/min(B,x[1]) > 1.e-3) {
		//	printf("Exact B = %.3lf,	Recovered B = %.3lf,	diff = %.3e\n\n", B, x[1], abs(B-x[1])/min(B,x[1]));
		//} else {
		//	printf("Exact B = %.3lf,	Recovered B = %.3lf\n\n", B, x[1]);
		//}
		
		//// Compute A e B (T and eta_nu) along CCSN 1D profile
		C[0] = ynu*nb;
		C[1] = znu*nb;
		p = {ynu*nb, znu*nb};


		if (i == 0) {
			x0[0] = A*0.98;
        		x0[1] = B*1.02;
			x[0] = x0[0];
			x[1] = x0[1];
		} else {
		        x0[0] = A*0.98;
                        x0[1] = B*1.02;
			//x0[0] = tmp[0];
			//x0[1] = tmp[1];
			x[0] = x0[0];
			x[1] = x0[1];
		}

		if (opt == 0) {
			FDF = {&thick_f, &thick_df, &thick_fdf_GSL, 2, &p};
			F = {&thick_f, 2, &p};
			mnewt2d(x,C,&thick_fdf);
		} else if (opt == 1) {
			FDF = {&thin_f, &thin_df, &thin_fdf_GSL, 2, &p};
			F = {&thin_f, 2, &p};
			mnewt2d(x,C,&thin_fdf);
		}


	        NRaphson_2D(x0, tmp, FDF);
	        //NRaphson_2D_findiff(x0, x, F);
	
		//tmp[0] = x[0];
		//tmp[1] = x[1];

		//x0[0] = A;
		//x0[1] = B;

	        //NRaphson_2D(x0, x, FDF);
		
		//printf("%.5e, %.5e, %.5e, %.5e, %.5e\n", r, d, mu_nu, x[1]/x[0], tmp[1]/tmp[0]);
		printf("%.5e, %.5e, %.5e, %.5e, %.5e\n",   r, d, T, x[0], x[1]);
		printf("%.5e, %.5e, %.5e, %.5e, %.5e\n\n", r, d, T, tmp[0], tmp[1]);
		//if (abs(x[0]-tmp[0]) > 1.e-4) printf("r = %.3e\n", r);
		//if (abs(x[1]-tmp[1]) > 1.e-4) printf("r = %.3e\n", r);
                std::getline(fin, finLine);
                std::getline(din, dinLine);
		i += 1;
        }


	din.close();
	fin.close();

	return 0;
}

