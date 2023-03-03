#ifndef SPEC_FUNCTION_H
#define SPEC_FUNCTION_H

// \file  spectral_function.hpp
// \brief Define energy-dependent functions for emissivity/absorptivity
//        of neutrino-matter interactions

struct function_struct
{
  double (* function) (double x, void * params);
  void * params;
};
typedef struct function_struct my_function;

struct FD_params { int k; double eta; double alpha; };
struct my_f_params { double nb; double T; double mL; double yp; double yn; double mu_L; double mu_hat; double dU; };
struct M1_values { double n_nu; double J_nu; };
struct fnu_params {double A; double B; };
struct K_params {struct my_f_params rate_pars; struct fnu_params f_pars; };

double j0_nue(double x, void *p);

double j0_anue(double x, void *p);

double j_nue(double x, void *p);

double j_anue(double x, void *p);

/*	OPTICALLY THICK CONDITIONS	*/
double k0_thick_nue(double x, void *p);

double k_thick_nue(double x, void *p);

double k0_thick_anue(double x, void *p);

double k_thick_anue(double x, void *p);

/*	OPTICALLY THIN CONDITIONS	*/
double k0_thin_nue(double x, void *p);
	
double k_thin_nue(double x, void *p);

double k0_thin_anue(double x, void *p);

double k_thin_anue(double x, void *p);

double k0_thin_anue_old(double x, struct my_f_params p1, struct fnu_params p2);

double FD_laguerre(double x, void *p);

double FD_legendre(double x, void *p);

double f_NRaphson(double x, void *p);

double df_NRaphson(double x, void *p);

void fdf_NRaphson(double x, void *p, double * f, double * df);

#endif
