#pragma once

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_multiroots.h>
#include "./tools/fermi_integrals.h"
#include "./nu_abs_em.h" //Header file containing all rates

using namespace nuabsem;

struct FD_params { int k; double eta; double alpha; };
struct function_struct
{
  double (* function) (double x, void * params);
  void * params;
};
typedef struct function_struct my_function;

struct my_f_params { double nb; double T; double mL; double yp; double yn; double mu_L; double mu_hat; double dU; };
struct M1_values { double n_nu; double J_nu; };
struct fnu_params {double A; double B; };
struct K_params {struct my_f_params rate_pars; struct fnu_params f_pars; };

double j0_nue(double x, void *p) {
        struct my_f_params * pars = (struct my_f_params *)p;
        double fEm, fAb;
        std::tie(fEm,fAb) = nu_n_abs(x, pars->nb, pars->T, pars->mL, pars->yp, pars->yn, pars->mu_L, pars->mu_hat, pars->dU);
        return 4. * pi * pow(x,2.) * fEm;
}


double j0_anue(double x, void *p) {
        struct my_f_params * pars = (struct my_f_params *)p;
        double fEm, fAb;
        std::tie(fEm,fAb) = nu_p_abs(x, pars->nb, pars->T, pars->mL, pars->yp, pars->yn, pars->mu_L, pars->mu_hat, pars->dU);
        return 4. * pi * pow(x,2.) * fEm;
}


double j_nue(double x, void *p) {
        //struct my_f_params * pars = (struct my_f_params *)p;
        //double fEm, fAb;
        //std::tie(fEm,fAb) = nu_n_abs(x, pars->nb, pars->T, pars->mL, pars->yp, pars->yn, pars->mu_L, pars->mu_hat, pars->dU);
        //return 4. * pi * pow(x,3.) * fEm;
	return j0_nue(x,p)*x;
}

double j_anue(double x, void *p) {
        //struct my_f_params * pars = (struct my_f_params *)p;
        //double fEm, fAb;
        //std::tie(fEm,fAb) = nu_p_abs(x, pars->nb, pars->T, pars->mL, pars->yp, pars->yn, pars->mu_L, pars->mu_hat, pars->dU);
        //return 4. * pi * pow(x,3.) * fEm;
	return j0_anue(x,p)*x;
}

/*	OPTICALLY THICK CONDITIONS	*/
double k0_thick_nue(double x, void *p) {
	struct K_params * pars = (struct K_params *)p;
        struct my_f_params p1 = pars->rate_pars;
        struct fnu_params  p2 = pars->f_pars;

        double fnu = 1.+exp(p2.A*x+p2.B);
        double fEm, fAb;
        std::tie(fEm,fAb) = nu_n_abs(x, p1.nb, p1.T, p1.mL, p1.yp, p1.yn, p1.mu_L, p1.mu_hat, p1.dU);

        return 4. * pi * pow(x,2.) * fAb / fnu;
}

double k_thick_nue(double x, void *p) {
        return k0_thick_nue(x,p)*x;
}

double k0_thick_anue(double x, void *p) {
        struct K_params * pars = (struct K_params *)p;
        struct my_f_params p1 = pars->rate_pars;
        struct fnu_params  p2 = pars->f_pars;
	
        double fnu = 1.+exp(p2.A*x+p2.B);
        double fEm, fAb;
        std::tie(fEm,fAb) = nu_p_abs(x, p1.nb, p1.T, p1.mL, p1.yp, p1.yn, p1.mu_L, p1.mu_hat, p1.dU);
        
	return 4. * pi * pow(x,2.) * fAb / fnu;
}

double k_thick_anue(double x, void *p) {
        return k0_thick_anue(x,p)*x;
}


/*	OPTICALLY THIN CONDITIONS	*/
double k0_thin_nue(double x, void *p) {
	struct K_params * pars = (struct K_params *)p;
        struct my_f_params p1 = pars->rate_pars;
        struct fnu_params  p2 = pars->f_pars;

        double fEm, fAb;
        std::tie(fEm,fAb) = nu_n_abs(x, p1.nb, p1.T, p1.mL, p1.yp, p1.yn, p1.mu_L, p1.mu_hat, p1.dU);

        return 4. * pi * pow(x,p2.A+2.) * fAb * exp(-p2.B*x);
}

	
double k_thin_nue(double x, void *p) {
	return k0_thin_nue(x,p)*x;
}


double k0_thin_anue(double x, void *p) {
        struct K_params * pars = (struct K_params *)p;
        struct my_f_params p1 = pars->rate_pars;
        struct fnu_params  p2 = pars->f_pars;

        double fEm, fAb;
        std::tie(fEm,fAb) = nu_p_abs(x, p1.nb, p1.T, p1.mL, p1.yp, p1.yn, p1.mu_L, p1.mu_hat, p1.dU);

        return 4. * pi * pow(x,p2.A+2.) * fAb * exp(-p2.B*x);
}

double k_thin_anue(double x, void *p) {
        return k0_thin_anue(x,p)*x;
}


double k0_thin_anue_old(double x, struct my_f_params p1, struct fnu_params p2) {
        double nb = (p1.nb);
        double T  = (p1.T);
        double me = (p1.mL);
        double yp = (p1.yp);
        double yn = (p1.yn);
        double mu_e = (p1.mu_L);
        double mu_hat = (p1.mu_hat);
        double dU = (p1.dU);

        double A = (p2.A);
        double B = (p2.B);

        double fEm, fAb;

        std::tie(fEm,fAb) = nu_p_abs(x, nb, T, me, yp, yn, mu_e, mu_hat, dU);

        return 4. * pi * pow(x,A+2.) * fAb * exp(-B*x);
}


int thick_f(const gsl_vector *x, void *p, gsl_vector *f) {
        struct M1_values * params = (struct M1_values *) p;
        const double n_nu = (params->n_nu);
        const double J_nu = (params->J_nu);

        const double x0 = gsl_vector_get(x,0);
        const double x1 = gsl_vector_get(x,1);

        gsl_vector_set(f, 0, n_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p2(x1)/pow(x0,3.));
        gsl_vector_set(f, 1, J_nu*pow(h*c,3.) - 4.*pi*Fermi_integral_p3(x1)/pow(x0,4.));

        return GSL_SUCCESS;
}


int thick_df(const gsl_vector *x, void *p, gsl_matrix *J) {
        struct M1_values * params = (struct M1_values *) p;

        const double x0 = gsl_vector_get(x,0);
        const double x1 = gsl_vector_get(x,1);

        gsl_matrix_set(J, 0, 0,  3. * 4.*pi*Fermi_integral_p2(x1)/pow(x0,4.));
        gsl_matrix_set(J, 0, 1, -2. * 4.*pi*Fermi_integral_p1(x1)/pow(x0,3.));
        gsl_matrix_set(J, 1, 0,  4. * 4.*pi*Fermi_integral_p3(x1)/pow(x0,5.));
        gsl_matrix_set(J, 1, 1, -3. * 4.*pi*Fermi_integral_p2(x1)/pow(x0,4.));

        return GSL_SUCCESS;
}

int thick_fdf_old (gsl_vector *x, void *p, gsl_vector *f, gsl_matrix *J) {
        struct M1_values * params = (struct M1_values *)p;
        const double n_nu = (params->n_nu);
        const double J_nu = (params->J_nu);

        const double x0 = gsl_vector_get(x,0);
        const double x1 = gsl_vector_get(x,1);

        gsl_vector_set(f, 0, n_nu - 4.*pi*Fermi_integral_p2(x1)/pow(x0,3.));
        gsl_vector_set(f, 1, J_nu - 4.*pi*Fermi_integral_p3(x1)/pow(x0,4.));

        gsl_matrix_set(J, 0, 0,  3. * 4.*pi*Fermi_integral_p2(x1)/pow(x0,4.));
        gsl_matrix_set(J, 0, 1, -2. * 4.*pi*Fermi_integral_p1(x1)/pow(x0,3.));
        gsl_matrix_set(J, 1, 0,  4. * 4.*pi*Fermi_integral_p3(x1)/pow(x0,5.));
        gsl_matrix_set(J, 1, 1, -3. * 4.*pi*Fermi_integral_p2(x1)/pow(x0,4.));

        return GSL_SUCCESS;
}

int thick_fdf_GSL(const gsl_vector * x, void *params, gsl_vector *f, gsl_matrix *J) {
        thick_f(x, params, f);
        thick_df(x, params, J);
        return GSL_SUCCESS;
}


int thin_f(const gsl_vector *x, void *p, gsl_vector *f) {
        struct M1_values * params = (struct M1_values *)p;
        const double n_nu = (params->n_nu);
        const double J_nu = (params->J_nu);


        const double x0 = gsl_vector_get(x,0);
        const double x1 = gsl_vector_get(x,1);

        gsl_vector_set(f, 0, n_nu*pow(h*c,3.) - 4.*pi*gsl_sf_gamma(x0+3.)/pow(x1,x0+3.));
        gsl_vector_set(f, 1, J_nu*pow(h*c,3.) - 4.*pi*gsl_sf_gamma(x0+4.)/pow(x1,x0+4.));

        return GSL_SUCCESS;
}

int thin_df(const gsl_vector *x, void *p, gsl_matrix *J) {
        struct M1_values * params = (struct M1_values *) p;

        const double x0 = gsl_vector_get(x,0);
        const double x1 = gsl_vector_get(x,1);

        gsl_matrix_set(J, 0, 0, 4.*pi*gsl_sf_gamma(x0+3.)/pow(x1,x0+3.)*(log(x1)-gsl_sf_psi(x0+3.)));
        gsl_matrix_set(J, 0, 1, 4.*pi*(x0+3.)*gsl_sf_gamma(x0+3.)/pow(x1,x0+4.));
        gsl_matrix_set(J, 1, 0, 4.*pi*gsl_sf_gamma(x0+4.)/pow(x1,x0+4.)*(log(x1)-gsl_sf_psi(x0+4.)));
        gsl_matrix_set(J, 1, 1, 4.*pi*(x0+4.)*gsl_sf_gamma(x0+4.)/pow(x1,x0+5.));

        return GSL_SUCCESS;
}

int thin_fdf_GSL(const gsl_vector * x, void *params, gsl_vector *f, gsl_matrix *J) {
        thin_f(x, params, f);
        thin_df(x, params, J);
        return GSL_SUCCESS;
}

double FD_laguerre(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        int k = (params->k);
        double eta = (params->eta);
        double alpha = (params->alpha);
        return exp(x) * pow(x,-alpha) * pow(x,k) / (exp(x-eta)+1.);
}

double FD_legendre(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        int k = (params->k);
        double eta = (params->eta);
        return pow(x,k) / (exp(x-eta)+1.);
}

double f_NRaphson(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        double eta = (params->eta);
        return 5.*exp(eta) - exp(x)*(x-5.);
}

double df_NRaphson(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        double eta = (params->eta);
        return exp(x)*(4.-x);
}

void fdf_NRaphson(double x, void *p, double * f, double * df) {
        struct FD_params * params = (struct FD_params *)p;
        double eta = (params->eta);
        double t  = 5.*exp(eta) - exp(x)*(x-5.);
        double dt = exp(x)*(4.-x);
        *f = t;
        *df = dt;   /* uses existing value */
}

