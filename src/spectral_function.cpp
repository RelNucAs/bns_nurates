// \file  spectral_function.cpp
// \brief Define energy-dependent functions for emissivity/absorptivity
//        of neutrino-matter interactions

#include <tuple>

#include "spectral_function.hpp" // Header file for energy-dependent functions
#include "physics/nu_abs_em.hpp" // Header file for nu_abs_em rates

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

