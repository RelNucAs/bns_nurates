#include "complete_FG.h"
#include "find_eta.h"

double n_part(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 0.5, w_leg, w_lag);
        f2 = compute_res(eta, theta, 1.5, w_leg, w_lag);

        return  K*pow(theta,1.5)*(f1+theta*f2);
}


double n_antipart(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 0.5, w_leg, w_lag);
        f2 = compute_res(-eta-2./theta, theta, 1.5, w_leg, w_lag);

        return  K*pow(theta,1.5)*(f1+theta*f2);
}

double P_part(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
	double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 1.5, w_leg, w_lag);
        f2 = compute_res(eta, theta, 2.5, w_leg, w_lag);

        return  K/3.*mLep*MeV*pow(theta,2.5)*(2.*f1+theta*f2); //erg/cm^3
}

double P_antipart(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 1.5, w_leg, w_lag);
        f2 = compute_res(-eta-2./theta, theta, 2.5, w_leg, w_lag);

        return  K/3.*mLep*MeV*pow(theta,2.5)*(2.*f1+theta*f2); //erg/cm^3
}


double e_part(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 1.5, w_leg, w_lag);
        f2 = compute_res(eta, theta, 2.5, w_leg, w_lag);

        return  (K/(mb*nLep))*mLep*MeV*pow(theta,2.5)*(f1+theta*f2); //erg/g
}

double e_antipart(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 1.5, w_leg, w_lag);
        f2 = compute_res(-eta-2./theta, theta, 2.5, w_leg, w_lag);

        return  (K/(mb*nLep))*mLep*MeV*pow(theta,2.5)*(f1+theta*f2) + 2.*mLep*MeV*n_antipart(nLep,T,eta,mLep,alpha,w_leg,w_lag)/(mb*nLep); //erg/g
}


double s_part(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 0.5, w_leg, w_lag);
        f2 = compute_res(eta, theta, 1.5, w_leg, w_lag);
        f3 = compute_res(eta, theta, 2.5, w_leg, w_lag);

        return  (K/(mb*nLep))*kB*MeV*pow(theta,1.5)*(-eta*f1+(5./3.-eta*theta)*f2+4./3.*theta*f3); //erg/K/g
}

double s_antipart(double nLep, double T, double eta, double mLep, gsl_integration_fixed_workspace *w_leg, gsl_integration_fixed_workspace *w_lag) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 0.5, w_leg, w_lag);
        f2 = compute_res(-eta-2./theta, theta, 1.5, w_leg, w_lag);
        f3 = compute_res(-eta-2./theta, theta, 2.5, w_leg, w_lag);

	eta_anti = -eta-2./theta;

        return  (K/(mb*nLep))*kB*MeV*pow(theta,1.5)*(-eta_anti*f1+(5./3.-eta_anti*theta)*f2+4./3.*theta*f3); //erg/K/g
}

