#include "complete_FG.h"

double n_part(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 0.5);
        f2 = compute_res(eta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1+theta*f2);
}


double n_antipart(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1+theta*f2);
}

double P_part(double T, double eta, double mLep) {
	double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 1.5);
        f2 = compute_res(eta, theta, 2.5);

        return  K/3.*mLep*MeV*pow(theta,2.5)*(2.*f1+theta*f2); //erg/cm^3
}

double P_antipart(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 1.5);
        f2 = compute_res(-eta-2./theta, theta, 2.5);

        return  K/3.*mLep*MeV*pow(theta,2.5)*(2.*f1+theta*f2); //erg/cm^3
}


double e_part(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 0.5);
        f2 = compute_res(eta, theta, 1.5);
        f3 = compute_res(eta, theta, 2.5);

        return  K*mLep*MeV*pow(theta,1.5)*(f1+2.*theta*f2+theta*theta*f3); //erg/cm^3
}

double e_antipart(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 1.5);
        f3 = compute_res(-eta-2./theta, theta, 2.5);

        //return  (K/(mb*nLep))*mLep*MeV*pow(theta,2.5)*(f1+theta*f2) + 2.*mLep*MeV*n_antipart(nLep,T,eta,mLep,alpha,w_leg,w_lag)/(mb*nLep); //erg/g
        return  K*mLep*MeV*pow(theta,1.5)*(f1+2.*theta*f2+theta*theta*f3); //erg/cm^3
}


double s_part(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(eta, theta, 0.5);
        f2 = compute_res(eta, theta, 1.5);
        f3 = compute_res(eta, theta, 2.5);

        return  K*kB*MeV*pow(theta,1.5)*(-eta*f1+(5./3.-eta*theta)*f2+4./3.*theta*f3); //erg/K/cm^3
}

double s_antipart(double T, double eta, double mLep) {
        double K = 8*sqrt(2)*pi*pow(mLep/(h*c),3.);
        double f1, f2, f3;
        double theta = T/mLep;

        f1 = compute_res(-eta-2./theta, theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 1.5);
        f3 = compute_res(-eta-2./theta, theta, 2.5);

	double eta_anti = -eta-2./theta;

        return  K*kB*MeV*pow(theta,1.5)*(-eta_anti*f1+(5./3.-eta_anti*theta)*f2+4./3.*theta*f3); //erg/K/cm^3
}

