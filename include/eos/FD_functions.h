#include <algorithm>
#include "../parameters.h"
#include "./gleg.h"
#include "./glag.h"
using namespace parameters;

struct FD_params { float k; double eta; double theta; double a; double b; };

double safe_FDexp(double x) {
	double out;
	double tmp = std::min(x,120.);
	if (tmp<20.){
		out = 1. - exp(tmp)/(1.+exp(tmp));	
	} else {
		out = pow(1.+exp(tmp),-1.);
	}
	return out;
}

void fermi_func1(double eta, double theta, float k, double a, double b, double f[]){
	double x, tmp;
	for (int i=0;i<ngle;i++) {
		x = (b-a)*xgle[i]+a;

		tmp = x*x-eta;
		tmp = safe_FDexp(tmp);
		f[i] = (b-a)*2.*pow(x,(2.*k+1.))*pow((1.+x*x*theta/2.),0.5)*tmp;
		//printf("x = %.10e, f = %.10e\n", x, f[i]);
	}
	//printf("\n");
	return;
}

void fermi_func2(double eta, double theta, float k, double a, double b, double f[]){
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

                tmp = x-eta;
		tmp = safe_FDexp(tmp);

                f[i] = (b-a)*pow(x,k)*pow((1.+x*theta/2.),0.5)*tmp;
		//printf("x = %.10e, f = %.10e\n", x, f[i]);
        }
	//printf("\n");
        return;
}
         


void fermi_func3(double eta, double theta, float k, double s3, double f[]){
	double x, tmp;
	for (int i=0;i<ngla;i++) {
		x = xgla[i];

		tmp = -x-s3+eta;
		tmp = safe_FDexp(tmp);

		f[i] = exp(eta-s3)*pow(x+s3,k)*pow((1.+(x+s3)*theta/2.),0.5) * tmp; /// (1.+exp(-x-s3+eta));
	}

	if (alpha != 0.) {
		for (int i=0;i<ngla;i++) f[i] = f[i]/pow(xgla[i],alpha);
	}
	return; 
}


void fermi_ed1(double eta, double theta, float k, double a, double b, double f[]) {
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

                tmp = x*x-eta;
                tmp = std::min(tmp,120.);
        	tmp = exp(tmp)/pow(1.+exp(tmp),2.);

                f[i] = (b-a)*2.*pow(x,(2.*k+1.))*pow((1.+x*x*theta/2.),0.5)*tmp;
        }
        return;
}

void fermi_ed2(double eta, double theta, float k, double a, double b, double f[]) {
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

        	tmp = x-eta;
		tmp = std::min(tmp,120.);
		tmp = exp(tmp)/pow(1.+exp(tmp),2.);
		
		f[i] = (b-a)*pow(x,k)*pow((1.+theta*x/2.),0.5)*tmp;
	}
	return;
}


void fermi_ed3(double eta, double theta, float k, double s3, double f[]) {
        double x, tmp, tmp1, tmp2;
        for (int i=0;i<ngla;i++) {
                x = xgla[i];

		tmp1 = 2.*x+s3-eta;
		tmp1 = std::min(tmp1,120.);
		tmp2 = x+s3-eta;
		tmp2 = std::min(tmp2,120.);
		tmp  = exp(tmp1)/pow((1.+exp(tmp2)),2.);

		f[i] = pow((x+s3),k)*pow((1.+(x+s3)*theta/2.),0.5)*tmp;
	}

	if (alpha != 0.) {
		for (int i=0;i<ngla;i++) f[i] = f[i]/pow(xgla[i],alpha);
	}
	return;
}
