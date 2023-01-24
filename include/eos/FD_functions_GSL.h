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

double fermi_func1_GSL(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
	float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
        double a = (params->a);
        double b = (params->b);
	
	double y = (b-a)*x+a;
	double tmp;

	tmp = y*y-eta;
        tmp = safe_FDexp(tmp);

	return (b-a)*2.*pow(y,(2.*k+1.))*pow((1.+y*y*theta/2.),0.5)*tmp;
}


double fermi_func2_GSL(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
        double a = (params->a);
        double b = (params->b);
        
	double y = (b-a)*x+a;
	double tmp; 
         
        tmp = y-eta;
        tmp = safe_FDexp(tmp);

        return (b-a)*pow(y,k)*pow((1.+y*theta/2.),0.5)*tmp;
}

double fermi_func3_GSL(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
	double s3 = (params->a);
	double tmp = -x-s3+eta;
	
	tmp = safe_FDexp(tmp);
	
	double f3 = exp(eta-s3)*pow(x+s3,k)*pow((1.+(x+s3)*theta/2.),0.5) * tmp; // / (1.+exp(-x-s3+eta));
	
	if (alpha != 0.) f3 = f3 / pow(x,alpha);

	return f3; 
}

double fermi_ed1_GSL(double x, void *p) {
	struct FD_params * params = (struct FD_params *)p;
        float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
        double tmp;
        double a = (params->a);
        double b = (params->b);

        double y = (b-a)*x+a;

	tmp = y*y-eta;
        tmp = std::min(tmp,120.);
        tmp = exp(tmp)/pow(1.+exp(tmp),2.);

	return (b-a)*2.*pow(y,(2.*k+1.))*pow((1.+theta*y*y/2.),0.5)*tmp;
}


double fermi_ed2_GSL(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
        double tmp;
        double a = (params->a);
        double b = (params->b);

        double y = (b-a)*x+a;

        tmp = y-eta;
        tmp = std::min(tmp,120.);
        tmp = exp(tmp)/pow(1.+exp(tmp),2.);

        return (b-a)*pow(y,k)*pow((1.+theta*y/2.),0.5)*tmp;
}


double fermi_ed3_GSL(double x, void *p) {
        struct FD_params * params = (struct FD_params *)p;
        float k = (params->k);
        double eta = (params->eta);
        double theta = (params->theta);
        double s3 = (params->a);
        double tmp,tmp1,tmp2;
	double f3;

        tmp1 = 2.*x+s3-eta;
        tmp1 = std::min(tmp1,120.);
        tmp2 = x+s3-eta;
        tmp2 = std::min(tmp2,120.);
        tmp  = exp(tmp1)/pow((1.+exp(tmp2)),2.);

	f3 = pow((x+s3),k)*pow((1.+(x+s3)*theta/2.),0.5)*tmp;

	if (alpha != 0.) f3 = f3/pow(x,alpha);
	
	return f3;
}


