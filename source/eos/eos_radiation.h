/* RADIATION EOS */
#include <cmath>
#include "../constants_bis.h"

using namespace constants;

double zeta3 = std::riemann_zeta(3.);

//number density
double n_rad(double T){ //T in MeV
      return 16.*pi*zeta3*pow(T/(h*c),3.); //cm^-3
}

//internal energy (per unit volume)
double e_rad(double T){ //T in MeV
	return (8.*pow(pi,5.)*MeV*pow(T,4.)) / (15*pow(h*c,3.)); //erg cm^-3	
}

//pressure 
double P_rad(double T){ //T in MeV
	return e_rad(T)/3.; //erg cm^-3
}
	
//entropy
double s_rad(double T){ //T in MeV
	return (32.*pow(pi,5.)/45.) * kB * MeV * pow(T/(h*c),3.); //erg cm^-3 K^-1
}

