#include <iostream>
#include <ostream>

#include "physics/NN_brem.hpp"

int main(){
	const double w = 10; //Neutrino energy [MeV]
	const double wp = 20; //Antineutrino energy [MeV]
	const double rho = 0.3730E+15; //Matter density [g/cm3]
	const double T = 1.204000e+01; //Temperature [MeV]
	const double xn = 6.866000e-01 ; //free neutron mass fraction
	const double xp = 3.134000e-01; //free proton mass fraction
	
	std::tuple<double,double> brem_kernel = bremcal_bruenn(w,wp,rho,T,xn,xp);

	cout << "w [MeV], wp [MeV], rho [g/cm3], T [MeV], Xn, Xp, em [1/MeV], abs [1/MeV]" << endl;
	cout << w << "\t" << wp << "\t" << rho << "\t" << T << "\t" << xn << "\t" << xp << "\t";
	cout << std::scientific << std::get<0>(brem_kernel) << "\t" << std::get<1>(brem_kernel) << endl;

}
