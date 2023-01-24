/* NEUTRINO EOS AT EQUILIBRIUM WITH RADIATION (neutrinos degeneracy parameter is equal to minus anti-neutrinos degeneracy parameter) */
#include "../constants_bis.h"
#include "../tools/fermi_integrals.h"

using namespace constants;

//Input parameters are matter density rho [g/cm^3], temperature t [MeV] and neutrinos
//degeneracy parameter eta. The density rho is meant as the baryon number density 
//multiplied by the atomic mass unit

//fraction of neutrinos (yn)
double nu_fraction(double rho, double T, double eta) {
	return (4.*pi*mu/rho) * pow(T/(h*c),3.) * Fermi_integral_p2(eta);
}

//fraction of antineutrinos (yan)
double anu_fraction(double rho, double T, double eta) {
	return (4.*pi*mu/rho) * pow(T/(h*c),3.) * Fermi_integral_p2(-eta);
}

//energy per baryon of neutrinos (en)
double nu_energy(double rho, double T, double eta) {
        return (4.*pi*mu/rho) * pow(T/(h*c),3.) * T * Fermi_integral_p3(eta); //MeV/baryon
}

//energy per baryon of antineutrinos (ean)
double anu_energy(double rho, double T, double eta) {
        return (4.*pi*mu/rho) * pow(T/(h*c),3.) * T * Fermi_integral_p3(-eta); //MeV/baryon
}
