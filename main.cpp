// This contains all the main functions which implement the various opactities due to neutral current interactions (elastic scattering)
// Based on Bruenn(1985) and Burrows et al(2004). List of interactions considered: 
// 1: electron neutrino absorption on neutrons
// 2: electron antineutrino absorption on protons

//Libraries to be loaded

#include <iostream> //Input&Output stream model based library
#include <cmath> //math library for basic operations, pi
#include <math.h> //math library for pow, trig
#include <stdio.h> //libabry for printf 

//Unit conversions

const double eV = 1.6*10e-12

//Global Constants (CGS)

const double h_bar = 1.05457266*(10e-27)
const double c = 29979245800
const double Gf = 1.436*(10e-49)
const double me = 9.10938356*(10e-28)
const double mn = 1.6749286*(10e-24)
const double mp = 1.6726231*(10e-24)
const double pi = 3.14159265358979323846
const double kb = 1.380658*(10e-16)

//Simplying constants
const double e_rm = me*pow(c,2)
const double n_rm = mn*pow(c,2)
const double p_rm = mp*pow(c,2)
const double sigma_o = (4*Gf*pow(Gf,2)*pow(e_rm,2))/(pi*pow(h_bar*c,4))
printf("Sigma_0",sigma_o);
const double gA = 1.23
const delta_np = 1.29332*eV*10e6

//corrections
//Fermi distributions

double Feqnu(double e_nu, double mu_e, double mu_np, double temp)
{
return pow(exp((e_nu-mu_e+mu_np)/(kb*temp))+1,-1)  //beta equilibrium assumed (mu_nu = mu_e-mu_np)
}

//blocking factor and stimulated absorption
double B_nu(double e_nu, double mu_e, double mu_np, double temp)
{
Bf = (Fn*(1-Fe)*(1-Fp)*(Feqnu-Fnu))/(1-F_tilde_nu)
return Bf
} 

//weak magnetism
double Wm(double e_nu)
{
return (1+1.1*e_nu)/n_rm
}

double Wm_bar(double e_nu_bar)
{
return (1-7.1*e_nu_bar)/n_rm
}

//Opacities

//electron neutrino absorption on neutrons
double nu_n_abs(double rho, double temp, double ye, double e_nu)
{
double sigma_nun_abs = sigma_o*((1+3*pow(gA,2))/4)*((e_nu+delta_np)/e_rm)*pow((1-pow((e_rm/(e_nu+delta_np)),2)),0.5)*Wm(e_nu)
return sigma_nun_abs 
}

//electron neutrino absorption on protonss
double nu_p_abs(double rho, double temp, double ye, double e_nu_bar)
{
double sigma_nup_abs = sigma_o*((1+3*pow(gA,2))/4)*((e_nu_bar-delta_np)/e_rm)*pow((1-pow((e_rm/(e_nu_bar-delta_np)),2)),0.5)*Wm_bar(e_nu_bar)
return sigma_nup_abs
}

int main (){
double rho =;
double temp =;
double ye=;
nun = nu_n_abs(rho, temp, ye);
nup = nu_p_abs(rho, temp, ye);
std::cout << "Neutrino absorption on neutrons for" << "rho =" << rho << "," << "temp =" << temp << "is" << nup <<std:endl;
std::cout << "Neutrino absorption on protonss for" << "rho =" << rho << "," << "temp =" << temp << "is" << nun <<std:endl;
}
