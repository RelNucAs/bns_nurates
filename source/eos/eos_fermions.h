#pragma once

#include "../constants_bis.h"
#include "complete_FG.h"
#include "find_eta.h"
#include "interp.h"

using namespace constants;

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

	//printf("F3/2 = %.20e\n", f1);
	//printf("F5/2 = %.20e\n", f2);
	//printf("test = %.20e\n", K/3.*mLep*MeV*pow(theta,2.5)*(2.*f1+theta*f2));
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
	
	//printf("%.12e\t", fabs(-eta*f1+(5./3.-eta*theta)*f2 + 4./3.*theta*f3)/std::max(std::max(fabs(-eta*f1),fabs((5./3.-eta*theta)*f2)),fabs(4./3.*theta*f3)));
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



std::array<double,9> eos_ferm_onthefly(double eta, double T, int species) {
        std::array<double,9> eos_array;
        double mLep = 0;
        double f12, f32, f52, a_f12, a_f32, a_f52;

//.......Select particle mass
        if (species == 1) {
          mLep = me;
        } else if (species == 2) {
          mLep = mmu;
        }

//.......Define theta (relativity parameter)
        const double theta = T/mLep;

	f12 = compute_res(eta, theta, 0.5);
	f32 = compute_res(eta, theta, 1.5);
	f52 = compute_res(eta, theta, 2.5);

	const double a_eta = -(eta+2./theta);

	a_f12 = compute_res(a_eta, theta, 0.5);
	a_f32 = compute_res(a_eta, theta, 1.5);
	a_f52 = compute_res(a_eta, theta, 2.5);
	
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.);
	const double K3 = K/3.;

	//.......Compute eos_array[0]
	//eos_array[0] = n_part(T, eta, mLep);
	eos_array[0] = K*pow(theta,1.5)*(f12+theta*f32);

	//.......Compute eos_array[1]
	//eos_array[1] = n_antipart(T, eta, mLep);
	eos_array[1] = K*pow(theta,1.5)*(a_f12+theta*a_f32);

	//.......Compute eos_array[2]
	//eos_array[2] = P_part(T, eta, mLep);
	eos_array[2] = K3*mLep*MeV*pow(theta,2.5)*(2.*f32+theta*f52); //erg/cm^3

	//.......Compute eos_array[3]
	//eos_array[3] = P_antipart(T, eta, mLep)
	eos_array[3] = K3*mLep*MeV*pow(theta,2.5)*(2.*a_f32+theta*a_f52); //erg/cm^3

	//.......Compute eos_array[4]
	//eos_array[4] = e_part(T, eta, mLep);
	eos_array[4] = K*mLep*MeV*pow(theta,1.5)*(f12+2.*theta*f32+theta*theta*f52); //erg/cm^3

	//.......Compute eos_array[5]
	//eos_array[5] = e_antipart(T, eta, mLep);
	eos_array[5] = K*mLep*MeV*pow(theta,1.5)*(a_f12+2.*theta*a_f32+theta*theta*a_f52); //erg/cm^3

	//.......Compute eos_array[6]
	//eos_array[6] = s_part(T, eta, mLep);
	eos_array[6] = K*kB*MeV*pow(theta,1.5)*(-eta*f12+(5./3.-eta*theta)*f32+4./3.*theta*f52); //erg/K/cm^3

	//.......Compute eos_array[7]
	//eos_array[7] = s_antipart(T, eta, mLep);
	eos_array[7] = K*kB*MeV*pow(theta,1.5)*(-a_eta*a_f12+(5./3.-a_eta*theta)*a_f32+4./3.*theta*a_f52); //erg/K/cm^3

	//.......Compute eos_array[8]
	eos_array[8] = eta;

	return eos_array;
}

template <int species>
std::array<double,9> eos_ferm_fromNR(const double nLep, const double temp) {
	struct n_net_FD FDout;
	const double guess = find_guess_eta<species>(1.e39*nLep, temp);
        const double eta = rtsafe_mod<species>(1.e39*nLep, temp, guess, FDout);

	double mLep = 0;
        double f12, f32, f52, a_f12, a_f32, a_f52;

//.......Select particle mass
        if (species == 1) {
          mLep = me;
        } else if (species == 2) {
          mLep = mmu;
        }

//.......Define theta (relativity parameter)
        const double theta = temp/mLep;

	f12 = FDout.f12; //compute_res(eta, theta, 0.5);
	f32 = FDout.f32; //compute_res(eta, theta, 1.5);
	f52 = compute_res(eta, theta, 2.5);

	const double a_eta = -(eta+2./theta);

	a_f12 = FDout.a_f12; //compute_res(a_eta, theta, 0.5);
	a_f32 = FDout.a_f32; //compute_res(a_eta, theta, 1.5);
	a_f52 = compute_res(a_eta, theta, 2.5);
	
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.);
	const double K3 = K/3.;

        std::array<double,9> eos_array;
	
	//.......Compute eos_array[0]
	//eos_array[0] = n_part(T, eta, mLep);
	eos_array[0] = K*pow(theta,1.5)*(f12+theta*f32);

	//.......Compute eos_array[1]
	//eos_array[1] = n_antipart(T, eta, mLep);
	eos_array[1] = K*pow(theta,1.5)*(a_f12+theta*a_f32);

	//.......Compute eos_array[2]
	//eos_array[2] = P_part(T, eta, mLep);
	eos_array[2] = K3*mLep*MeV*pow(theta,2.5)*(2.*f32+theta*f52); //erg/cm^3

	//.......Compute eos_array[3]
	//eos_array[3] = P_antipart(T, eta, mLep)
	eos_array[3] = K3*mLep*MeV*pow(theta,2.5)*(2.*a_f32+theta*a_f52); //erg/cm^3

	//.......Compute eos_array[4]
	//eos_array[4] = e_part(T, eta, mLep);
	eos_array[4] = K*mLep*MeV*pow(theta,1.5)*(f12+2.*theta*f32+theta*theta*f52); //erg/cm^3

	//.......Compute eos_array[5]
	//eos_array[5] = e_antipart(T, eta, mLep);
	eos_array[5] = K*mLep*MeV*pow(theta,1.5)*(a_f12+2.*theta*a_f32+theta*theta*a_f52); //erg/cm^3

	//.......Compute eos_array[6]
	//eos_array[6] = s_part(T, eta, mLep);
	eos_array[6] = K*kB*MeV*pow(theta,1.5)*(-eta*f12+(5./3.-eta*theta)*f32+4./3.*theta*f52); //erg/K/cm^3

	//.......Compute eos_array[7]
	//eos_array[7] = s_antipart(T, eta, mLep);
	eos_array[7] = K*kB*MeV*pow(theta,1.5)*(-a_eta*a_f12+(5./3.-a_eta*theta)*a_f32+4./3.*theta*a_f52); //erg/K/cm^3

	//.......Compute eos_array[8]
	eos_array[8] = eta;

	return eos_array;
}


template<int eos_method, int species>
std::array<double,9> eos_ferm_array(double nLep, double t, struct EOSeta &eta_table, struct EOScomplete &EOS_table) {
//!........................INPUT........................................................
//!...... matter density rho [g/cm^3]; net particle fraction y_tilde; temperature T [MeV];
//!.......species 1 for electrons, 2 for muons; interval for non relat. degen. parameter [eta1,eta2]
//!.......el_interp = .true. to interpolate electrons, el_interp = .false. to use newrap 1D
//!.......mu_interp = .true. to interpolate muons, mu_interp = .false. to use newrap 1D
//!.....................................................................................
//!.......................OUTPUT........................................................
//!.......eos_array(1) = number density of particles [1/cm^3]...........................
//!.......eos_array(2) = number density of antiparticles [1/cm^3].......................
//!.......eos_array(3) = pressure of particles [erg/cm^3]...............................
//!.......eos_array(4) = pressure of antiparticles [erg/cm^3]...........................
//!.......eos_array(5) = energy of particles [erg/g]....................................
//!.......eos_array(6) = energy of antiparticles [erg/g]................................
//!.......eos_array(7) = entropy of particles [erg/K/g].................................
//!.......eos_array(8) = entropy of antiparticles [erg/K/g].............................
//!.......eos_array(9) = non relat. degeneracy parameter of particles...................
//!.......For a complete description of output variables, see subroutines above.........

	double mLep = 0;
	//double eta1_Lep, eta2_Lep;
	//bool guess_from_interp = false;

//.......Select particle mass
        if constexpr(species == 1) {
          mLep = me;
	  //eta1_Lep = eta1_e;
	  //eta2_Lep = eta2_e;
	} else if constexpr(species == 2) {
	  //eta1_Lep = eta1_mu;
	  //eta2_Lep = eta2_mu;
	  mLep = mmu;
	}

	if constexpr(eos_method == 3) {
		std::array<double,9> eos_array = eos_interp(nLep, t, EOS_table);
		eos_array[8] = (eos_array[8] - mLep) / t;
		return eos_array;
	} else if constexpr(eos_method == 1) {
		double guess = find_guess_eta<species>(1.e39*nLep, t);
		double eta = rtsafe<species>(1.e39*nLep, t, guess);
		return eos_ferm_onthefly(eta, t, species);
		//return eos_ferm_fromNR<species>(nLep, t);
	} else if constexpr(eos_method == 2) {
        	double eta = eos_tintep(nLep, t, eta_table.nL, eta_table.t, eta_table.eta);
		//if (guess_from_interp == true) {
		//	guess = eta;
		//	eta = rtsafe(1.e39*nLep, T, mLep, guess, eta1_Lep, eta2_Lep);i
		//}
		return eos_ferm_onthefly(eta, t, species);
	}
}
