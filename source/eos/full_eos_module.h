#pragma once
//!.......Module to build the complete EOS at finite T by assembling the following contributions:
//!.......baryons, electrons, positrons, muons, antimuons, photons and neutrinos
//!.......Baryon contribution: tables in Shen98 format (see Matthias Hempel eos), while baryon_eos_module reads tables and interpolates them	
//!.......Electron and Positron contributions: defined in eos_fermions.f90
//!.......Muon and Antimuon contributions: defined in eos_fermions.f90
//!.......Photon contributions: defined in eos_radiation.f90
//!.......Neutrinos contribution: defined in eos_neutrinos.f90


//!.......Warning 1: fermionic eos does not include coulomb correction; but the coulomb correction
//!.......should be included in the baryon pressure (see Manual from Matthias Hempel)

//!.......Warning 2: in this code the mass density is the baryon number density times a.m.u.
//!.......The physical parameter is the baryon number density. The mass density is the mass 
//!.......density of baryons. If muons are not created, then the mass density is the physical one

//!.......The module contains two subroutines: eos_extended and eos_full. Notice that
//!.......eos_extended DOES NOT contain neutrinos and constraints, while eos_full
//!.......contains neutrinos and constrains. For a complete description see each subroutine.

#include "interp.h"
#include "eos_baryons.h" 
#include "eos_fermions.h"
#include "eos_radiation.h"
#include "eos_neutrinos.h"
#include "../constants_bis.h"
#include "../parameters.h"

//.......Define parameters..................................................................
	
bool e_interp = true;    //True: activate electrons table interpolation; false: newt-rap to find eta
bool mu_interp = true;   //True: activate muons table interpolation; false: newt-rap to find eta

//const double c = 1.23984193e3; //h*c in [MeV*fm]

std::array<double,26> eos_extended(double den, double temp, double ye, double ymu, double yp) { 

//!.......Input...........................................................................
//!.......density [g/cm**3] which is the baryon number density times the a.m.u. (range [1.d3,1.d16] g/cm^3)
//!.......temperature [MeV] (range [0.1,110] MeV)
//!.......proton fraction yp (range [0.01,0.6])
//!.......electronic net fraction ye[] (range (0,0.6])
//!.......muonic net fraction ymu[] (range (-0.01,0.2))

//!.......Notice that the baryonic table has to be called in the main program

//!.......Output: array eos_out(26)...............................................................
//!.......1: protons non relativistic degeneracy parameter [](proton mass subtracted)
//!.......2: neutrons non relativistic degeneracy parameter [](neutron mass subtracted)
//!.......3: electrons non relativistic degeneracy parameter []
//!.......4: muons non relativistic degeneracy parameter []
//!.......5: baryons pressure [MeV/fm**3]
//!.......6: electrons pressure [MeV/fm**3]
//!.......7: positrons pressure [MeV/fm**3]
//!.......8: muons pressure [MeV/fm**3]
//!.......9: anti-muons pressure [MeV/fm**3]
//!......10: radiation pressure [MeV/fm**3]
//!......11: total pressure [MeV/fm**3]
//!......12: baryons energy [MeV/baryon] (rest mass amu subtracted)
//!......13: electrons energy(rest mass not included) [MeV/baryon]
//!......14: positrons energy(rest mass included 2 times) [MeV/baryon]
//!......15: muons energy(rest mass not included) [MeV/baryon]
//!......16: antimuons energy(rest mass included two times) [MeV/baryon]
//!......17: radiation energy [MeV/baryon]
//!......18: total energy [MeV/baryon]
//!......19: baryons entropy [kb/baryon]
//!......20: electrons entropy [kb/baryon]
//!......21: positrons entropy [kb/baryon]
//!......22: muons entropy [kb/baryon]
//!......23: antimuons entropy [kb/baryon]
//!......24: photons entropy [kb/baryon]
//!......25: total entropy per baryon [kb/baryon]
//!......26: fraction of electrons []
	
	double n_bar; //, t_kelv
	double pres, ener, entr;
	double guess_e, guess_mu;
	std::array<double,entries> b_eos; //Eos of baryons, entries defined in baryon_eos_module
	std::array<double,9>       e_eos, m_eos; //Eos of electrons and muons, see eos_ferm_array in eos_fermions
	std::array<double,26>      eos_out;
	int s_el, s_mu; //to select electrons or muons
	int i;
	
//.......Define baryon number density [fm^-3]
	n_bar = den/(mu*1.e39);
//.......Define electronic species
	s_el = 1;
//.......Define muonic species
	s_mu = 2;

        eos_type dd2_eos;
        dd2_eos = eos_readin();

	std::array<double,t_entries>  t_array  = read_t_baryon(); //temperature array
        std::array<double,y_entries>  y_array  = read_y_baryon(); //nominal proton fraction array
        std::array<double,nb_entries> nb_array = read_nb_baryon(); //nominal baryon number density array [fm^-3]

//.......Call baryon eos (interpolated in temperature)
        b_eos = eos_tinterp_bar(n_bar,yp,temp,t_array,y_array,nb_array,dd2_eos);

//.......Call electronic eos
//.......Provide guess
//	guess_e = find_guess_e(1.e39*n_bar*ye,temp);
//      eta = rtsafe(1.e39*n_bar*ye, temp, me, guess_e, eta1_e, eta2_e);
	//eta1 = -1.3d4  !minimum non relat degeneracy parameter
	//eta2 = 1.3d4 !maximum non relat degeneracy parameter
	//guess_e = (eta1+eta2)*0.5
	//if(ye < 0.d0) then
	//  !e_interp_temp = .false.
	//  e_interp_temp = e_interp
	//else
	//  e_interp_temp = e_interp	  
	//end if
        
	e_eos = eos_ferm_array(den, ye, temp, s_el, 2);

	//e_interp_temp = e_interp

//.......Call muonic eos
	if (fabs(ymu) > 1.e-10) {
	// eta3 = -2.d4 !-5.d1  !minimum non relat degeneracy parameter
	// eta4 = 5.d4 !4.d3 !maximum non relat degeneracy parameter
	// guess_m = (eta3+eta4)*0.5
	// if(ymu .le. 0.01) then
        //   nq = (2.*pi*m_m*temp/(c**2.))**1.5
	//   guess_m = -log(abs(2.*nq/(ymu*n_bar)))
	// end if
		m_eos = eos_ferm_array(den, ymu, temp, s_mu, 2);
	} else {
		for (i=0;i<8;i++) {
			m_eos[i] = 0.;
		}
		m_eos[8] = -mmu/temp;
	}

//.......Protons non relat. degeneracy parameter	
	eos_out[0] = b_eos[16]/temp;
//.......Neutrons non relat. degen. parameter
	eos_out[1] = b_eos[15]/temp;
//.......Electrons non relat. degen. parameter in the interval [eta1,eta2]
	eos_out[2] = e_eos[8];	
//.......Muons non relat. degen. parameter in the interval [eta3,eta4]
	eos_out[3] = m_eos[8];
//.......Baryons pressure [MeV/fm**3]
	eos_out[4] = b_eos[14];
//.......Electrons pressure [MeV/fm**3]
	pres = e_eos[2]; 
	eos_out[5] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Positrons pressure [MeV/fm**3]
	pres = e_eos[3];
	eos_out[6] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Muons pressure [MeV/fm**3]
	pres = m_eos[2];
	eos_out[7] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Anti-muons pressure [MeV/fm**3]
	pres = m_eos[3];
	eos_out[8] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Radiation pressure [MeV/fm**3]
	pres = P_rad(temp);
	eos_out[9] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Total pressure [MeV/fm**3]
	pres = 0.;
	for (i=4;i<10;i++) pres = pres + eos_out[i];
	eos_out[10] = pres;
//.......Baryons energy [MeV/baryon]
	eos_out[11] = b_eos[5];
//.......Electrons energy [MeV/baryon]
	ener = e_eos[4];
	eos_out[12] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Positrons energy [MeV/baryon]
	ener = e_eos[5];
	eos_out[13] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Muons energy [MeV/baryon]
	ener = m_eos[4];
	eos_out[14] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Anti-muons energy [MeV/baryon]
	ener = m_eos[5];
	eos_out[15] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Radiation energy [MeV/baryon]
	ener = e_rad(temp);
	eos_out[16] = (ener*MeV/1.e39)/n_bar;     //change units from erg/cm**3 to MeV/baryon
//.......Total energy [Mev/baryon]
	ener = 0.;
	for (i=11;i<17;i++) ener = ener+eos_out[i];
	eos_out[17] = ener;
//.......Baryons entropy [kb/baryon]
	eos_out[18] = b_eos[6];
//.......Electrons entropy [kb/baryon]
	entr = e_eos[6];
	eos_out[19] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Positrons entropy [kb/baryon]
	entr = e_eos[7];
	eos_out[20] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Muons entropy [kb/baryon]
	entr = m_eos[6];
	eos_out[21] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Antimuons entropy [kb/baryon]
	entr = m_eos[7];
	eos_out[22] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//......Radiation entropy [kb/baryon]
	entr = s_rad(temp);
	eos_out[23] = (entr/1.e39)*MeV/(n_bar*kB); //change units from erg/K/cm**3 to kb/baryon
//.......Total entropy [kb/baryon]
	entr = 0.;
	for (i=18;i<24;i++) entr = entr+eos_out[i];
	eos_out[24] = entr;
//.......Fraction of electrons
	eos_out[25] = e_eos[0]*mu/den;
	
	return eos_out;
}


std::array<double,53> eos_full(double den, double temp, double ye, double ymu) {
//.......Input...........................................................................
//.......density [g/cm**3] which is the baryon number density times the a.m.u. (range [1.d3,1.d16] g/cm^3)
//.......temperature [MeV] (range [0.1,110] MeV)
//.......electronic net fraction ye[] (range (0,0.6])
//.......muonic net fraction ymu[] (range (-0.01,0.2])


//!.......Output: array eos_out(53)...............................................................
//!.......1: protons non relativistic degeneracy parameter [](proton mass subtracted)
//!.......2: neutrons non relativistic degeneracy parameter [](neutron mass subtracted)
//!.......3: electrons non relativistic degeneracy parameter []
//!.......4: muons non relativistic degeneracy parameter []
//!.......5: electronic neutrinos degeneracy parameter[]
//!.......6: muonic neutrinos degeneracy parameter[]
//!.......7: baryons pressure [MeV/fm**3]
//!.......8: electrons pressure [MeV/fm**3]
//!.......9: positrons pressure [MeV/fm**3]
//!......10: muons pressure [MeV/fm**3]
//!......11: anti-muons pressure [MeV/fm**3]
//!......12: radiation pressure [MeV/fm**3]
//!......13: electr. neutrinos pressure [MeV/fm**3]
//!......14: electr. anti-neutrinos pressure [MeV/fm**3]
//!......15: muonic neutrinos pressure [MeV/fm**3]
//!......16: muonic anti-neutrinos pressure [MeV/fm**3]
//!......17: tau neutrinos+anti-neutrinos pressure [MeV/fm**3]
//!......18: total pressure [MeV/fm**3]
//!......19: baryons energy [MeV/baryon] (rest mass amu subtracted)
//!......20: electrons energy(rest mass not included) [MeV/baryon]
//!......21: positrons energy(rest mass included 2 times) [MeV/baryon]
//!......22: muons energy(rest mass not included) [MeV/baryon]
//!......23: antimuons energy(rest mass included two times) [MeV/baryon]
//!......24: radiation energy [MeV/baryon]
//!......25: elect. neutrinos energy [MeV/baryon]
//!......26: elect. anti-neutrinos energy [MeV/baryon]
//!......27: muonic neutrinos energy [MeV/baryon]
//!......28: muonic anti-neutrinos energy [MeV/baryon]
//!......29: tau neutrinos+anti-neutrinos energy [MeV/baryon]
//!......30: total energy [MeV/baryon]
//!......31: fraction of electrons []
//!......32: fraction of positrons []
//!......33: fraction of muons []
//!......34: fraction of anti-muons []
//!......35: electr. neutrinos fraction []
//!......36: electr. anti-neutrinos fraction []
//!......37: muonic neutrinos fraction []
//!......38: muonic anti-neutrinos fraction[]
//!......39: tau neutrinos fraction []
//!......40: tau anti-neutrinos fraction []
//!......41: fraction of protons []
//!......42: baryons entropy [kb/baryon]
//!......43: electrons entropy [kb/baryon]
//!......44: positrons entropy [kb/baryon]
//!......45: muons entropy [kb/baryon]
//!......46: antimuons entropy [kb/baryon]
//!......47: photons entropy [kb/baryon]
//!......48: elect. neutrinos entropy [kb/baryon]
//!......49: elect. anti-neutrinos entropy [kb/baryon]
//!......50: muonic neutrinos entropy [kb/baryon]
//!......51: muonic anti-neutrinos entropy [kb/baryon]
//!......52: tau neutrinos+antineutrinos entropy [kb/baryon]
//!......53: total entropy per baryon [kb/baryon]
//
//!.......Constraints: 1)protons fraction is the sum of electrons and muons fraction, that is yp=ye+ymu
//!.......2)electronic and muonic neutrinos coupled with electrons and muons via beta-equilibrium
//!.......3)only trapped neutrinos are considered (see limit densities on Perego,Bernuzzi,Radice paper)

//.......den is density in g/cm**3, temp is temperature in MeV
//.......ye is the net electronic fraction, ymu is the net muonic fraction

	double n_bar, yp; //, t_kelv
        double pres, ener, entr;
        double eta_ne, eta_nmu; //Degeneracy parameters of electronic and muonic neutrinos
	double yne, yane, ynmu, yanmu, yntau, yantau; //Fractions of neutrinos and antineutrinos (3 flavors)
	double ener_ne, ener_ane, ener_nmu, ener_anmu; //Energies of electronic and muonic neutrinos and anti-neu 
	double den_lime, den_limx;
        std::array<double,entries> b_eos; //Eos of baryons, entries defined in baryon_eos_module
        std::array<double,9>       e_eos, m_eos; //Eos of electrons and muons, see eos_ferm_array in eos_fermions
        std::array<double,53>      eos_out;
        int s_el, s_mu; //to select electrons or muons
        int i;
        bool e_interp_temp, mu_interp_temp;

//.......Define baryon number density [fm^-3]
        n_bar = den/(mu*1.e39);
//.......Define electronic species
        s_el = 1;
//.......Define muonic species
        s_mu = 2;

//.......Define temperature in K (kb is boltz. constant in [MeV/K] defined in eos_fermions)
	//t_kelv = temp/kb 
//.......Define the proton fraction
	yp = ye + ymu;
//.......Define limit density [g/cm**3] of electronic and muonic neutrinos
	den_lime = 1.e11;
//.......Define limit density [g/cm**3] of tauonic neutrinos
	den_limx = 1.e12;	

	eos_type dd2_eos;
        dd2_eos = eos_readin();

	std::array<double,t_entries>  t_array  = read_t_baryon(); //temperature array
	std::array<double,y_entries>  y_array  = read_y_baryon(); //nominal proton fraction array
	std::array<double,nb_entries> nb_array = read_nb_baryon(); //nominal baryon number density array [fm^-3]

//.......Call baryon eos (interpolated in temperature)
        b_eos = eos_tinterp_bar(n_bar,yp,temp,t_array,y_array,nb_array,dd2_eos);

//.......Call electronic eos
//.......Provide guess
	//guess_e = find_guess_e(1.e39*n_bar*ye,temp)
	//eta1 = -1.3d4  !minimum non relat degeneracy parameter
	//eta2 = 1.3d4 !maximum non relat degeneracy parameter
	//guess_e = (eta1+eta2)*0.5

	//if(ye < 0.d0) then
	  //e_interp_temp = .false.
	//else
	  //e_interp_temp = e_interp  
	//end if

        e_eos = eos_ferm_array(den, ye, temp, s_el, 2);

        //e_interp_temp = e_interp

	
//.......Call muonic eos
//	if (abs(ymu) < 1.d-10) then
	if (fabs(ymu) <= 1.e-17) {
		for (i=0;i<8;i++) m_eos[i] = 0.;
	  	m_eos[8] = -mmu/temp;
	} else {
//	  //eta3 = -500.  !minimum non relat degeneracy parameter (funziona per t>5 mev)
//	  //eta4 = 1.d3 !maximum non relat degeneracy parameter   (funziona per t>5 mev)
	  //eta3 = -1.d4  !minimum non relat degeneracy parameter
	  //eta4 = 1.3d4 !maximum non relat degeneracy parameter
	  //guess_m = -100.d0
//	  //guess_m = (eta3+eta4)*0.5
//	  //if(temp .ge. 20.) then
//            //nq = (2.*pi*m_m*temp/(c**2.))**1.5
//	    //guess_m = -log(abs(2.*nq/(ymu*n_bar)))	    
//	  //end if
	//if (abs(ymu) < 1.d-10 .and. abs(ymu) > 1.d-17) then
	  //mu_interp_temp = .false.
	//else
      //mu_interp_temp = mu_interp 
	//end if
		m_eos = eos_ferm_array(den, ymu, temp, s_mu, 2);
	//mu_interp_temp = mu_interp
	//end if
        }

	
//.......Protons non relat. degeneracy parameter	
	eos_out[0] = b_eos[16]/temp;
//.......Neutrons non relat. degen. parameter
	eos_out[1] = b_eos[15]/temp;
//.......Electrons non relat. degen. parameter in the interval [eta1,eta2]
	eos_out[2] = e_eos[8];
//.......Muons non relat. degen. parameter in the interval [eta3,eta4]
	eos_out[3] = m_eos[8];
//.......Electronic neutrinos degeneracy parameter (beta equilibrium)
//.......m_e is electron mass defined in eos_fermions
	eta_ne = eos_out[0]-eos_out[1]+eos_out[2]+(mp-mn+me)/temp;
	eos_out[4] = eta_ne;
//.......Muonic neutrinos degeneracy parameter (beta equilibrium)
//.......m_m is the muonic mass defined in eos_fermions
	if(fabs(ymu) < 1.e-10) {
         eta_nmu = 0.;
	} else {
	 eta_nmu = eos_out[0]-eos_out[1]+eos_out[3]+(mp-mn+mmu)/temp;
	}
	eos_out[5] = eta_nmu;
//.......Baryons pressure [MeV/fm**3]
	eos_out[6] = b_eos[14];
//.......Electrons pressure [MeV/fm**3]
	pres = e_eos[2];
	eos_out[7] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Positrons pressure [MeV/fm**3]
	pres = e_eos[3];
	eos_out[8] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Muons pressure [MeV/fm**3]
	pres = m_eos[2];
	eos_out[9] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Anti-muons pressure [MeV/fm**3]
	pres = m_eos[3];
	eos_out[10] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Radiation pressure [MeV/fm**3]
	pres = P_rad(temp);
	eos_out[11] = (pres/1.e39)*MeV; //convert pres from erg/cm^3 to MeV/fm^3
//.......Electr. neutrinos pressure [MeV/fm**3]
	ener_ne  = nu_energy(den,  temp, eta_ne);
	ener_ane = anu_energy(den, temp, eta_ne);
	eos_out[12] = n_bar*ener_ne*exp(-den_lime/den)/3.;
//.......Electr. anti-neutrinos pressure [MeV/fm**3]
	eos_out[13] = n_bar*ener_ane*exp(-den_lime/den)/3.;
//.......Muonic neutrinos pressure [MeV/fm**3]
	ener_nmu  = nu_energy(den,  temp, eta_nmu);
	ener_anmu = anu_energy(den, temp, eta_nmu);
	eos_out[14] = n_bar*ener_nmu*exp(-den_limx/den)/3.;
//.......Muonic anti-neutrinos pressure [MeV/fm**3]
	eos_out[15] = n_bar*ener_anmu*exp(-den_limx/den)/3.;
//.......Tau neutrinos+antineutrinos pressure (use sum relations)
	eos_out[16] = (7./45.)*pow(pi,5.)*pow(temp,4.)*exp(-den_limx/den)/pow(c,3.);
//.......Total pressure [MeV/fm**3]
	pres = 0.;
	for (i=6;i<17;i++) pres = pres+eos_out[i];
	eos_out[17] = pres;
//.......Baryons energy [MeV/baryon]
	eos_out[18] = b_eos[5];
//.......Electrons energy [MeV/baryon]
	ener = e_eos[4];
	eos_out[19] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Positrons energy [MeV/baryon]
	ener = e_eos[5];
	eos_out[20] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Muons energy [MeV/baryon]
	ener = m_eos[4];
	eos_out[21] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Anti-muons energy [MeV/baryon]
	ener = m_eos[5];
	eos_out[22] = (den/1.e39)*MeV*ener/n_bar; //change units from erg/g to MeV/baryon
//.......Radiation energy [MeV/baryon]
	ener = e_rad(temp);
	eos_out[23] = (ener*MeV/1.e39)/n_bar; //change units from erg/cm**3 to MeV/baryon
//.......Electr. neutrinos energy [MeV/baryon]
	eos_out[24] = ener_ne*exp(-den_lime/den);
//.......Electr. anti-neutrinos energy [MeV/baryon]
	eos_out[25] = ener_ane*exp(-den_lime/den);
//.......Muonic neutrinos energy [MeV/baryon]
	eos_out[26] = ener_nmu*exp(-den_limx/den);
//.......Muonic anti-neutrinos energy [MeV/baryon]
	eos_out[27] = ener_anmu*exp(-den_limx/den);
//.......Tau neutrinos+anti-neutrinos energy [MeV/baryon] (use sum relations)
	eos_out[28] = ((7./15.)*pow(pi,5.)*pow(temp,4.)/(n_bar*pow(c,3.)))*exp(-den_limx/den);
//.......Total energy [Mev/baryon]
	ener = 0.;
	for (i=18;i<29;i++) ener = ener+eos_out[i];
	eos_out[29] = ener;
//.......Fraction of electrons
	eos_out[30] = e_eos[0]*mu/den;
//.......Fraction of positrons
	eos_out[31] = e_eos[1]*mu/den;
//.......Fraction of muons
	eos_out[32] = m_eos[0]*mu/den;
//.......Fraction of anti-muons
	eos_out[33] = m_eos[1]*mu/den;
//.......Fraction of electronic neutrinos
	yne  = nu_fraction(den,  temp, eta_ne);
	yane = anu_fraction(den, temp, eta_ne);
	eos_out[34] = yne*exp(-den_lime/den);
//.......Fraction of electronic anti-neutrinos
	eos_out[35] = yane*exp(-den_lime/den);
//.......Fraction of muonic neutrinos
	ynmu  = nu_fraction(den,  temp, eta_nmu);
	yanmu = anu_fraction(den, temp, eta_nmu);
	eos_out[36] = ynmu*exp(-den_limx/den);
//.......Fraction of muonic anti-neutrinos	
	eos_out[37] = yanmu*exp(-den_limx/den);
//.......Fraction of tauonic neutrinos
	yntau  = nu_fraction(den,  temp, (double) 0.);
	yantau = anu_fraction(den, temp, (double) 0.);
	eos_out[38] = yntau*exp(-den_limx/den);
//.......Fraction of tauonic anti-neutrinos
	eos_out[39] = yantau*exp(-den_limx/den);
//.......Fraction of protons
	eos_out[40] = yp;
//.......Baryons entropy [kb/baryon]
	eos_out[41] = b_eos[6];
//.......Electrons entropy [kb/baryon]
	entr = e_eos[6];
	eos_out[42] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Positrons entropy [kb/baryon]
	entr = e_eos[7];
	eos_out[43] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Muons entropy [kb/baryon]
	entr = m_eos[6];
	eos_out[44] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Antimuons entropy [kb/baryon]
	entr = m_eos[7];
	eos_out[45] = (den/(1.e39*kB))*(MeV/n_bar)*entr; //change units from erg/K/g to kb/baryon
//.......Radiation entropy [kb/baryon]
	entr = s_rad(temp);
	eos_out[46] = (entr/1.e39)*MeV/(n_bar*kB); //change units from erg/K/cm**3 to kb/baryon 
//.......Electr. neutrinos entropy [kb/baryon]
	eos_out[47] = (eos_out[24]/temp)+(eos_out[12]/temp/n_bar)-(eta_ne*eos_out[34]);
//.......Electr. anti-neutrinos entropy [kb/baryon]
	eos_out[48] = (eos_out[25]/temp)+(eos_out[13]/temp/n_bar)-(-eta_ne*eos_out[35]);
//.......Muonic neutrinos entropy [kb/baryon]
	eos_out[49] = (eos_out[26]/temp)+(eos_out[14]/temp/n_bar)-(eta_nmu*eos_out[36]);
//.......Muonic anti-neutrinos entropy [kb/baryon]
	eos_out[50] = (eos_out[27]/temp)+(eos_out[15]/temp/n_bar)-(-eta_nmu*eos_out[37]);
//.......Tau neutrinos + anti-neutr. entropy [kb/baryon]
	eos_out[51] = (eos_out[28]/temp)+(eos_out[16]/temp/n_bar);
//.......Total entropy [kb/baryon]
	entr = 0.;
	for (i=41;i<52;i++) entr = entr+eos_out[i];

	eos_out[52] = entr;

	return eos_out;
}
