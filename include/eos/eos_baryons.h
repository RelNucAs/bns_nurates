#pragma once 

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>

#include "../parameters.h"

//Baryon EOS module to read in baryonic tables (Shen format) and to interpolate them
//Shen 98 format, which does not include electrons, positrons and photons
//WARNING: for DD2 choose nb_entries = 326, for BLh choose nb_entries = 303, for SHFo choose nb_entries = 326

using namespace std;

const int t_entries = 81; //number of temperature entries
const int y_entries = 60; //number of proton fraction entries
const int nb_entries = 326; //303 number of nb entries
const int entries = 19; //number of EOS entries per grid-point

//i=0, t_entries: temperature index
//j=0, y_entries: proton fraction index
//k=0, nb_entries: baryon number density index
//l=0, entries: EOS entry index

typedef std::vector<std::vector<std::vector<std::vector<double>>>> eos_type;
        
//std::array<double,t_entries> t; //temperature array
//std::array<double,y_entries> y; //nominal proton fraction array
//std::array<double,nb_entries> nb; //nominal baryon number density array [fm^-3]


//......define the constants of the t, y, and nb array
const double nbin = 1.0e-12;
const double tin = 0.1e0;
const double yin = 0.0e0;
const double dnb = 0.04e0;
const double dt = 0.04e0;
const double dy = 0.01e0;

std::array<double,t_entries> read_t_baryon() {
//.....Uncomment the following part, if you do not want to use the stored
//.....values of t (round-off error) which are read-in above, but the correctly calculated ones:
//.....calculate the temperature array
	std::array<double,t_entries> t; //temperature array
	for (int i=0;i<t_entries;i++) t[i] = tin * pow(10.0, (double) (i-1)*dt);
	return t;
}

std::array<double,y_entries> read_y_baryon() {
//.....calculate the electron fraction array
	std::array<double,y_entries> y; //nominal proton fraction array
	for (int j=0;j<y_entries;j++) y[j] = yin + (double) j * dy;
	return y;
}

std::array<double,nb_entries> read_nb_baryon() {
	std::array<double,nb_entries> nb; //nominal baryon number density array [fm^-3]
//.....calculate the baryon number density array
	for (int k=0;k<nb_entries;k++) nb[k] = nbin * pow(10.0, (double) (k-1)*dnb);
	return nb;
}

//.....eos(i,j,k,l): EOS array

//.....brief list of the entries:
//.....==========================
//.....1	logarithm of baryon mass density: log(rho) [g/cm^3]
//.....2	baryon number density: nb [fm^-3]
//.....3	logarithm of total proton fraction: log(yp)
//.....4	total proton fraction: yp 
//.....5	free energy per baryon wrt 938 MeV: f/nb-938 [MeV]
//.....6	internal energy per baryon wrt amu=931.49432 MeV: e/nb-amu [MeV]
//.....7	entropy per baryon: s/nb [kB]
//.....8	averaged mass number of heavy nuclei: <A>
//.....9	averaged charge number of heavy nuclei: <Z>
//.....10	effective mass: m^* [MeV]
//.....11	unbound neutron mass fraction: X_n
//.....12	unbound proton mass fraction: X_p
//.....13	light nuclei mass fraction: X_s
//.....14	heavy nuclei mass fraction: X_A
//.....15	pressure p: [MeV/fm^3]
//.....16	neutron chemical potential wrt m_n: mu_n-m_n [MeV]
//.....17	proton chemical potential wrt m_n: mu_p-m_p [MeV]
//.....18	averaged mass number of light nuclei: <a>
//.....19	averaged charge number of light nuclei: <z>



//This function reads in the ASCII EOS table and calculates the nominal proton fraction and nominal baryon number density arrays
//Input: eos_filename.tab
//std::array<std::array<std::array<std::array<double,entries>,nb_entries>,y_entries>,t_entries> eos_readin(std::string eos_filename) {
eos_type eos_readin() {
	double number;

	eos_type eos;
	std::vector<std::vector<std::vector<double>>> eos_nb_y;
	std::vector<std::vector<double>> eos_nb;
	std::vector<double> single_eos;
	int i, j, k, l;

        if (nb_entries==326) {
		printf("\n");
		printf("WARNING from baryon_eos_module\n");
		printf("nb_entries compatible with DD2 EOS or SFHo EOS\n");
		printf("\n");
	}  

        if (nb_entries==303) {
		printf("\n");
		printf("WARNING from baryon_eos_module\n");
		printf("nb_entries compatible with BLh EOS\n");
		printf("\n");
	}

        // open input table
	std::ifstream EOStable;
	std::string eos_filename = parameters::abs_path + "eos_table/dd2_frdm_eos_shen98format_v1.02.tab";
	EOStable.open(eos_filename);  // open
        
	std::string EOSline;

        if (!EOStable) {
                cout << "EOS table not found!" << endl;
                std::exit(EXIT_FAILURE);
        }

//.....t loop
        for (i=0;i<t_entries;i++) {
          getline(EOStable, EOSline); //read(10,*)
          getline(EOStable, EOSline); //read(10,21) dummy,t(i)
	  getline(EOStable, EOSline); //read(10,*)
//.....ye loop
          eos_nb_y.clear();
          for (j=0;j<y_entries;j++) {
//.....nb loop
            eos_nb.clear();
            for(k=0;k<nb_entries;k++) {
		getline(EOStable, EOSline);
		std::stringstream s(EOSline);
		single_eos.clear();
		for (l=0;l<entries;l++) {
			s >> number;
			single_eos.push_back(number);
	    	}
		eos_nb.push_back(single_eos);
	    }
	    eos_nb_y.push_back(eos_nb);
	    getline(EOStable, EOSline); //read(10,*)
	  }
	  eos.push_back(eos_nb_y);
	}


	EOStable.close();
        
	return eos;
}


//.....This function interpolates the table in the temperature
std::array<double,entries> eos_tinterp_bar(double inb, double iy, double it, std::array<double,t_entries> t, std::array<double,y_entries> y, std::array<double,nb_entries> nb, eos_type eos) { //status)
//-----------------------------------------------------------------------
//
//     input:
//     inb ... baryon number density [fm^-3]
//     iy ... proton fraction []
//     it ... temperature [MeV]
//
//     output:
//     ieos... array which stores the interpolated eos output
//             check the above list of the entries
//     status ... >0 on failure
//
//-----------------------------------------------------------------------

       int jnb, jy, jt;
       std::array<double,3> xmin, dx;
       std::array<double,entries> ieos;
       double lognb, rnb, ry, rt, logt;

	
//.....find number density index.........................................
//.....rnb is the weight (to be used in the interpolation)
//.....jnb is the table index of the maximum nb lower than inb
      lognb = log10(inb);
      xmin[0] = log10(nb[0]);
      dx[0] = 0.04;  //increment as in manual doc
      jnb = floor( (lognb - xmin[0])/dx[0] );
      if (jnb < 0) {
	     jnb = 0;
	     rnb = 0.;
      } else if (jnb >= nb_entries-1) {
	     jnb = nb_entries-2;
	     rnb = 1.;
      } else {
	     rnb = ( lognb - log10(nb[jnb]) )/dx[0];
      }

//.....find y index.....................................................
//.....ry is the weight (to be used in the interpolation)
//.....jy is the table index of the maximum y lower than iy
      xmin[1] = y[0];
      dx[1] = 0.01; //increment as in the manual
      jy = floor( (iy - xmin[1])/dx[1] );
      if (jy < 0) {
	     jy = 0;
	     ry = 0.;
	} else if (jy >= y_entries-1) {
	       jy = y_entries-2;
	       ry = 1.;
	} else {
		ry = (iy - y[jy])/dx[1];
	}

//.....find temperature index.........................................
//.....rt is the weight (to be used in the interpolation)
//.....jt is the table index of the maximum t lower than it
      logt = log10(it);
      xmin[2] = log10(t[0]);
      dx[2] = 0.04; //increment as in the manual
      jt = floor( (logt - xmin[2])/dx[2] );
      if (jt < 0) {
	      jt = 0;
	      rt = 0.;
	} else if (jt >= t_entries-1) {
		jt = t_entries-2;
		rt = 1.;
	} else {
		rt = ( logt - log10(t[jt]) )/dx[2];
	}

//.....interpolate and assign results................................
      for (int l=0;l<entries;l++) {
	      ieos[l] = (1.-rnb)*( (1.-ry)*( (1.-rt)*eos[jt][jy][jnb][l]
       +                         rt *eos[jt+1][jy][jnb][l] )
       +               ry *( (1.-rt)*eos[jt][jy+1][jnb][l]
       +                         rt *eos[jt+1][jy+1][jnb][l] ) )
       +     rnb *( (1.-ry)*( (1.-rt)*eos[jt][jy][jnb+1][l]
       +                         rt *eos[jt+1][jy][jnb+1][l] )
       +               ry *( (1.-rt)*eos[jt][jy+1][jnb+1][l]
       +                         rt *eos[jt+1][jy+1][jnb+1][l] ) );
      }

      return ieos;
}
