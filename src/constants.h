#ifndef BNS_NURATES_SRC_CONSTANTS_H_
#define BNS_NURATES_SRC_CONSTANTS_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file constants.h
//  \brief constants used for the rest of the code

// Units conversion
const double MeV = 1.602176634e-6;    //convert MeV to CGS (erg)
const double cm = 10e13;              //convert fm to cm

// Electron mass
const double me = 0.510998928;        // MeV

// Muon mass
const double mmu = 105.6583745;       // MeV

// Neutron mass
const double mn = 939.5654133;        // MeV

// Proton mass
const double mp = 938.2720813;        // MeV
const double Q = 1.2935;              // MeV

// Atomic mass unit
const double mb = 1.674e-24;          // g
const double mu = 1.66054e-24;        //g

// Physical constants
const double h_bar = 6.582 * (10e-22);      // MeV*s
const double h = 4.1356943e-21;             // MeV*s
const double c = 2.997924562e+10;           // cm/s
const double GF = 8.957e-44;                // MeV*cm^3
const double pi = 3.1415926535898;          //acos(static_cast<double>(-1.));
const double kB = 8.617333262e-11;          // MeV/K

// Coupling constants
const double gA = 1.23;
const double gV = 1.;
const double gS = 0.;
const double sinsqthetaw = 0.2325;

// Neutral current nucleon form factors (Q^2=0)
const double hnv = -0.5;
const double hna = -0.5 * gA;
const double hpv = 0.5 - 2. * sinsqthetaw;
const double hpa = 0.5 * gA;

// Constants for the pair reaction
const double G_sqr = 1.55 * 1E-33;                  // [cm^3 MeV^-2 s^-1], from Bruenn Eqn. (C51)
const double alpha_1_e = 1. + 2. * sinsqthetaw;     // Table 1, Pons et. al. (1998)
const double alpha_2_e = 2. * sinsqthetaw;          // Table 1, Pons et. al. (1998)
const double alpha_1_x = -1. + 2. * sinsqthetaw;    // Table 1, Pons et. al. (1998)
const double alpha_2_x = 2. * sinsqthetaw;          // Table 1, Pons et. al. (1998)
const double alpha_1[2] = {alpha_1_e, alpha_1_x};
const double alpha_2[2] = {alpha_2_e, alpha_2_x};

#endif //BNS_NURATES_SRC_CONSTANTS_H_
