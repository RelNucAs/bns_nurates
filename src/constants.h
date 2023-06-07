#ifndef BNS_NURATES_SRC_CONSTANTS_H_
#define BNS_NURATES_SRC_CONSTANTS_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file constants.h
//  \brief constants used throughout the code

// Units conversion
const double kMeV = 1.602176634e-6;       // [erg]
const double kCm = 10e13;                 // [fm]

// Electron mass
const double kMe = 0.510998928;           // [MeV]

// Muon mass
const double kMmu = 105.6583745;          // [MeV]

// Neutron mass
const double kMn = 939.5654133;           // [MeV]

// Proton mass
const double kMp = 938.2720813;           // [MeV]
const double kQ = 1.2935;                 // [MeV]

// Atomic mass unit
const double kMb = 1.674e-24;             // [g]
const double kMu = 1.66054e-24;           // [g]

// Physical constants
const double kHbar = 6.582119569e-22;     // [MeV s]
const double kH = 4.1356943e-21;          // [MeV s]
const double kClight = 2.997924562e+10;   // [cm s^-1]
const double kGf = 8.957e-44;             // [MeV cm^3]
const double kPi = 3.1415926535898;       //  acos(static_cast<double>(-1.));
const double kKB = 8.617333262e-11;       // [MeV K^-1]

// Coupling constants
const double kGa = 1.23;
const double kGv = 1.;
const double kGs = 0.;
const double kSinsqthetaw = 0.2325;

// Neutral current nucleon form factors (kQ^2=0)
const double kHnv = -0.5;
const double kHna = -0.5 * kGa;
const double kHpv = 0.5 - 2. * kSinsqthetaw;
const double kHpa = 0.5 * kGa;

// Constants for the pair reaction
const double kGSqr = 1.55 * 1E-33;                  // kCm^3 kMeV^-2 s^-1, from Bruenn Eqn. (C51)
const double kAlpha1E = 1. + 2. * kSinsqthetaw;     // Table 1, Pons et. al. (1998)
const double kAlpha2E = 2. * kSinsqthetaw;          // Table 1, Pons et. al. (1998)
const double kAlpha1X = -1. + 2. * kSinsqthetaw;    // Table 1, Pons et. al. (1998)
const double kAlpha2X = 2. * kSinsqthetaw;          // Table 1, Pons et. al. (1998)
const double kAlpha1[2] = {kAlpha1E, kAlpha1X};
const double kAlpha2[2] = {kAlpha2E, kAlpha2X};

// Constants for the bremsstrahlung reaction
const double kBremXmin = 1.e-10;
const double kBremYmin = 1.e-10;
const double kBremEtamin = 1.e-10;
#endif //BNS_NURATES_SRC_CONSTANTS_H_
