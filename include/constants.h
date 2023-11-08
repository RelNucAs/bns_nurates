#ifndef BNS_NURATES_SRC_CONSTANTS_H_
#define BNS_NURATES_SRC_CONSTANTS_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file constants.h
//  \brief constants used throughout the code

// @TODO: increase precision in definition of constants when possible

// Units conversion
static const double kMeV = 1.602176634e-6;       // [erg]
static const double kCm = 10e13;                 // [fm]

// Electron mass
static const double kMe = 0.510998928;           // [MeV]

// Muon mass
static const double kMmu = 105.6583745;          // [MeV]

// Neutron mass
static const double kMn = 939.5654133;           // [MeV]

// Proton mass
static const double kMp = 938.2720813;           // [MeV]
static const double kQ = 1.2935;                 // [MeV]

// Atomic mass unit
static const double kMb = 1.674e-24;             // [g]
static const double kMu = 1.66054e-24;           // [g]

// Physical constants
static const double kHbar = 6.582119569e-22;       // [MeV s]
static const double kH = 4.1356943e-21;            // [MeV s]
static const double kClight = 2.997924562e+10;     // [cm s^-1]
static const double kHClight = kH * kClight;       // [MeV cm]
static const double kHbarClight = kHbar * kClight; // [MeV cm]
static const double kGf0 = 1.1663787E-11;          // [MeV^-2]
// @TODO: decide how to define the Fermi constant
static const double kGf = kGf0 * kHbarClight * kHbarClight * kHbarClight;               // [MeV cm^3]
//static const double kGf = 8.957e-44;               // [MeV cm^3]
static const double kPi = 3.1415926535898;         //  acos(static_cast<double>(-1.));
static const double kKB = 8.617333262e-11;         // [MeV K^-1]

// Coupling constants
static const double kGa = 1.23;
static const double kGv = 1.;
static const double kGs = 0.;
static const double kSinsqthetaw = 0.2325;

// Neutral current nucleon form factors (kQ^2=0)
static const double kHnv = -0.5;
static const double kHna = -0.5 * kGa;
static const double kHpv = 0.5 - 2. * kSinsqthetaw;
static const double kHpa = 0.5 * kGa;

// Constants for the pair reaction
static const double kGSqr = 1.55 * 1E-33;                  // kCm^3 kMeV^-2 s^-1, from Bruenn Eqn. (C51)
static const double kAlpha1E = 1. + 2. * kSinsqthetaw;     // Table 1, Pons et. al. (1998)
static const double kAlpha2E = 2. * kSinsqthetaw;          // Table 1, Pons et. al. (1998)
static const double kAlpha1X = -1. + 2. * kSinsqthetaw;    // Table 1, Pons et. al. (1998)
static const double kAlpha2X = 2. * kSinsqthetaw;          // Table 1, Pons et. al. (1998)
static const double kAlpha1[2] = {kAlpha1E, kAlpha1X};
static const double kAlpha2[2] = {kAlpha2E, kAlpha2X};

#endif //BNS_NURATES_SRC_CONSTANTS_H_
