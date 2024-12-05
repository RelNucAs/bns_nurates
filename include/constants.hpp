#ifndef BNS_NURATES_SRC_CONSTANTS_H_
#define BNS_NURATES_SRC_CONSTANTS_H_

#include "bns_nurates.hpp"

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file constants.h
//  \brief constants used throughout the code

////////////////////////////////////////////////////////////////////////////////
// All values are consistent with the CODATA Internationally recommended 2022 //
// values of the Fundamental Physical Constants                               //
// (https://pml.nist.gov/cuu/Constants/index.html)                            //
////////////////////////////////////////////////////////////////////////////////

//////////////////////
// UNIT CONVERSIONS //
//////////////////////

// MeV to nm^2 g s^-2
inline constexpr double kBS_MeV_double = 1.6021766341182763e+8;
inline constexpr BS_REAL kBS_MeV       = kBS_MeV_double;


////////////////////////
// PHYSICAL CONSTANTS //
////////////////////////

// Electron mass
inline constexpr double kBS_Me = 0.51099895069; // [MeV]

// Muon mass
inline constexpr double kBS_Mmu = 105.6583755; // [MeV]

// Neutron mass
inline constexpr double kBS_Mn = 939.56542194; // [MeV]

// Proton mass
inline constexpr double kBS_Mp = 938.27208943; // [MeV]

// Neutron-proton mass difference
inline constexpr double kBS_Q = 1.29333251; // [MeV]

// Atomic mass unit
inline constexpr double kBS_Mu = 1.66053906892e-24; // [g]

// Average baryon mass
inline constexpr double kBS_Mb = 1.6737747126723517e-24; // [g]

// pi
inline constexpr double kBS_Pi = 3.14159265358979323846264338328;

// Planck's constant
inline constexpr double kBS_H = 4.135667696e-21; // [MeV s]

// Reduced Planck's constant
inline constexpr double kBS_Hbar = 0.5 * kBS_H / kBS_Pi; // [MeV s]

// Speed of light
inline constexpr double kBS_Clight = 2.99792458e+17; // [nm s^-1]

// Boltzmann's constant
inline constexpr double kBS_KB = 8.617333262e-11; // [MeV K^-1]

// Reduced Fermi constant and Fermi constant
inline constexpr double kBS_Gf0 = 1.1663787e-11; // [MeV^-2]
inline constexpr double kBS_Gf =
    kBS_Gf0 * POW3(kBS_Hbar * kBS_Clight); // [MeV nm^3]

// Nuclear saturation number density
inline constexpr BS_REAL kBS_Saturation_n = 0.15e+18; // [nm^-3]


////////////////////////////////////////
// COUPLING CONSTANTS AND FORM FACTORS//
////////////////////////////////////////

inline constexpr double kBS_Ga      = 1.23;
inline constexpr double kBS_Gv      = 1.;
inline constexpr double kBS_Gs      = 0.;
inline constexpr double kBS_SinThW2 = 0.2325;

// Neutral current nucleon form factors (kBS_Q^2=0)
inline constexpr double kBS_Hnv = -0.5;
inline constexpr double kBS_Hna = -0.5 * kBS_Ga;
inline constexpr double kBS_Hpv = 0.5 - 2. * kBS_SinThW2;
inline constexpr double kBS_Hpa = 0.5 * kBS_Ga;


//////////////////////////////
// BETA-PROCESSES CONSTANTS //
//////////////////////////////

// Constant in Eq.(C13,C15,C19,C20) of Bruenn
inline constexpr BS_REAL kBS_Beta_Const =
    (POW3(kBS_Clight) * POW2(kBS_Gf0 * kBS_Hbar)) *
    (POW2(kBS_Gv) + 3. * POW2(kBS_Ga)) / kBS_Pi;


///////////////////////////////////////////
// NEUTRINO-NUCLEON SCATTERING CONSTANTS //
///////////////////////////////////////////

// Constant in Fermi energy calculation
// The constant 4.78... is 0.5 * (3 * pi^2)^(2./3.)
inline constexpr BS_REAL kBS_Iso_eF =
    4.785390000313652 * POW2(kBS_Hbar) / kBS_Mb * kBS_MeV_double;

inline constexpr double kBS_Iso_Const =
    2. * kBS_Pi * POW2(kBS_Gf) / kBS_Hbar / POW3(kBS_H * kBS_Clight);

inline constexpr double kBS_Iso_h0_p = POW2(kBS_Hpv) + 3. * POW2(kBS_Hpa);
inline constexpr double kBS_Iso_h1_p = POW2(kBS_Hpv) - POW2(kBS_Hpa);
inline constexpr double kBS_Iso_h0_n = POW2(kBS_Hnv) + 3. * POW2(kBS_Hna);
inline constexpr double kBS_Iso_h1_n = POW2(kBS_Hnv) - POW2(kBS_Hna);

// 0th Legendre coefficient (protons)
inline constexpr BS_REAL kBS_Iso_c0_p = kBS_Iso_Const * kBS_Iso_h0_p;
// 1st Legendre coefficient (protons)
inline constexpr BS_REAL kBS_Iso_c1_p = kBS_Iso_Const * kBS_Iso_h1_p;
// 0th Legendre coefficient (neutrons)
inline constexpr BS_REAL kBS_Iso_c0_n = kBS_Iso_Const * kBS_Iso_h0_n;
// 1st Legendre coefficient (neutrons)
inline constexpr BS_REAL kBS_Iso_c1_n = kBS_Iso_Const * kBS_Iso_h1_n;


//////////////////////////////////////////////
// NUCLEON-NUCLEON BREMSSTRAHLUNG CONSTANTS //
//////////////////////////////////////////////

inline constexpr BS_REAL kBS_Brem_Xmin   = 1.0e-10;
inline constexpr BS_REAL kBS_Brem_Ymin   = 1.0e-10;
inline constexpr BS_REAL kBS_Brem_Etamin = 1.0e-10;

inline constexpr BS_REAL kBS_Brem_Aux1 =
    0.25 * POW2(kBS_H) * kBS_MeV_double / POW2(kBS_Pi);

// -1.26 / 2
inline constexpr BS_REAL kBS_Brem_Ca = -0.63;

inline constexpr BS_REAL kBS_Brem_Const =
    POW2(kBS_Brem_Ca) * POW2(kBS_Gf) / kBS_Hbar;

// (3*5*7*11)/2^{11} / 4 (C in BRT06 Eq.143 divided by four)
inline constexpr BS_REAL kBS_Brem_C4BRT06 = 0.1409912109375;


////////////////////////////
// PAIR PROCESS CONSTANTS //
////////////////////////////

// Constant in Bruenn Eq. (C51)
inline constexpr BS_REAL kGSqr = 1.55e-12; // [nm^3 MeV^-2 s^-1]

inline constexpr BS_REAL kBS_Pair_Phi = kGSqr / kBS_Pi;

// Table 1, Pons et. al. (1998)
inline constexpr BS_REAL kBS_Pair_Alpha1E = 1. + 2. * kBS_SinThW2;
// Table 1, Pons et. al. (1998)
inline constexpr BS_REAL kBS_Pair_Alpha2E = 2. * kBS_SinThW2;
// Table 1, Pons et. al. (1998)
inline constexpr BS_REAL kBS_Pair_Alpha1X = -1. + 2. * kBS_SinThW2;
// Table 1, Pons et. al. (1998)
inline constexpr BS_REAL kBS_Pair_Alpha2X = 2. * kBS_SinThW2;

inline constexpr BS_REAL kBS_Pair_Alpha1_0 = kBS_Pair_Alpha1E;
inline constexpr BS_REAL kBS_Pair_Alpha1_1 = kBS_Pair_Alpha1X;
inline constexpr BS_REAL kBS_Pair_Alpha2_0 = kBS_Pair_Alpha2E;
inline constexpr BS_REAL kBS_Pair_Alpha2_1 = kBS_Pair_Alpha2X;


/////////////////////////////////////////////////////
// NEUTRINO-ELECTRON/POSITRON SCATTERING CONSTANTS //
/////////////////////////////////////////////////////

inline constexpr BS_REAL kBS_NEPS_Const = 2. * POW2(kBS_Gf0) * kBS_Clight *
                                          POW2(kBS_Hbar * kBS_Clight) /
                                          (3. * kBS_Pi);

inline constexpr BS_REAL kBS_NEPS_BPlus  = POW2(2. * kBS_SinThW2 + 1.);
inline constexpr BS_REAL kBS_NEPS_BMinus = POW2(2. * kBS_SinThW2 - 1.);
inline constexpr BS_REAL kBS_NEPS_BZero  = 4. * POW2(kBS_SinThW2);


/////////////////////////////
// MISCELLANEOUS CONSTANTS //
/////////////////////////////

// Cut-off values for safe exponential evaluation
inline constexpr BS_REAL kBS_ExpUppLim = std::is_same_v<BS_REAL, float>  ? 80 :
                                         std::is_same_v<BS_REAL, double> ? 700 :
                                                                           80;
inline constexpr BS_REAL kBS_ExpLowLim = -kBS_ExpUppLim;

inline constexpr BS_REAL kBS_HClight6FourPiSquared =
    POW6(kBS_H * kBS_Clight) / (16. * POW2(kBS_Pi));

// Energy scale for weak magnetism
inline constexpr BS_REAL kBS_WM_e_scale =
    kBS_MeV_double / (kBS_Mb * POW2(kBS_Clight));

static const BS_REAL kBS_FourPi_hc3 =
    (4. * kBS_Pi) / POW3(kBS_H * kBS_Clight); // [MeV^-3 nm^-3]
static const BS_REAL kBS_FourPi_hc3_sqr =
    POW2(kBS_FourPi_hc3); // [MeV^-6 nm^-6]

// Neutron mass in grams
inline constexpr BS_REAL kBS_MnGrams =
    kBS_Mn * kBS_MeV_double / POW2(kBS_Clight);

// Proton mass in grams
inline constexpr BS_REAL kBS_MpGrams =
    kBS_Mp * kBS_MeV_double / POW2(kBS_Clight);

// Geometric mean of nucleon masses in grams
inline constexpr BS_REAL kBS_MAvgGrams =
    kBS_MeV_double / POW2(kBS_Clight) * 938.918532994116257276238085018813;

// sqrt(pi)
inline constexpr BS_REAL kBS_SqrtPi = 1.77245385090551602729816748334106;
// kBS_Pi * pi
inline constexpr BS_REAL kBS_PiSquared = 9.86960440108935861883449099987572;
// 4 * piSquared
inline constexpr BS_REAL kBS_FourPiSquared = 39.4784176043574344753379639995029;
// pi^(1/8)
inline constexpr BS_REAL kBS_Pi2OneEighth = 1.15383506784998943054096521314981;
// (pi / 2)^2.5
inline constexpr BS_REAL kBS_PiHalfToFiveHalves = 3.09242868139914350627854469835251;


#endif // BNS_NURATES_SRC_CONSTANTS_H_
