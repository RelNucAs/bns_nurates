//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
#define BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_

#include "bns_nurates.hpp"
#include "constants.hpp"

/*===========================================================================*/

// nucfrmfac.c
//  \brief Calculation of single nucleon form factors as in C.J. Horowitz, 2002
//         (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001).
//         These are needed to compute the recoil and weak magnetism correction
//         for (anti)neutrino absorption on nucleons and elastic scattering on
//         nucleons.

// Reactions are distinguished using the following indices:
//   reacflag = 0: (anti)neutrino scattering on proton  (nu p -> nu p)
//   reacflag = 1: (anti)neutrino scattering on neutron (nu n -> nu n)
//   reacflag = 2: (anti)neutrino absorption on nucleon (nue n -> e- p, anue p
//   -> e+ n)

// Nucleon constants
constexpr BS_REAL lamp = 1.793;  //  proton magnetic moment?
constexpr BS_REAL lamn = -1.913; // neutron magnetic moment

// Computation of single nucleon form factors for reaction reacflag,
// given the (anti)neutrino energy

/*
 * Input:
 * 	- E : (anti)neutrino energy [MeV]
 * 	- reacflag : index defining the reaction (see above)
 *
 * Output:
 * 	- cv : vector form factor
 * 	- ca : axial vector form factor
 * 	- F2 : tensor/Pauli form factor
 *
 */

KOKKOS_INLINE_FUNCTION
void NucFrmFac(const BS_REAL E, BS_REAL* cv, BS_REAL* ca, BS_REAL* F2,
               const int reacflag)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL two  = 2;
    constexpr BS_REAL half = 0.5;

    constexpr BS_REAL c1 = 5.6;
    constexpr BS_REAL c2 = 4.97;
    constexpr BS_REAL c3 = 3.53;

    constexpr BS_REAL sinthw2 = kBS_SinThW2;
    constexpr BS_REAL ga      = kBS_Ga;
    constexpr BS_REAL gs      = kBS_Gs;

    // (Anti)neutrino energy rescaled by the nucleon mass, Eq. 4
    const BS_REAL ehor = E * kBS_WM_e_scale; // dimensionless

    const BS_REAL tau = half * POW2(ehor) / (one + ehor); // Eq.(B10)
    const BS_REAL eta = one / (one + c1 * tau);           // Eq.(B16)
    const BS_REAL G   = one / pow(one + c2 * tau, two);   // Eq.(B17)
    const BS_REAL Fp1 =
        (one + tau * (one + lamp)) * G / (one + tau);               // Eq.(B11)
    const BS_REAL Fp2 = lamp * G / (one + tau);                     // Eq.(B12)
    const BS_REAL Fn1 = tau * lamn * (one - eta) * G / (one + tau); // Eq.(B13)
    const BS_REAL Fn2 = lamn * (one + tau * eta) * G / (one + tau); // Eq.(B14)

    BS_REAL frm1, frm2, frm3;

    /* Different parametrization depending on the reaction */
    if (reacflag == 1)
    {
        frm1 = (half - two * sinthw2) * Fp1 - half * Fn1; // Eq.(B1)
        frm2 = half * (ga - gs) / POW2(one + c3 * tau);   // Eq.(B2)
        frm3 = (half - two * sinthw2) * Fp2 - half * Fn2; // Eq.(B3)
    }
    else if (reacflag == 2)
    {
        frm1 = (half - two * sinthw2) * Fn1 - half * Fp1; // Eq.(B4)
        frm2 = -half * (ga + gs) / POW2(one + c3 * tau);  // Eq.(B5)
        frm3 = (half - two * sinthw2) * Fn2 - half * Fp2; // Eq.(B6)
    }
    else if (reacflag == 3)
    {
        frm1 = Fp1 - Fn1;                 // Eq.(B7)
        frm2 = ga / POW2(one + c3 * tau); // Eq.(B8)
        frm3 = Fp2 - Fn2;                 // Eq.(B9)
    }
    else
    {
        BS_PRINTF("Error: reacflag out of range in NucFrmFac\n");
    }

    *cv = frm1;
    *ca = frm2;
    *F2 = frm3;

    return;
}


/*===========================================================================*/

// weak_magnetism.c

// !\file weak_magnetism.c
// \brief Evaluation of phase space/recoil/weak magnetism correction for
// (anti)neutrino
//        emission/absorption on nucleons and elastic scattering on nucleons
//        Ref: Horowitz, 2002
//        (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.65.043001)


// R    -> Correction for electron     neutrino absorption on neutron (nu_l + n
// -> l- + p) Rbar -> Correction for electron antineutrino absorption on proton
// (anu_l + p -> l+ + n) reacflag = 3 (for nuclear form factors) Input: omega ->
// neutrino energy [MeV]
KOKKOS_INLINE_FUNCTION
void WMAbsEm(const BS_REAL omega, BS_REAL* R, BS_REAL* Rbar)
{
    constexpr BS_REAL one   = 1;
    constexpr BS_REAL two   = 2;
    constexpr BS_REAL three = 3;
    constexpr BS_REAL four  = 4;

    constexpr BS_REAL four_thirds    = 4. / 3.;
    constexpr BS_REAL five_thirds    = 5. / 3.;
    constexpr BS_REAL eight_thirds   = 8. / 3.;
    constexpr BS_REAL sixteen_thirds = 16. / 3.;
    constexpr BS_REAL two_fifth      = 2. / 5.;

    constexpr BS_REAL ga = kBS_Ga;
    constexpr BS_REAL gv = kBS_Gv;

    BS_REAL cv, ca, F2;

    NucFrmFac(omega, &cv, &ca, &F2, 3); // nuclear form factors

    const BS_REAL ehor = omega * kBS_WM_e_scale;

    const BS_REAL tmp1 =
        POW2(cv) * (one + four * ehor + sixteen_thirds * POW2(ehor)) +
        three * ca * ca * POW2(one + four_thirds * ehor) +
        eight_thirds * cv * F2 * POW2(ehor) +
        five_thirds * POW2(ehor) * (one + two_fifth * ehor) * POW2(F2);
    const BS_REAL tmp2 =
        four * (cv + F2) * ca * ehor * (one + four_thirds * ehor);
    // const BS_REAL tmp3 = (cv*cv+3.0*ca*ca)*POW3(one+two*ehor);
    const BS_REAL tmp3 = (POW2(gv) + three * POW2(ga)) * POW3(one + two * ehor);

    *R    = (tmp1 + tmp2) / tmp3; // Eq.(22)
    *Rbar = (tmp1 - tmp2) / tmp3; // Eq.(22)

    return;
}

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N):
// reacflag = 1 | 2 Input: omega -> neutrino energy [MeV] Output: correction to
// zeroth (R0) and first Legendre (R1) coefficients of scattering kernel
KOKKOS_INLINE_FUNCTION
void WMScatt(const BS_REAL omega, BS_REAL* R0, BS_REAL* R1, const int reacflag)
{
    constexpr BS_REAL two          = 2;
    constexpr BS_REAL three        = 3;
    constexpr BS_REAL four         = 4;
    constexpr BS_REAL three_halves = 1.5;

    constexpr BS_REAL hpv = kBS_Hpv;
    constexpr BS_REAL hpa = kBS_Hpa;
    constexpr BS_REAL hnv = kBS_Hnv;
    constexpr BS_REAL hna = kBS_Hna;

    BS_REAL cv, ca, F2;
    BS_REAL h0, h1;

    NucFrmFac(omega, &cv, &ca, &F2, reacflag); // nuclear form factors
    // NucFrmFac(0., &cv_0, &ca_0, &F2_0, reacflag); //nuclear form factors at
    // Q^2=0

    // @TODO: evaluate this at compile time
    if (reacflag == 1)
    {
        h0 = POW2(hpv) + three * POW2(hpa);
        h1 = POW2(hpv) - POW2(hpa);
    }
    else if (reacflag == 2)
    {
        h0 = POW2(hnv) + three * POW2(hna);
        h1 = POW2(hnv) - POW2(hna);
    }

    const BS_REAL ehor = omega * kBS_WM_e_scale;

    /* Low-energy limit derived from Eq.(12) */
    // correction to zeroth coefficient
    *R0 = (POW2(cv) + three * POW2(ca) + three_halves * POW2(ehor * F2) +
           two * POW2(ehor) * cv * F2) /
          h0;
    // correction to first coefficient
    *R1 = (POW2(cv) - POW2(ca) - two * POW2(ehor * F2) -
           four * POW2(ehor) * cv * F2) /
          h1;

    return;
}

/*===========================================================================*/

#endif // BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
