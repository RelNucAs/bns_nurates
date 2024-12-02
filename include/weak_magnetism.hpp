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
const BS_REAL lamp = 1.793;  //  proton magnetic moment?
const BS_REAL lamn = -1.913; // neutron magnetic moment

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
    // (Anti)neutrino energy rescaled by the nucleon mass, Eq. 4
    const BS_REAL ehor = E * kBS_WM_e_scale; // dimensionless

    const BS_REAL tau = 0.5 * POW2(ehor) / (1. + ehor);            // Eq.(B10)
    const BS_REAL eta = 1. / (1. + 5.6 * tau);                     // Eq.(B16)
    const BS_REAL G   = 1. / pow(1. + 4.97 * tau, 2.);             // Eq.(B17)
    const BS_REAL Fp1 = (1. + tau * (1. + lamp)) * G / (1. + tau); // Eq.(B11)
    const BS_REAL Fp2 = lamp * G / (1. + tau);                     // Eq.(B12)
    const BS_REAL Fn1 = tau * lamn * (1. - eta) * G / (1. + tau);  // Eq.(B13)
    const BS_REAL Fn2 = lamn * (1. + tau * eta) * G / (1. + tau);  // Eq.(B14)

    BS_REAL frm1, frm2, frm3;

    /* Different parametrization depending on the reaction */
    if (reacflag == 1)
    {
        frm1 = (0.5 - 2. * kBS_SinThW2) * Fp1 - 0.5 * Fn1;      // Eq.(B1)
        frm2 = 0.5 * (kBS_Ga - kBS_Gs) / POW2(1. + 3.53 * tau); // Eq.(B2)
        frm3 = (0.5 - 2. * kBS_SinThW2) * Fp2 - 0.5 * Fn2;      // Eq.(B3)
    }
    else if (reacflag == 2)
    {
        frm1 = (0.5 - 2. * kBS_SinThW2) * Fn1 - 0.5 * Fp1;       // Eq.(B4)
        frm2 = -0.5 * (kBS_Ga + kBS_Gs) / POW2(1. + 3.53 * tau); // Eq.(B5)
        frm3 = (0.5 - 2. * kBS_SinThW2) * Fn2 - 0.5 * Fp2;       // Eq.(B6)
    }
    else if (reacflag == 3)
    {
        frm1 = Fp1 - Fn1;                      // Eq.(B7)
        frm2 = kBS_Ga / POW2(1. + 3.53 * tau); // Eq.(B8)
        frm3 = Fp2 - Fn2;                      // Eq.(B9)
    }
    else
    {
        printf("Error: reacflag out of range in NucFrmFac\n");
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
    BS_REAL cv, ca, F2;

    NucFrmFac(omega, &cv, &ca, &F2, 3); // nuclear form factors

    const BS_REAL ehor = omega * kBS_WM_e_scale;

    const BS_REAL tmp1 =
        POW2(cv) * (1. + 4. * ehor + 16. / 3. * POW2(ehor)) +
        3. * ca * ca * POW2(1. + 4. / 3. * ehor) +
        8. / 3. * cv * F2 * POW2(ehor) +
        5. / 3. * POW2(ehor) * (1. + 2. / 5. * ehor) * POW2(F2);
    const BS_REAL tmp2 = 4. * (cv + F2) * ca * ehor * (1. + 4. / 3. * ehor);
    // const BS_REAL tmp3 = (cv*cv+3.0*ca*ca)*POW3(1.+2.*ehor);
    const BS_REAL tmp3 =
        (POW2(kBS_Gv) + 3.0 * POW2(kBS_Ga)) * POW3(1. + 2. * ehor);

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
    BS_REAL cv, ca, F2;
    BS_REAL h0, h1;

    NucFrmFac(omega, &cv, &ca, &F2, reacflag); // nuclear form factors
    // NucFrmFac(0., &cv_0, &ca_0, &F2_0, reacflag); //nuclear form factors at
    // Q^2=0

    // @TODO: evaluate this at compile time
    if (reacflag == 1)
    {
        h0 = POW2(kBS_Hpv) + 3. * POW2(kBS_Hpa);
        h1 = POW2(kBS_Hpv) - POW2(kBS_Hpa);
    }
    else if (reacflag == 2)
    {
        h0 = POW2(kBS_Hnv) + 3. * POW2(kBS_Hna);
        h1 = POW2(kBS_Hnv) - POW2(kBS_Hna);
    }

    const BS_REAL ehor = omega * kBS_WM_e_scale;

    /* Low-energy limit derived from Eq.(12) */
    // correction to zeroth coefficient
    *R0 = (POW2(cv) + 3. * POW2(ca) + 1.5 * POW2(ehor * F2) +
           2. * POW2(ehor) * cv * F2) /
          h0;
    // correction to first coefficient
    *R1 = (POW2(cv) - POW2(ca) - 2. * POW2(ehor * F2) -
           4. * POW2(ehor) * cv * F2) /
          h1;

    return;
}

/*===========================================================================*/

#endif // BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
