//
// Created by maitraya on 6/2/23.
//

#ifndef BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
#define BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_

#include "constants.hpp"

#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))

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
const double lamp = 1.793;  //  proton magnetic moment?
const double lamn = -1.913; // neutron magnetic moment

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
void NucFrmFac(const double E, double* cv, double* ca, double* F2,
               const int reacflag)
{
    /* (Anti)neutrino energy rescaled by the nucleon mass */
    const double ehor =
        E * kMeV / (kMb * kClight * kClight); // Eq.(4), dimensionless

    const double tau = 0.5 * ehor * ehor / (1. + ehor);           // Eq.(B10)
    const double eta = 1. / (1. + 5.6 * tau);                     // Eq.(B16)
    const double G   = 1. / pow(1. + 4.97 * tau, 2.);             // Eq.(B17)
    const double Fp1 = (1. + tau * (1. + lamp)) * G / (1. + tau); // Eq.(B11)
    const double Fp2 = lamp * G / (1. + tau);                     // Eq.(B12)
    const double Fn1 = tau * lamn * (1. - eta) * G / (1. + tau);  // Eq.(B13)
    const double Fn2 = lamn * (1. + tau * eta) * G / (1. + tau);  // Eq.(B14)

    double frm1, frm2, frm3;

    /* Different parametrization depending on the reaction */
    if (reacflag == 1)
    {
        frm1 = (0.5 - 2. * kSinsqthetaw) * Fp1 - 0.5 * Fn1;  // Eq.(B1)
        frm2 = 0.5 * (kGa - kGs) / pow(1. + 3.53 * tau, 2.); // Eq.(B2)
        frm3 = (0.5 - 2. * kSinsqthetaw) * Fp2 - 0.5 * Fn2;  // Eq.(B3)
    }
    else if (reacflag == 2)
    {
        frm1 = (0.5 - 2. * kSinsqthetaw) * Fn1 - 0.5 * Fp1;   // Eq.(B4)
        frm2 = -0.5 * (kGa + kGs) / pow(1. + 3.53 * tau, 2.); // Eq.(B5)
        frm3 = (0.5 - 2. * kSinsqthetaw) * Fn2 - 0.5 * Fp2;   // Eq.(B6)
    }
    else if (reacflag == 3)
    {
        frm1 = Fp1 - Fn1;                      // Eq.(B7)
        frm2 = kGa / pow(1. + 3.53 * tau, 2.); // Eq.(B8)
        frm3 = Fp2 - Fn2;                      // Eq.(B9)
    }
    else
    {
        printf("Error: reacflag out of range in NucFrmFac\n");
        exit(EXIT_FAILURE);
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
void WMAbsEm(const double omega, double* R, double* Rbar)
{
    double cv, ca, F2;

    NucFrmFac(omega, &cv, &ca, &F2, 3); // nuclear form factors

    const double ehor = omega * kMeV / (kMb * kClight * kClight); // Eq.(4)
    const double tmp1 = cv * cv * (1. + 4. * ehor + 16. / 3. * ehor * ehor) +
                        3. * ca * ca * POW2(1. + 4. / 3. * ehor) +
                        8. / 3. * cv * F2 * ehor * ehor +
                        5. / 3. * ehor * ehor * (1. + 2. / 5. * ehor) * F2 * F2;
    const double tmp2 = 4. * (cv + F2) * ca * ehor * (1. + 4. / 3. * ehor);
    // const double tmp3 = (cv*cv+3.0*ca*ca)*POW3(1.+2.*ehor);
    const double tmp3 = (kGv * kGv + 3.0 * kGa * kGa) * POW3(1. + 2. * ehor);

    *R    = (tmp1 + tmp2) / tmp3; // Eq.(22)
    *Rbar = (tmp1 - tmp2) / tmp3; // Eq.(22)

    return;
}

// Correction for (anti)neutrino scattering on nucleons (nu + N -> nu + N):
// reacflag = 1 | 2 Input: omega -> neutrino energy [MeV] Output: correction to
// zeroth (R0) and first Legendre (R1) coefficients of scattering kernel
KOKKOS_INLINE_FUNCTION
void WMScatt(const double omega, double* R0, double* R1, const int reacflag)
{
    double cv, ca, F2;
    double h0, h1;

    NucFrmFac(omega, &cv, &ca, &F2, reacflag); // nuclear form factors
    // NucFrmFac(0., &cv_0, &ca_0, &F2_0, reacflag); //nuclear form factors at
    // Q^2=0

    // @TODO: evaluate this at compile time
    if (reacflag == 1)
    {
        h0 = kHpv * kHpv + 3. * kHpa * kHpa;
        h1 = kHpv * kHpv - kHpa * kHpa;
    }
    else if (reacflag == 2)
    {
        h0 = kHnv * kHnv + 3. * kHna * kHna;
        h1 = kHnv * kHnv - kHna * kHna;
    }

    const double ehor = omega * kMeV / (kMb * kClight * kClight);

    /* Low-energy limit derived from Eq.(12) */
    *R0 = (cv * cv + 3. * ca * ca + 1.5 * ehor * ehor * F2 * F2 +
           2. * ehor * ehor * cv * F2) /
          h0; // correction to zeroth coefficient
    *R1 = (cv * cv - ca * ca - 2. * ehor * ehor * F2 * F2 -
           4. * ehor * ehor * cv * F2) /
          h1; // correction to first coefficient

    return;
}

/*===========================================================================*/

#endif // BNS_NURATES_SRC_OPACITIES_WEAK_MAGNETISM_WEAK_MAGNETISM_H_
