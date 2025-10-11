//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_brem_BRT06.hpp
//  \brief contains bremsstrahlung kernels and associated helper functions
//
// Computation of nucleon-nucleon bremsstrahlung kernel using the approach of
// A. Burrows et al. Nuclear Physics A 777 (2006) 356-394

#ifndef BNS_NURATES_INCLUDE_KERNEL_BREM_BRT06_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_BREM_BRT06_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"


// Bremsstrahlung fitting formula described in
// A. Burrows et al. Nuclear Physics A 777 (2006) 356-394
// * The factor 2.0778 is different from the paper 1.04 to account
//   for the nuclear matrix element for one-pion exchange
//   (Adam Burrows, private comm)
KOKKOS_INLINE_FUNCTION
BS_REAL QBrem_BRT06(const BS_REAL nb, const BS_REAL T, const BS_REAL xn,
                    const BS_REAL xp)
{
    constexpr BS_REAL half               = 0.5;
    constexpr BS_REAL twentyeight_thirds = 28. / 3.;
    constexpr BS_REAL eleven_halves      = 5.5;
    constexpr BS_REAL mb                 = kBS_Mb;

    constexpr BS_REAL kBS_Brem_BRT06_Const = 2.0778e+02;
    const BS_REAL rho                      = nb * mb; // mass density [g nm-3]
    return kBS_Brem_BRT06_Const * half * kBS_MeV *
           (POW2(xn) + POW2(xp) + twentyeight_thirds * xn * xp) * POW2(rho) *
           pow(T, eleven_halves); // [MeV nm-3 s-1]
}

// Bremsstrahlung kernel from BRT06 Eq.(143) rewritten consistently
// to fit within the framework of the present library
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsBRT06(BremKernelParams* kernel_params,
                                MyEOSParams* eos_pars)
{
    constexpr BS_REAL half = 0.5;

    const BS_REAL omega       = kernel_params->omega;
    const BS_REAL omega_prime = kernel_params->omega_prime;
    const BS_REAL temp        = eos_pars->temp;

    const BS_REAL x = half * (omega + omega_prime) / temp;
    const BS_REAL q_nb =
        QBrem_BRT06(eos_pars->nb, temp, eos_pars->yn, eos_pars->yp);

    const BS_REAL tmp = kBS_HClight6FourPiSquared * kBS_Brem_C4BRT06 *
                        (q_nb / POW7(temp)) * bessk1(x) / x;
    const BS_REAL s_em  = tmp * SafeExp(-x);
    const BS_REAL s_abs = tmp * SafeExp(x);

    MyKernelOutput brem_kernel;
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_BREM_BRT06_HPP_
