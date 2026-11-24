//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nms_direct_interp.hpp
//  \brief contains kernels for inelastic neutrino scattering
//         off muons
//
// Computation of inelastic neutrino-muon scattering 
// via kernel grid interpolation.
//
// @TODO: ADD REFERENCE 


#ifndef BNS_NURATES_INCLUDE_KERNEL_NMAS_DIRECT_INTERP_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_NMAS_DIRECT_INTERP_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"
#include "scatt_off_muon_kernel_table.hpp"
#include "scatt_off_antimuon_kernel_table.hpp"


/* 4D linear interpolation of the numerical table of the three parts of the kernel
 * for neutrino scattering off (anti)muons
 *
 *      Inputs:
 *          T:            temperature                [MeV]
 *          mu:           muon chemical potential    [MeV]
 *          w :           incoming neutrino energy   [MeV]
 *          wp :          outgoing neutrino energy   [MeV]
 *          *table:       RTable-type struct      
 *
 *      Output:
 *          R1, R2, R3:   Interpolated values        []
 */

#if !defined(__CUDA_ARCH__) && !defined(__SYCL_DEVICE_ONLY__)  // NMAS direct-interp table data lives in CPU memory only

KOKKOS_INLINE_FUNCTION
NMAS_KernelResult NMAS_direct_interpolator(const BS_REAL T, const BS_REAL mu, const BS_REAL w,
                                          const BS_REAL wp)
{
    constexpr BS_REAL one = 1;

    NMAS_KernelResult res = {};

    int i0, i1, j0, j1, k0, k1, l0, l1;
    BS_REAL tx, ty, tz, tw;

    if (MuonLinear_find_bracketing_indices(T, NMAS_T_axis, NMAS_T_dims, &i0, &i1, &tx) < 0 or
        MuonLinear_find_bracketing_indices(mu, NMAS_mu_axis, NMAS_mu_dims, &j0, &j1, &ty) < 0 or
        MuonLinear_find_bracketing_indices(w, NMAS_w_axis, NMAS_w_dims, &k0, &k1, &tz) < 0 or
        MuonLinear_find_bracketing_indices(wp, NMAS_wp_axis, NMAS_wp_dims, &l0, &l1, &tw) < 0)
    {
        return res;
    }

// macro to access the 1D array
#define NMAS_IDX(i, j, k, l)                                                        \
    (((i * NMAS_mu_dims + j) * NMAS_w_dims + k) * NMAS_wp_dims + l)

    BS_REAL a0 = (one - tx), a1 = tx;
    BS_REAL b0 = (one - ty), b1 = ty;
    BS_REAL c0 = (one - tz), c1 = tz;
    BS_REAL d0 = (one - tw), d1 = tw;

    // extract the 16 vertexes of the hypercube for the three parts of the kernel
    // R1, muons
    BS_REAL R1_mu_c0000 = R1_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R1_mu_c0001 = R1_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R1_mu_c0010 = R1_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R1_mu_c0011 = R1_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R1_mu_c0100 = R1_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R1_mu_c0101 = R1_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R1_mu_c0110 = R1_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R1_mu_c0111 = R1_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R1_mu_c1000 = R1_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R1_mu_c1001 = R1_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R1_mu_c1010 = R1_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R1_mu_c1011 = R1_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R1_mu_c1100 = R1_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R1_mu_c1101 = R1_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R1_mu_c1110 = R1_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R1_mu_c1111 = R1_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R1_mu = R1_mu_c0000 * a0 * b0 * c0 * d0 + R1_mu_c0001 * a0 * b0 * c0 * d1 +
                    R1_mu_c0010 * a0 * b0 * c1 * d0 + R1_mu_c0011 * a0 * b0 * c1 * d1 +
                    R1_mu_c0100 * a0 * b1 * c0 * d0 + R1_mu_c0101 * a0 * b1 * c0 * d1 +
                    R1_mu_c0110 * a0 * b1 * c1 * d0 + R1_mu_c0111 * a0 * b1 * c1 * d1 +
                    R1_mu_c1000 * a1 * b0 * c0 * d0 + R1_mu_c1001 * a1 * b0 * c0 * d1 +
                    R1_mu_c1010 * a1 * b0 * c1 * d0 + R1_mu_c1011 * a1 * b0 * c1 * d1 +
                    R1_mu_c1100 * a1 * b1 * c0 * d0 + R1_mu_c1101 * a1 * b1 * c0 * d1 +
                    R1_mu_c1110 * a1 * b1 * c1 * d0 + R1_mu_c1111 * a1 * b1 * c1 * d1;

    // R2, muons
    BS_REAL R2_mu_c0000 = R2_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R2_mu_c0001 = R2_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R2_mu_c0010 = R2_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R2_mu_c0011 = R2_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R2_mu_c0100 = R2_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R2_mu_c0101 = R2_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R2_mu_c0110 = R2_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R2_mu_c0111 = R2_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R2_mu_c1000 = R2_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R2_mu_c1001 = R2_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R2_mu_c1010 = R2_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R2_mu_c1011 = R2_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R2_mu_c1100 = R2_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R2_mu_c1101 = R2_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R2_mu_c1110 = R2_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R2_mu_c1111 = R2_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R2_mu = R2_mu_c0000 * a0 * b0 * c0 * d0 + R2_mu_c0001 * a0 * b0 * c0 * d1 +
                    R2_mu_c0010 * a0 * b0 * c1 * d0 + R2_mu_c0011 * a0 * b0 * c1 * d1 +
                    R2_mu_c0100 * a0 * b1 * c0 * d0 + R2_mu_c0101 * a0 * b1 * c0 * d1 +
                    R2_mu_c0110 * a0 * b1 * c1 * d0 + R2_mu_c0111 * a0 * b1 * c1 * d1 +
                    R2_mu_c1000 * a1 * b0 * c0 * d0 + R2_mu_c1001 * a1 * b0 * c0 * d1 +
                    R2_mu_c1010 * a1 * b0 * c1 * d0 + R2_mu_c1011 * a1 * b0 * c1 * d1 +
                    R2_mu_c1100 * a1 * b1 * c0 * d0 + R2_mu_c1101 * a1 * b1 * c0 * d1 +
                    R2_mu_c1110 * a1 * b1 * c1 * d0 + R2_mu_c1111 * a1 * b1 * c1 * d1;

    // R3, muons
    BS_REAL R3_mu_c0000 = R3_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R3_mu_c0001 = R3_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R3_mu_c0010 = R3_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R3_mu_c0011 = R3_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R3_mu_c0100 = R3_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R3_mu_c0101 = R3_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R3_mu_c0110 = R3_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R3_mu_c0111 = R3_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R3_mu_c1000 = R3_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R3_mu_c1001 = R3_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R3_mu_c1010 = R3_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R3_mu_c1011 = R3_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R3_mu_c1100 = R3_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R3_mu_c1101 = R3_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R3_mu_c1110 = R3_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R3_mu_c1111 = R3_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R3_mu = R3_mu_c0000 * a0 * b0 * c0 * d0 + R3_mu_c0001 * a0 * b0 * c0 * d1 +
                    R3_mu_c0010 * a0 * b0 * c1 * d0 + R3_mu_c0011 * a0 * b0 * c1 * d1 +
                    R3_mu_c0100 * a0 * b1 * c0 * d0 + R3_mu_c0101 * a0 * b1 * c0 * d1 +
                    R3_mu_c0110 * a0 * b1 * c1 * d0 + R3_mu_c0111 * a0 * b1 * c1 * d1 +
                    R3_mu_c1000 * a1 * b0 * c0 * d0 + R3_mu_c1001 * a1 * b0 * c0 * d1 +
                    R3_mu_c1010 * a1 * b0 * c1 * d0 + R3_mu_c1011 * a1 * b0 * c1 * d1 +
                    R3_mu_c1100 * a1 * b1 * c0 * d0 + R3_mu_c1101 * a1 * b1 * c0 * d1 +
                    R3_mu_c1110 * a1 * b1 * c1 * d0 + R3_mu_c1111 * a1 * b1 * c1 * d1;

    // R1, anti-muons
    BS_REAL R1_amu_c0000 = R1_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R1_amu_c0001 = R1_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R1_amu_c0010 = R1_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R1_amu_c0011 = R1_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R1_amu_c0100 = R1_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R1_amu_c0101 = R1_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R1_amu_c0110 = R1_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R1_amu_c0111 = R1_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R1_amu_c1000 = R1_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R1_amu_c1001 = R1_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R1_amu_c1010 = R1_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R1_amu_c1011 = R1_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R1_amu_c1100 = R1_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R1_amu_c1101 = R1_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R1_amu_c1110 = R1_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R1_amu_c1111 = R1_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R1_amu = R1_amu_c0000 * a0 * b0 * c0 * d0 + R1_amu_c0001 * a0 * b0 * c0 * d1 +
                     R1_amu_c0010 * a0 * b0 * c1 * d0 + R1_amu_c0011 * a0 * b0 * c1 * d1 +
                     R1_amu_c0100 * a0 * b1 * c0 * d0 + R1_amu_c0101 * a0 * b1 * c0 * d1 +
                     R1_amu_c0110 * a0 * b1 * c1 * d0 + R1_amu_c0111 * a0 * b1 * c1 * d1 +
                     R1_amu_c1000 * a1 * b0 * c0 * d0 + R1_amu_c1001 * a1 * b0 * c0 * d1 +
                     R1_amu_c1010 * a1 * b0 * c1 * d0 + R1_amu_c1011 * a1 * b0 * c1 * d1 +
                     R1_amu_c1100 * a1 * b1 * c0 * d0 + R1_amu_c1101 * a1 * b1 * c0 * d1 +
                     R1_amu_c1110 * a1 * b1 * c1 * d0 + R1_amu_c1111 * a1 * b1 * c1 * d1;

    // R2, anti-muons
    BS_REAL R2_amu_c0000 = R2_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R2_amu_c0001 = R2_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R2_amu_c0010 = R2_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R2_amu_c0011 = R2_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R2_amu_c0100 = R2_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R2_amu_c0101 = R2_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R2_amu_c0110 = R2_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R2_amu_c0111 = R2_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R2_amu_c1000 = R2_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R2_amu_c1001 = R2_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R2_amu_c1010 = R2_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R2_amu_c1011 = R2_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R2_amu_c1100 = R2_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R2_amu_c1101 = R2_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R2_amu_c1110 = R2_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R2_amu_c1111 = R2_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R2_amu = R2_amu_c0000 * a0 * b0 * c0 * d0 + R2_amu_c0001 * a0 * b0 * c0 * d1 +
                     R2_amu_c0010 * a0 * b0 * c1 * d0 + R2_amu_c0011 * a0 * b0 * c1 * d1 +
                     R2_amu_c0100 * a0 * b1 * c0 * d0 + R2_amu_c0101 * a0 * b1 * c0 * d1 +
                     R2_amu_c0110 * a0 * b1 * c1 * d0 + R2_amu_c0111 * a0 * b1 * c1 * d1 +
                     R2_amu_c1000 * a1 * b0 * c0 * d0 + R2_amu_c1001 * a1 * b0 * c0 * d1 +
                     R2_amu_c1010 * a1 * b0 * c1 * d0 + R2_amu_c1011 * a1 * b0 * c1 * d1 +
                     R2_amu_c1100 * a1 * b1 * c0 * d0 + R2_amu_c1101 * a1 * b1 * c0 * d1 +
                     R2_amu_c1110 * a1 * b1 * c1 * d0 + R2_amu_c1111 * a1 * b1 * c1 * d1;

    // R3, muons
    BS_REAL R3_amu_c0000 = R3_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l0)];
    BS_REAL R3_amu_c0001 = R3_anti_muon_scattering[NMAS_IDX(i0, j0, k0, l1)];
    BS_REAL R3_amu_c0010 = R3_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l0)];
    BS_REAL R3_amu_c0011 = R3_anti_muon_scattering[NMAS_IDX(i0, j0, k1, l1)];
    BS_REAL R3_amu_c0100 = R3_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l0)];
    BS_REAL R3_amu_c0101 = R3_anti_muon_scattering[NMAS_IDX(i0, j1, k0, l1)];
    BS_REAL R3_amu_c0110 = R3_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l0)];
    BS_REAL R3_amu_c0111 = R3_anti_muon_scattering[NMAS_IDX(i0, j1, k1, l1)];
    BS_REAL R3_amu_c1000 = R3_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l0)];
    BS_REAL R3_amu_c1001 = R3_anti_muon_scattering[NMAS_IDX(i1, j0, k0, l1)];
    BS_REAL R3_amu_c1010 = R3_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l0)];
    BS_REAL R3_amu_c1011 = R3_anti_muon_scattering[NMAS_IDX(i1, j0, k1, l1)];
    BS_REAL R3_amu_c1100 = R3_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l0)];
    BS_REAL R3_amu_c1101 = R3_anti_muon_scattering[NMAS_IDX(i1, j1, k0, l1)];
    BS_REAL R3_amu_c1110 = R3_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l0)];
    BS_REAL R3_amu_c1111 = R3_anti_muon_scattering[NMAS_IDX(i1, j1, k1, l1)];

    BS_REAL R3_amu = R3_amu_c0000 * a0 * b0 * c0 * d0 + R3_amu_c0001 * a0 * b0 * c0 * d1 +
                     R3_amu_c0010 * a0 * b0 * c1 * d0 + R3_amu_c0011 * a0 * b0 * c1 * d1 +
                     R3_amu_c0100 * a0 * b1 * c0 * d0 + R3_amu_c0101 * a0 * b1 * c0 * d1 +
                     R3_amu_c0110 * a0 * b1 * c1 * d0 + R3_amu_c0111 * a0 * b1 * c1 * d1 +
                     R3_amu_c1000 * a1 * b0 * c0 * d0 + R3_amu_c1001 * a1 * b0 * c0 * d1 +
                     R3_amu_c1010 * a1 * b0 * c1 * d0 + R3_amu_c1011 * a1 * b0 * c1 * d1 +
                     R3_amu_c1100 * a1 * b1 * c0 * d0 + R3_amu_c1101 * a1 * b1 * c0 * d1 +
                     R3_amu_c1110 * a1 * b1 * c1 * d0 + R3_amu_c1111 * a1 * b1 * c1 * d1;

    res.R1_mu = R1_mu;
    res.R2_mu = R2_mu;
    res.R3_mu = R3_mu;
    res.R1_amu = R1_amu;
    res.R2_amu = R2_amu;
    res.R3_amu = R3_amu;

#undef NMAS_IDX

    return res;
}


//==============================================================//
// --- Inelastic Neutrino Scattering off (anti)muons Kernel --- //
// ---        via DIRECT KERNEL INTERPOLATION               --- //
//==============================================================//
/* Interpolate the MonteCarlo grid of the kernel via 4D linear interpolation
 *
 * NOTE: this function is just for points inside the ranges of the 4D table!
 *
 *      Inputs:
 *          T:          temperature                [MeV]
 *          mu:         muon chemical potential    [MeV]
 *          w:          incoming neutrino energy   [MeV]
 *          wp:         outgoing neutrino energy   [MeV]
 *
 *      Output:
 *           kernel:    Inelastic Neutrino Scattering off (anti)muons Kernel  [nm^3 s^-1]
 */


// Calculates and saves the neutrino muon scattering in and out kernel
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticNMASKernels_DirectInterp(InelasticScattKernelParams* kernel_params,
                                        MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega;
    const BS_REAL wp                   = kernel_params->omega_prime;
    const BS_REAL mu_mu                = eos_params->mu_mu;
    const BS_REAL det_bal_arg          = (wp - w) / T;
    const BS_REAL det_bal              = SafeExp(det_bal_arg);

    MyKernelOutput nmas_kernel;

    // Interpolation of the 4D table
    const NMAS_KernelResult R123_values = NMAS_direct_interpolator(T, mu_mu, w, wp);

    // Neutrino-flavor dependent kernel
    // nu_e
    nmas_kernel.abs[id_nue] = kBS_NMAS_Conv_Const * (
                             kBS_NMAS_C1_nue * (R123_values.R1_mu + R123_values.R2_amu) +
                             kBS_NMAS_C2_nue * (R123_values.R2_mu + R123_values.R1_amu) +
                             kBS_NMAS_C3_nue * (R123_values.R3_mu + R123_values.R3_amu));
    
    // anu_e
    nmas_kernel.abs[id_anue] = kBS_NMAS_Conv_Const * (
                              kBS_NMAS_C2_nue * (R123_values.R1_mu + R123_values.R2_amu) +
                              kBS_NMAS_C1_nue * (R123_values.R2_mu + R123_values.R1_amu) +
                              kBS_NMAS_C3_nue * (R123_values.R3_mu + R123_values.R3_amu));
    
    // num
    nmas_kernel.abs[id_num] = kBS_NMAS_Conv_Const * (
                             kBS_NMAS_C1_numu * (R123_values.R1_mu + R123_values.R2_amu) +
                             kBS_NMAS_C2_numu * (R123_values.R2_mu + R123_values.R1_amu) +
                             kBS_NMAS_C3_numu * (R123_values.R3_mu + R123_values.R3_amu));
    
    // anum
    nmas_kernel.abs[id_anum] = kBS_NMAS_Conv_Const * (
                              kBS_NMAS_C2_numu * (R123_values.R1_mu + R123_values.R2_amu) +
                              kBS_NMAS_C1_numu * (R123_values.R2_mu + R123_values.R1_amu) + 
                              kBS_NMAS_C3_numu * (R123_values.R3_mu + R123_values.R3_amu));

    // nu_t
    nmas_kernel.abs[id_nut] = nmas_kernel.abs[id_nue];
    
    // anu_t
    nmas_kernel.abs[id_anut] = nmas_kernel.abs[id_anue];

    // Absorption and Emission kernels [nm^3 s^-1]
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        nmas_kernel.em[idx]  = nmas_kernel.abs[idx] * det_bal;
    }

    return nmas_kernel;
}

#else  // __CUDA_ARCH__ || __SYCL_DEVICE_ONLY__ — NMS direct-interp is host-only; provide a stub so device code links

KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticNMASKernels_DirectInterp(InelasticScattKernelParams* /*kernel_params*/,
                                        MyEOSParams* /*eos_params*/)
{
    MyKernelOutput nmas_kernel = {0};
    return nmas_kernel;
}

#endif  // !__CUDA_ARCH__ && !__SYCL_DEVICE_ONLY__


#endif // BNS_NURATES_INCLUDE_KERNEL_NMAS_DIRECT_INTERP_HPP_
