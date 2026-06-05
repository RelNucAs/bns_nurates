//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nms_direct_interp.hpp
//  \brief contains kernels for inelastic neutrino scattering
//         off muons
//
// Computation of inelastic neutrino-muon scattering 
// via direct grid interpolation.
//
// @TODO: ADD REFERENCE 


#ifndef BNS_NURATES_INCLUDE_KERNEL_NMS_DIRECT_INTERP_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_NMS_DIRECT_INTERP_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"
#include "nms_kernel_table.hpp"


/* 4D linear interpolation of the numerical table of the three parts of the kernel
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

// Struct to store the three terms of the kernel
struct NMS_KernelResult{
    BS_REAL R1;
    BS_REAL R2;
    BS_REAL R3;
};

KOKKOS_INLINE_FUNCTION
NMS_KernelResult NMS_direct_interpolator(const BS_REAL T, const BS_REAL mu, const BS_REAL w,
                          const BS_REAL wp)
{
    constexpr BS_REAL one = 1;

    NMS_KernelResult res = {};

    int i0, i1, j0, j1, k0, k1, l0, l1;
    BS_REAL tx, ty, tz, tw;

    if (NMS_find_bracketing_indices(T, NMS_T_axis, NMS_T_dims, &i0, &i1, &tx) < 0 or
        NMS_find_bracketing_indices(mu, NMS_mu_axis, NMS_mu_dims, &j0, &j1, &ty) < 0 or
        NMS_find_bracketing_indices(w, NMS_w_axis, NMS_w_dims, &k0, &k1, &tz) < 0 or
        NMS_find_bracketing_indices(wp, NMS_wp_axis, NMS_wp_dims, &l0, &l1, &tw) < 0)
    {
        return res;
    }

// macro to access the 1D array
#define NMS_IDX(i, j, k, l)                                                        \
    (((i * NMS_mu_dims + j) * NMS_w_dims + k) * NMS_wp_dims + l)

    // extract the 16 vertexes of the hypercube for the three parts of the kernel
    // R1
    BS_REAL R1_c0000 = R1_data[NMS_IDX(i0, j0, k0, l0)];
    BS_REAL R1_c0001 = R1_data[NMS_IDX(i0, j0, k0, l1)];
    BS_REAL R1_c0010 = R1_data[NMS_IDX(i0, j0, k1, l0)];
    BS_REAL R1_c0011 = R1_data[NMS_IDX(i0, j0, k1, l1)];
    BS_REAL R1_c0100 = R1_data[NMS_IDX(i0, j1, k0, l0)];
    BS_REAL R1_c0101 = R1_data[NMS_IDX(i0, j1, k0, l1)];
    BS_REAL R1_c0110 = R1_data[NMS_IDX(i0, j1, k1, l0)];
    BS_REAL R1_c0111 = R1_data[NMS_IDX(i0, j1, k1, l1)];
    BS_REAL R1_c1000 = R1_data[NMS_IDX(i1, j0, k0, l0)];
    BS_REAL R1_c1001 = R1_data[NMS_IDX(i1, j0, k0, l1)];
    BS_REAL R1_c1010 = R1_data[NMS_IDX(i1, j0, k1, l0)];
    BS_REAL R1_c1011 = R1_data[NMS_IDX(i1, j0, k1, l1)];
    BS_REAL R1_c1100 = R1_data[NMS_IDX(i1, j1, k0, l0)];
    BS_REAL R1_c1101 = R1_data[NMS_IDX(i1, j1, k0, l1)];
    BS_REAL R1_c1110 = R1_data[NMS_IDX(i1, j1, k1, l0)];
    BS_REAL R1_c1111 = R1_data[NMS_IDX(i1, j1, k1, l1)];

    // R2
    BS_REAL R2_c0000 = R2_data[NMS_IDX(i0, j0, k0, l0)];
    BS_REAL R2_c0001 = R2_data[NMS_IDX(i0, j0, k0, l1)];
    BS_REAL R2_c0010 = R2_data[NMS_IDX(i0, j0, k1, l0)];
    BS_REAL R2_c0011 = R2_data[NMS_IDX(i0, j0, k1, l1)];
    BS_REAL R2_c0100 = R2_data[NMS_IDX(i0, j1, k0, l0)];
    BS_REAL R2_c0101 = R2_data[NMS_IDX(i0, j1, k0, l1)];
    BS_REAL R2_c0110 = R2_data[NMS_IDX(i0, j1, k1, l0)];
    BS_REAL R2_c0111 = R2_data[NMS_IDX(i0, j1, k1, l1)];
    BS_REAL R2_c1000 = R2_data[NMS_IDX(i1, j0, k0, l0)];
    BS_REAL R2_c1001 = R2_data[NMS_IDX(i1, j0, k0, l1)];
    BS_REAL R2_c1010 = R2_data[NMS_IDX(i1, j0, k1, l0)];
    BS_REAL R2_c1011 = R2_data[NMS_IDX(i1, j0, k1, l1)];
    BS_REAL R2_c1100 = R2_data[NMS_IDX(i1, j1, k0, l0)];
    BS_REAL R2_c1101 = R2_data[NMS_IDX(i1, j1, k0, l1)];
    BS_REAL R2_c1110 = R2_data[NMS_IDX(i1, j1, k1, l0)];
    BS_REAL R2_c1111 = R2_data[NMS_IDX(i1, j1, k1, l1)];

    // R3
    BS_REAL R3_c0000 = R3_data[NMS_IDX(i0, j0, k0, l0)];
    BS_REAL R3_c0001 = R3_data[NMS_IDX(i0, j0, k0, l1)];
    BS_REAL R3_c0010 = R3_data[NMS_IDX(i0, j0, k1, l0)];
    BS_REAL R3_c0011 = R3_data[NMS_IDX(i0, j0, k1, l1)];
    BS_REAL R3_c0100 = R3_data[NMS_IDX(i0, j1, k0, l0)];
    BS_REAL R3_c0101 = R3_data[NMS_IDX(i0, j1, k0, l1)];
    BS_REAL R3_c0110 = R3_data[NMS_IDX(i0, j1, k1, l0)];
    BS_REAL R3_c0111 = R3_data[NMS_IDX(i0, j1, k1, l1)];
    BS_REAL R3_c1000 = R3_data[NMS_IDX(i1, j0, k0, l0)];
    BS_REAL R3_c1001 = R3_data[NMS_IDX(i1, j0, k0, l1)];
    BS_REAL R3_c1010 = R3_data[NMS_IDX(i1, j0, k1, l0)];
    BS_REAL R3_c1011 = R3_data[NMS_IDX(i1, j0, k1, l1)];
    BS_REAL R3_c1100 = R3_data[NMS_IDX(i1, j1, k0, l0)];
    BS_REAL R3_c1101 = R3_data[NMS_IDX(i1, j1, k0, l1)];
    BS_REAL R3_c1110 = R3_data[NMS_IDX(i1, j1, k1, l0)];
    BS_REAL R3_c1111 = R3_data[NMS_IDX(i1, j1, k1, l1)];

    BS_REAL a0 = (one - tx), a1 = tx;
    BS_REAL b0 = (one - ty), b1 = ty;
    BS_REAL c0 = (one - tz), c1 = tz;
    BS_REAL d0 = (one - tw), d1 = tw;

    // R1
    BS_REAL R1 = R1_c0000 * a0 * b0 * c0 * d0 + R1_c0001 * a0 * b0 * c0 * d1 +
                R1_c0010 * a0 * b0 * c1 * d0 + R1_c0011 * a0 * b0 * c1 * d1 +
                R1_c0100 * a0 * b1 * c0 * d0 + R1_c0101 * a0 * b1 * c0 * d1 +
                R1_c0110 * a0 * b1 * c1 * d0 + R1_c0111 * a0 * b1 * c1 * d1 +
                R1_c1000 * a1 * b0 * c0 * d0 + R1_c1001 * a1 * b0 * c0 * d1 +
                R1_c1010 * a1 * b0 * c1 * d0 + R1_c1011 * a1 * b0 * c1 * d1 +
                R1_c1100 * a1 * b1 * c0 * d0 + R1_c1101 * a1 * b1 * c0 * d1 +
                R1_c1110 * a1 * b1 * c1 * d0 + R1_c1111 * a1 * b1 * c1 * d1;

    // R2
    BS_REAL R2 = R2_c0000 * a0 * b0 * c0 * d0 + R2_c0001 * a0 * b0 * c0 * d1 +
                R2_c0010 * a0 * b0 * c1 * d0 + R2_c0011 * a0 * b0 * c1 * d1 +
                R2_c0100 * a0 * b1 * c0 * d0 + R2_c0101 * a0 * b1 * c0 * d1 +
                R2_c0110 * a0 * b1 * c1 * d0 + R2_c0111 * a0 * b1 * c1 * d1 +
                R2_c1000 * a1 * b0 * c0 * d0 + R2_c1001 * a1 * b0 * c0 * d1 +
                R2_c1010 * a1 * b0 * c1 * d0 + R2_c1011 * a1 * b0 * c1 * d1 +
                R2_c1100 * a1 * b1 * c0 * d0 + R2_c1101 * a1 * b1 * c0 * d1 +
                R2_c1110 * a1 * b1 * c1 * d0 + R2_c1111 * a1 * b1 * c1 * d1;

    // R3
    BS_REAL R3 = R3_c0000 * a0 * b0 * c0 * d0 + R3_c0001 * a0 * b0 * c0 * d1 +
                R3_c0010 * a0 * b0 * c1 * d0 + R3_c0011 * a0 * b0 * c1 * d1 +
                R3_c0100 * a0 * b1 * c0 * d0 + R3_c0101 * a0 * b1 * c0 * d1 +
                R3_c0110 * a0 * b1 * c1 * d0 + R3_c0111 * a0 * b1 * c1 * d1 +
                R3_c1000 * a1 * b0 * c0 * d0 + R3_c1001 * a1 * b0 * c0 * d1 +
                R3_c1010 * a1 * b0 * c1 * d0 + R3_c1011 * a1 * b0 * c1 * d1 +
                R3_c1100 * a1 * b1 * c0 * d0 + R3_c1101 * a1 * b1 * c0 * d1 +
                R3_c1110 * a1 * b1 * c1 * d0 + R3_c1111 * a1 * b1 * c1 * d1;

    res.R1 = R1;
    res.R2 = R2;
    res.R3 = R3;

#undef NMS_IDX

    return res;
}


//========================================================//
// --- Inelastic Neutrino Scattering off muons Kernel --- //
// ---        via DIRECT KERNEL INTERPOLATION         --- //
//========================================================//
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
 *           kernel:    Inelastic Neutrino Scattering off muons Kernel  [nm^3 s^-1]
 */


// Calculates and saves the neutrino muon scattering in and out kernel
// @TODO: GENERALIZE TO 6 NEUTRINO SPECIES (nue, nut are equal)
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticNMSKernels_DirectInterp(InelasticScattKernelParams* kernel_params,
                                        MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega;
    const BS_REAL wp                   = kernel_params->omega_prime;
    const BS_REAL mu_mu                = eos_params->mu_mu;    
    const BS_REAL det_bal_arg          = (wp - w) / T;
    const BS_REAL det_bal              = SafeExp(det_bal_arg);

    MyKernelOutput nms_kernel;

    // Interpolation of the 4D table
    const NMS_KernelResult R123_values = NMS_direct_interpolator(T, mu_mu, w, wp);

    // Neutrino-flavor dependent kernel
    //nu_e, nu_t
    nms_kernel.abs[id_nue] = kBS_NMS_Conv_Const * (kBS_NMS_C1_nue * R123_values.R1 + 
             kBS_NMS_C2_nue * R123_values.R2 + 
             kBS_NMS_C3_nue * R123_values.R3);
    
    //anu_e, anu_t
    nms_kernel.abs[id_anue] = kBS_NMS_Conv_Const * (kBS_NMS_C2_nue * R123_values.R1 + 
             kBS_NMS_C1_nue * R123_values.R2 + 
             kBS_NMS_C3_nue * R123_values.R3);
    
    //nu_mu==nu_x
    nms_kernel.abs[id_nux] = kBS_NMS_Conv_Const * (kBS_NMS_C1_numu * R123_values.R1 + 
             kBS_NMS_C2_numu * R123_values.R2 + 
             kBS_NMS_C3_numu * R123_values.R3);
    
    //anu_mu==anu_x
    nms_kernel.abs[id_anux] = kBS_NMS_Conv_Const * (kBS_NMS_C2_numu * R123_values.R1 + 
             kBS_NMS_C1_numu * R123_values.R2 + 
             kBS_NMS_C3_numu * R123_values.R3);

    // Absorption and Emission kernels [nm^3 s^-1]
    // nu_e(idx=0), nu_ae(idx=1), nu_x=nu_mu(idx=2), anu_x=anu_mu(idx=3)
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        nms_kernel.em[idx]  = nms_kernel.abs[idx] * det_bal;
    }

    return nms_kernel;
}


#endif // BNS_NURATES_INCLUDE_KERNEL_NMS_DIRECT_INTERP_HPP_
