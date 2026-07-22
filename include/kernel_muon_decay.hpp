//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_muon_decay.hpp
//  \brief contains kernel for muon decay
//
// Computation of muon decay kernel via u_bar grid
// interpolation, neglecting muon blocking factor.
//
// @TODO: ADD REFERENCE 


#ifndef BNS_NURATES_INCLUDE_KERNEL_MUON_DECAY_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_MUON_DECAY_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"
#include "muon_decay_ubar_table.hpp"


// Analytical Expression of Inverse Muon Decay kernel [MeV^-2]
BS_REAL InverseMuonDecayKernel_Analytical(const BS_REAL T, const BS_REAL mu_el, 
                                          const BS_REAL u_bar, const BS_REAL wnu, const BS_REAL wanu)
{

    BS_REAL x2 = wnu * wanu;
    BS_REAL x5 = 1. - u_bar;
    BS_REAL x15 = x2 * x5;
    BS_REAL x17 = kBS_squared_mass_diff / x15 - 1.;
    if (x17 < 0){
        return 0.;
    }
    BS_REAL x0 = POW2(wanu);
    BS_REAL x1 = POW2(wnu);
    BS_REAL x3 = 2. * u_bar * x2 + x1;
    BS_REAL x4 = x0 + x3;
    BS_REAL x6 = POW2(x5);
    BS_REAL x7 = POW3(wnu);
    BS_REAL x8 = POW3(wanu);
    BS_REAL x9 = u_bar + 2.;
    BS_REAL x10 = wnu * x0;
    BS_REAL x11 = x1 * wanu;
    BS_REAL x12 = kBS_squared_mass_diff * x5;
    BS_REAL x13 = u_bar * x1;
    BS_REAL x14 = POW2(u_bar);
    BS_REAL x16 = wnu + wanu;
    BS_REAL x18 = x16 * x17;
    BS_REAL x19 = Kokkos::sqrt(x4 * Kokkos::fabs(POW2(x17) - 2. * POW2(kBS_Me) / x15));
    BS_REAL E_plus = 0.5 * (x18 + x19);
    BS_REAL E_minus = 0.5 * (x18 - x19);
    BS_REAL e_plus = E_plus / T;
    BS_REAL e_minus = E_minus / T;
    BS_REAL eta = mu_el / T;
    BS_REAL arg_plus = eta - e_plus;
    BS_REAL arg_minus = eta - e_minus;
    BS_REAL F0_minus = FDI_0(arg_minus);
    BS_REAL F0_plus = FDI_0(arg_plus);
    BS_REAL F1_minus = FDI_p1(arg_minus);
    BS_REAL F1_plus = FDI_p1(arg_plus);
    BS_REAL F2_minus = FDI_p2(arg_minus);
    BS_REAL F2_plus = FDI_p2(arg_plus);

    return (kBS_MuonDecay_C * (
                POW3(T) * x2 * x6 * (wnu * wanu * (u_bar + 3.) - x0 - x1) * 
                    (F0_minus * POW2(e_minus) - F0_plus * POW2(e_plus) 
                        + 2. * F1_minus * e_minus - 2. * F1_plus * e_plus + F2_minus - F2_plus) 
                + POW2(T) * (x11 * x6 * (x0 * (3. * u_bar + 1.) - 2. * x1 + x2 * (3. - u_bar)) 
                            + x12 * (-x10 * x9 - x11 * x9 + x7 + x8)) * 
                         (F0_minus * e_minus - F0_plus * e_plus + F1_minus - F1_plus) 
                + T * (F0_minus - F0_plus) * (
                    POW2(kBS_squared_mass_diff) * (u_bar * x0 + x13 + x2 * (0.5 * x14 + 1.5)) 
                    + 0.5 * POW2(kBS_Me) * x2 * x4 * (1. - x14) 
                    + wnu * x12 * (-u_bar * x8 + x10 * (x14 - 2.) + x13 * wanu + x7) 
                    - x6 * x7 * wanu * (x0 * (1.5 * x14 - 0.5) + x3)
                    )
                ) / Kokkos::pow(x4, 2.5)
            );
}

/* 4D linear interpolation of the numerical table of u_bar(int16)
 *
 *      Inputs:
 *          T:            temperature                     [MeV]
 *          mu:           electron chemical potential         [MeV]
 *          w_numu :      muon neutrino energy            [MeV]
 *          w_anue :      electron anti-neutrino energy   [MeV]
 *          *table:       RTable-type struct      
 *
 *      Output:
 *          u_bar:        cosine of the angle between neutrinos (int16)        []
 */

KOKKOS_INLINE_FUNCTION
BS_REAL muon_decay_ubar_int16_interpolator(const BS_REAL T, const BS_REAL mu, const BS_REAL w_numu,
                                                const BS_REAL w_anue)
{

    int i0, i1, j0, j1, k0, k1, l0, l1;
    BS_REAL tx, ty, tz, tw;

    if (MuonLinear_find_bracketing_indices(T, MuonDecay_T_axis, MuonDecay_T_dims, &i0, &i1, &tx) < 0 or
        MuonLinear_find_bracketing_indices(mu, MuonDecay_mu_axis, MuonDecay_mu_dims, &j0, &j1, &ty) < 0 or
        MuonLinear_find_bracketing_indices(w_numu, MuonDecay_wnumu_axis, MuonDecay_wnumu_dims, &k0, &k1, &tz) < 0 or
        MuonLinear_find_bracketing_indices(w_anue, MuonDecay_wanue_axis, MuonDecay_wanue_dims, &l0, &l1, &tw) < 0)
    {
        return -2.;
    }

// macro to access the 1D array
#define MuonDecay_IDX(i, j, k, l)                                                        \
    (((i * MuonDecay_mu_dims + j) * MuonDecay_wnumu_dims + k) * MuonDecay_wanue_dims + l)

    // extract the 16 vertexes of the hypercube for ubar_int16
    BS_REAL ubar_int16_c0000 = ubar_int16_data[MuonDecay_IDX(i0, j0, k0, l0)];
    BS_REAL ubar_int16_c0001 = ubar_int16_data[MuonDecay_IDX(i0, j0, k0, l1)];
    BS_REAL ubar_int16_c0010 = ubar_int16_data[MuonDecay_IDX(i0, j0, k1, l0)];
    BS_REAL ubar_int16_c0011 = ubar_int16_data[MuonDecay_IDX(i0, j0, k1, l1)];
    BS_REAL ubar_int16_c0100 = ubar_int16_data[MuonDecay_IDX(i0, j1, k0, l0)];
    BS_REAL ubar_int16_c0101 = ubar_int16_data[MuonDecay_IDX(i0, j1, k0, l1)];
    BS_REAL ubar_int16_c0110 = ubar_int16_data[MuonDecay_IDX(i0, j1, k1, l0)];
    BS_REAL ubar_int16_c0111 = ubar_int16_data[MuonDecay_IDX(i0, j1, k1, l1)];
    BS_REAL ubar_int16_c1000 = ubar_int16_data[MuonDecay_IDX(i1, j0, k0, l0)];
    BS_REAL ubar_int16_c1001 = ubar_int16_data[MuonDecay_IDX(i1, j0, k0, l1)];
    BS_REAL ubar_int16_c1010 = ubar_int16_data[MuonDecay_IDX(i1, j0, k1, l0)];
    BS_REAL ubar_int16_c1011 = ubar_int16_data[MuonDecay_IDX(i1, j0, k1, l1)];
    BS_REAL ubar_int16_c1100 = ubar_int16_data[MuonDecay_IDX(i1, j1, k0, l0)];
    BS_REAL ubar_int16_c1101 = ubar_int16_data[MuonDecay_IDX(i1, j1, k0, l1)];
    BS_REAL ubar_int16_c1110 = ubar_int16_data[MuonDecay_IDX(i1, j1, k1, l0)];
    BS_REAL ubar_int16_c1111 = ubar_int16_data[MuonDecay_IDX(i1, j1, k1, l1)];

    BS_REAL a0 = (one - tx), a1 = tx;
    BS_REAL b0 = (one - ty), b1 = ty;
    BS_REAL c0 = (one - tz), c1 = tz;
    BS_REAL d0 = (one - tw), d1 = tw;

    BS_REAL ubar_int16 = ubar_int16_c0000 * a0 * b0 * c0 * d0 + ubar_int16_c0001 * a0 * b0 * c0 * d1 +
                         ubar_int16_c0010 * a0 * b0 * c1 * d0 + ubar_int16_c0011 * a0 * b0 * c1 * d1 +
                         ubar_int16_c0100 * a0 * b1 * c0 * d0 + ubar_int16_c0101 * a0 * b1 * c0 * d1 +
                         ubar_int16_c0110 * a0 * b1 * c1 * d0 + ubar_int16_c0111 * a0 * b1 * c1 * d1 +
                         ubar_int16_c1000 * a1 * b0 * c0 * d0 + ubar_int16_c1001 * a1 * b0 * c0 * d1 +
                         ubar_int16_c1010 * a1 * b0 * c1 * d0 + ubar_int16_c1011 * a1 * b0 * c1 * d1 +
                         ubar_int16_c1100 * a1 * b1 * c0 * d0 + ubar_int16_c1101 * a1 * b1 * c0 * d1 +
                         ubar_int16_c1110 * a1 * b1 * c1 * d0 + ubar_int16_c1111 * a1 * b1 * c1 * d1;

#undef MuonDecay_IDX

    return ubar_int16;
}


//==========================================//
// ---         Muon Decay Kernel        --- //
// ---      via ubar interpolation      --- //
//==========================================//
/* Interpolate ubar grid + Analytical expression
 *
 * NOTE: this function is just for points inside the ranges of the 4D table!
 *
 *      Inputs:
 *          T:          temperature                    [MeV]
 *          mu:         electron chemical potential    [MeV]
 *          w_numu:     muon neutrino energy           [MeV]
 *          w_anue:     electron anti-neutrino energy  [MeV]
 *
 *      Output:
 *          kernel:     muon-decay Kernel  [nm^3 s^-1]
 */


// Calculates and saves the muon decay in and out kernel
KOKKOS_INLINE_FUNCTION
MyKernelOutput MuonDecayKernels(MuonDecayKernelParams* kernel_params,
                                MyEOSParams* eos_params)
{

    constexpr BS_REAL scale_factor = 1. / 32767.;

    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w_numu               = kernel_params->omega_numu;
    const BS_REAL w_anue               = kernel_params->omega_anue;
    const BS_REAL mu_el                = eos_params->mu_e;    
    const BS_REAL mu_mu                = eos_params->mu_mu;   
    const BS_REAL det_bal_arg          = (mu_mu - mu_el - w_numu - w_anue) / T;
    const BS_REAL det_bal              = SafeExp(det_bal_arg);

    MyKernelOutput mudec_kernel = {0.};

    // Considering that the integral is performed along anue (anue is the varying energy), 
    // to compute the quantities related to anue, we need to exchange numu<-->anue, so that numu
    // becomes the varying quantity along which we integrate.

    // numu
    BS_REAL ubar_numu = scale_factor * muon_decay_ubar_int16_interpolator(T, mu_el, w_numu, w_anue);
    BS_REAL InverseMuonDecay_Kernel_numu = (ubar_numu < -1.0) 
                                            ? 0. 
                                            : InverseMuonDecayKernel_Analytical(T, mu_el, 
                                                                            ubar_numu, w_numu, w_anue);

    // anue
    BS_REAL ubar_anue = scale_factor * muon_decay_ubar_int16_interpolator(T, mu_el, w_anue, w_numu);
    BS_REAL InverseMuonDecay_Kernel_anue = (ubar_anue < -1.0) 
                                            ? 0. 
                                            : InverseMuonDecayKernel_Analytical(T, mu_el, 
                                                                            ubar_anue, w_anue, w_numu);

    // Absorption and Emission kernels [nm^3 s^-1]
    // anu_e
    mudec_kernel.abs[id_anue] = kBS_NMS_Conv_Const * InverseMuonDecay_Kernel_anue;
    mudec_kernel.em[id_anue] = mudec_kernel.abs[id_anue] * det_bal;
    
    // num
    mudec_kernel.abs[id_num] = kBS_NMS_Conv_Const * InverseMuonDecay_Kernel_numu;
    mudec_kernel.em[id_num] = mudec_kernel.abs[id_num] * det_bal;

    return mudec_kernel;
}


#endif // BNS_NURATES_INCLUDE_KERNEL_NMS_DIRECT_INTERP_HPP_
