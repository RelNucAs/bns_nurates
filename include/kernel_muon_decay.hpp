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
KOKKOS_INLINE_FUNCTION
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


// Analytical Expression of Inverse Anti-Muon Decay kernel [MeV^-2]
KOKKOS_INLINE_FUNCTION
BS_REAL InverseAntiMuonDecayKernel_Analytical(const BS_REAL T, const BS_REAL mu_pos,
                                          const BS_REAL u_bar, const BS_REAL w_anumu, const BS_REAL w_nue)
{

    BS_REAL x2 = w_anumu * w_nue;
    BS_REAL x5 = 1. - u_bar;
    BS_REAL x15 = x2 * x5;
    BS_REAL x17 = kBS_squared_mass_diff / x15 - 1.;
    if (x17 < 0){
        return 0.;
    }
    BS_REAL x0 = POW2(w_nue);
    BS_REAL x1 = POW2(w_anumu);
    BS_REAL x3 = 2. * u_bar * x2 + x1;
    BS_REAL x4 = x0 + x3;
    BS_REAL x6 = POW2(x5);
    BS_REAL x7 = POW3(w_anumu);
    BS_REAL x8 = POW3(w_nue);
    BS_REAL x9 = u_bar + 2.;
    BS_REAL x10 = w_anumu * x0;
    BS_REAL x11 = x1 * w_nue;
    BS_REAL x12 = kBS_squared_mass_diff * x5;
    BS_REAL x13 = u_bar * x1;
    BS_REAL x14 = POW2(u_bar);
    BS_REAL x16 = w_anumu + w_nue;
    BS_REAL x18 = x16 * x17;
    BS_REAL x19 = Kokkos::sqrt(x4 * Kokkos::fabs(POW2(x17) - 2. * POW2(kBS_Me) / x15));
    BS_REAL E_plus = 0.5 * (x18 + x19);
    BS_REAL E_minus = 0.5 * (x18 - x19);
    BS_REAL e_plus = E_plus / T;
    BS_REAL e_minus = E_minus / T;
    BS_REAL eta = mu_pos / T;
    BS_REAL arg_plus = eta - e_plus;
    BS_REAL arg_minus = eta - e_minus;
    BS_REAL F0_minus = FDI_0(arg_minus);
    BS_REAL F0_plus = FDI_0(arg_plus);
    BS_REAL F1_minus = FDI_p1(arg_minus);
    BS_REAL F1_plus = FDI_p1(arg_plus);
    BS_REAL F2_minus = FDI_p2(arg_minus);
    BS_REAL F2_plus = FDI_p2(arg_plus);

    return (kBS_MuonDecay_C * (
                POW3(T) * x2 * x6 * (w_anumu * w_nue * (u_bar + 3.) - x0 - x1) * 
                    (F0_minus * POW2(e_minus) - F0_plus * POW2(e_plus) 
                        + 2. * F1_minus * e_minus - 2. * F1_plus * e_plus + F2_minus - F2_plus) 
                + POW2(T) * (x11 * x6 * (x0 * (3. * u_bar + 1.) - 2. * x1 + x2 * (3. - u_bar)) 
                            + x12 * (-x10 * x9 - x11 * x9 + x7 + x8)) * 
                         (F0_minus * e_minus - F0_plus * e_plus + F1_minus - F1_plus) 
                + T * (F0_minus - F0_plus) * (
                    POW2(kBS_squared_mass_diff) * (u_bar * x0 + x13 + x2 * (0.5 * x14 + 1.5)) 
                    + 0.5 * POW2(kBS_Me) * x2 * x4 * (1. - x14) 
                    + w_anumu * x12 * (-u_bar * x8 + x10 * (x14 - 2.) + x13 * w_nue + x7) 
                    - x6 * x7 * w_nue * (x0 * (1.5 * x14 - 0.5) + x3)
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
 * 
 * ATTENTION:
 * 
 *      We need to give the right values of w_mu and w_e:
 *          muon decay      : w_mu = w_numu,  w_e = w_anue
 *          anti-muon decay : w_mu = w_anumu, w_e = w_nue
 */

KOKKOS_INLINE_FUNCTION
Decay_ubar_int16 ubar_int16_interpolator(const BS_REAL T, const BS_REAL mu, 
                                         const BS_REAL w_mu, const BS_REAL w_e,
                                        const short* tab_mu, const short* tab_amu)
{

    Decay_ubar_int16 out;

// Keep mwe.cpp and CPU tests working
#if !defined(__CUDA_ARCH__) && !defined(__SYCL_DEVICE_ONLY__) && !defined(__HIP_DEVICE_COMPILE__)
    if (tab_mu  == nullptr) tab_mu  = muon_decay_ubar_int16_data;
    if (tab_amu == nullptr) tab_amu = anti_muon_decay_ubar_int16_data;
#endif

    // If there are no tables, return an interpolation error
    if (tab_mu == nullptr || tab_amu == nullptr){
        out.ubar_int16_muon_decay       = -65536.;
        out.ubar_int16_anti_muon_decay  = -65536.;
        return out;
    }

    int i0, j0, k0, l0;
    BS_REAL tx, ty, tz, tw;

    MuonDecay_LogAxisBracket(T, MuonDecay_T_min, MuonDecay_T_max, MuonDecay_T_dims, i0, tx);
    MuonDecay_LinAxisBracket(mu, MuonDecay_mu_min, MuonDecay_mu_max, MuonDecay_mu_dims, j0, ty);
    MuonDecay_LogAxisBracket(w_mu, MuonDecay_w_mu_min, MuonDecay_w_mu_max, MuonDecay_w_mu_dims, k0, tz);
    MuonDecay_LogAxisBracket(w_e, MuonDecay_w_e_min, MuonDecay_w_e_max, MuonDecay_w_e_dims, l0, tw);

    BS_REAL res_mu = 0., res_amu = 0.;

    // Consider 16 vertices of the hypercube
    for (int v = 0; v < 16; ++v){
        // Use binary representation to extract cube vertices (eg. 1=0001)
        const int a = (v >> 3) & 1; //fourth bit, T dimension
        const int b = (v >> 2) & 1; //third bit, mu dimension
        const int c = (v >> 1) & 1; //second bit, w_mu dimension
        const int d = v & 1;        //first bit, w_e dimension

        // Depending on (a,b,c,d), compute the weight of the vertex
        const BS_REAL wgt = (a ? tx : 1. - tx) * (b ? ty : 1. - ty)
                          * (c ? tz : 1. - tz) * (d ? tw : 1. - tw);
        
        // Index for the ubar vector
        const int idx = (((i0 + a) * MuonDecay_mu_dims + (j0 + b)) * MuonDecay_w_mu_dims
                        + (k0 + c)) * MuonDecay_w_e_dims + (l0 + d);

        // Sum weighted values at each vertex
        res_mu += wgt * MuonDecay_Load(tab_mu, idx);
        res_amu += wgt * MuonDecay_Load(tab_amu, idx);
    }

    // Return interpolated values
    out.ubar_int16_muon_decay      = res_mu;
    out.ubar_int16_anti_muon_decay = res_amu;

    return out;
}


//================================================//
// ---         (Anti)Muon Decay Kernel        --- //
// ---          via ubar interpolation        --- //
//================================================//
/* Interpolate ubar grid + Analytical expression
 *
 * NOTE: points outside the ranges of the 4D table are mapped to the borders!
 *
 *      Inputs:
 *          T:          temperature                    [MeV]
 *          mu:         electron chemical potential    [MeV]
 *          w_1:        muon (anti)neutrino energy     [MeV]
 *          w_2:        electron (anti)neutrino energy [MeV]
 *
 *      Output:
 *          kernel:     (anti)muon-decay Kernel  [nm^3 s^-1]
 *
 * ATTENTION:
 *    muon decay      : w_mu = w_numu,  w_e = w_anue
 *    anti-muon decay : w_mu = w_anumu, w_e = w_nue
 * 
*/

// Calculates and saves the (anti)muon decay in and out kernel
KOKKOS_INLINE_FUNCTION
MyKernelOutput MuonDecayKernels(MuonDecayKernelParams* kernel_params,
                                MyEOSParams* eos_params)
{

    constexpr BS_REAL scale_factor = 1. / 32767.;

    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w_mu                 = kernel_params->omega_mu;
    const BS_REAL w_e                  = kernel_params->omega_e;
    const BS_REAL mu_el                = eos_params->mu_e;  
    const BS_REAL mu_pos               = - mu_el;  
    const BS_REAL mu_mu                = eos_params->mu_mu;
    const BS_REAL mu_antimu            = - mu_mu;   
    const BS_REAL det_bal_arg_mudec    = (mu_mu - mu_el - w_mu - w_e) / T;
    const BS_REAL det_bal_mudec        = SafeExp(det_bal_arg_mudec);
    const BS_REAL det_bal_arg_amudec   = (mu_antimu - mu_pos - w_mu - w_e) / T;
    const BS_REAL det_bal_amudec       = SafeExp(det_bal_arg_amudec);

    const short* tab_mu = kernel_params->ubar_mu_table;
    const short* tab_amu = kernel_params->ubar_amu_table;

    MyKernelOutput mudec_kernel = {0.};

    // Spectral integrals related to numu are performed along anue 
    // (anue is the varying energy). But this is wrong when computing 
    // spectral quantities related to anue. Thus, we need to exchange 
    // numu<-->anue in the calculation of the kernel, 
    // so that numu becomes the varying quantity along which we integrate.
    // The same is valid regarding anti-muon decay: in this case we need to exchange
    // anumu<-->nue.

    // Instead, when computing gray quantities this change is not needed,
    // and we consider only the id_num(id_anum for anti-muon decay) term.

    // Interpolate u_bar in both cases
    Decay_ubar_int16 out_wmu_we = ubar_int16_interpolator(T, mu_el, w_mu, w_e, tab_mu, tab_amu);
    Decay_ubar_int16 out_we_wmu = ubar_int16_interpolator(T, mu_el, w_e, w_mu, tab_mu, tab_amu); 

    // Muon Decay (mu <--> el + numu + anue)
    // nu_mu
    BS_REAL ubar_numu = scale_factor * out_wmu_we.ubar_int16_muon_decay;
    BS_REAL InverseMuonDecay_Kernel_numu = (ubar_numu < -1.0) 
                                            ? 0. 
                                            : InverseMuonDecayKernel_Analytical(T, mu_el, 
                                                                            ubar_numu, w_mu, w_e);
    // anu_e
    BS_REAL ubar_anue = scale_factor * out_we_wmu.ubar_int16_muon_decay;
    BS_REAL InverseMuonDecay_Kernel_anue = (ubar_anue < -1.0) 
                                            ? 0. 
                                            : InverseMuonDecayKernel_Analytical(T, mu_el, 
                                                                            ubar_anue, w_e, w_mu);

    // Anti-muon Decay  (mu^+ <--> e^+ + anumu + nue)
    // anu_mu
    BS_REAL ubar_anumu = scale_factor * out_wmu_we.ubar_int16_anti_muon_decay;
    BS_REAL InverseMuonDecay_Kernel_anumu = (ubar_anumu < -1.0) 
                                            ? 0. 
                                            : InverseAntiMuonDecayKernel_Analytical(T, mu_pos, 
                                                                            ubar_anumu, w_mu, w_e);
    // nu_e
    BS_REAL ubar_nue = scale_factor * out_we_wmu.ubar_int16_anti_muon_decay;
    BS_REAL InverseMuonDecay_Kernel_nue = (ubar_nue < -1.0) 
                                            ? 0. 
                                            : InverseAntiMuonDecayKernel_Analytical(T, mu_pos, 
                                                                            ubar_nue, w_e, w_mu);

    // Absorption and Emission kernels [nm^3 s^-1]
    // nu_e
    mudec_kernel.abs[id_nue] = kBS_NMAS_Conv_Const * InverseMuonDecay_Kernel_nue;
    mudec_kernel.em[id_nue] = mudec_kernel.abs[id_nue] * det_bal_amudec;
    // anu_e
    mudec_kernel.abs[id_anue] = kBS_NMAS_Conv_Const * InverseMuonDecay_Kernel_anue;
    mudec_kernel.em[id_anue] = mudec_kernel.abs[id_anue] * det_bal_mudec;
    // nu_m
    mudec_kernel.abs[id_num] = kBS_NMAS_Conv_Const * InverseMuonDecay_Kernel_numu;
    mudec_kernel.em[id_num] = mudec_kernel.abs[id_num] * det_bal_mudec;
    // anu_m
    mudec_kernel.abs[id_anum] = kBS_NMAS_Conv_Const * InverseMuonDecay_Kernel_anumu;
    mudec_kernel.em[id_anum] = mudec_kernel.abs[id_anum] * det_bal_amudec;

    return mudec_kernel;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_MUON_DECAY_HPP_
