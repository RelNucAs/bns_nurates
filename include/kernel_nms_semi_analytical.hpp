//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nms_semi_analytical.hpp
//  \brief contains kernels for inelastic neutrino scattering
//         off muons
//
// Computation of inelastic neutrino-muon scattering
// via semi-analytical formula and parameters interpolation.
//
// @TODO: ADD REFERENCE 


#ifndef BNS_NURATES_INCLUDE_KERNEL_NMS_SEMI_ANALYTICAL_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_NMS_SEMI_ANALYTICAL_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"
#include "nms_semi_analytical_table.hpp"
#include <iostream>


/* 3D nearest interpolation of the numerical table of the semi-analytical parameters
 *
 *      Inputs:
 *          T:            temperature                [MeV]
 *          mu:           muon chemical potential    [MeV]
 *          w :           incoming neutrino energy   [MeV]
 *          *table:       RTable-type struct      
 *
 *      Output:
 *          Interpolated values:
 *          alpha: []
 *          beta: []
 *          gamma: []
 *          delta: []
 *          wp_bar: [MeV]
 *          wp_zero: [MeV]
 */

// Struct to store the different parameters
struct NMS_Parameters{
    BS_REAL alpha;
    BS_REAL beta;
    BS_REAL gamma;
    BS_REAL delta;
    BS_REAL wpbar;
    BS_REAL wpzero;
};

KOKKOS_INLINE_FUNCTION
NMS_Parameters NMS_parameters_interpolator(const BS_REAL T, const BS_REAL mu, const BS_REAL w)
{
    NMS_Parameters parameters = {0};

    int i_nearest, j_nearest, k_nearest;

    if (NMSNearest_find_nearest_index(T, NMSParams_T_axis, NMSParams_T_dims, &i_nearest) < 0 or
        NMSNearest_find_nearest_index(mu, NMSParams_mu_axis, NMSParams_mu_dims, &j_nearest) < 0 or
        NMSNearest_find_nearest_index(w, NMSParams_w_axis, NMSParams_w_dims, &k_nearest) < 0)
    {
        return parameters;
    }

// macro to access the 1D array
#define NMSParams_IDX(i, j, k)                                                        \
    ((i * NMSParams_mu_dims + j) * NMSParams_w_dims + k)

    parameters.alpha = alpha_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];
    parameters.beta = beta_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];
    parameters.gamma = gamma_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];
    parameters.delta = delta_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];
    parameters.wpbar = wpbar_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];
    parameters.wpzero = wpzero_numu_data[NMSParams_IDX(i_nearest, j_nearest, k_nearest)];

#undef NMSParams_IDX

    return parameters;
}


// Analytical approximation of numu-muon scattering kernel   [MeV^-2]
BS_REAL NMS_analytical_kernel_nu_mu(const BS_REAL T, const BS_REAL mu,
                                        const BS_REAL w, const BS_REAL wp_in, const BS_REAL u)
{

    // Protection for the elastic case (w=wp), to avoid divergences
    const BS_REAL wp = (std::abs(w - wp_in) <= 1e-6) ? (w * 0.99999) : wp_in;

    const BS_REAL x0 = POW2(w);
    const BS_REAL x1 = POW2(wp);
    const BS_REAL x2 = w + wp;
    const BS_REAL x3 = 2. * w;
    const BS_REAL x4 = x3 * wp;
    const BS_REAL x5 = x0 + x1;
    const BS_REAL x6 = x4 + x5;
    const BS_REAL x7 = POW2(kBS_Mmu);
    const BS_REAL x8 = 2. * kBS_SinThW2 + 1.;
    const BS_REAL x9 = Kokkos::fabs(w - wp);
    const BS_REAL x10 = Kokkos::fabs(x2);
    const BS_REAL x11 = w * wp;
    const BS_REAL x12 = x10 * x11;
    const BS_REAL x13 = x0 * x9;
    const BS_REAL x14 = -x0 * x10 - x1 * x10 + x1 * x9 + x13;
    
    const BS_REAL E_ = 0.5 * (
                            wp - w + Kokkos::sqrt( 
                                (x5 - 2. * x11 * u) * ( 1. + ((2. * POW2(kBS_Mmu)) / (x11 * (1. - u))) )
                            )
                        );
    const BS_REAL eta = mu / T;
    const BS_REAL arg1 = eta - E_ / T;
    const BS_REAL arg2 = eta + (wp - w) / T - E_ / T;
    const BS_REAL F0eta = FDI_0(arg1);
    const BS_REAL F1eta = FDI_p1(arg1);
    const BS_REAL F2eta = FDI_p2(arg1);
    const BS_REAL F0eta2 = FDI_0(arg2);
    const BS_REAL F1eta2 = FDI_p1(arg2);
    const BS_REAL F2eta2 = FDI_p2(arg2);

    const BS_REAL x15 = F0eta - F0eta2;
    const BS_REAL x16 = T * (F1eta - F1eta2);
    const BS_REAL x17 = x2 * x6 * (POW2(E_) * x15 + 2. * E_ * x16 + POW2(T) * (F2eta - F2eta2)) * 
                        (x11 * x9 + x12 + x14);
    const BS_REAL x18 = POW4(wp);
    const BS_REAL x19 = x10 * x18;
    const BS_REAL x20 = x18 * x9;
    const BS_REAL x21 = POW4(w);
    const BS_REAL x22 = x10 * x21;
    const BS_REAL x23 = x21 * x9;
    const BS_REAL x24 = 2. * x1;
    const BS_REAL x25 = x13 * x24;
    const BS_REAL x26 = POW3(w);
    const BS_REAL x27 = 4. * x26;
    const BS_REAL x28 = x10 * wp;
    const BS_REAL x29 = POW3(wp);
    const BS_REAL x30 = x29 * x9;
    const BS_REAL x31 = 2. * wp;
    const BS_REAL x32 = x26 * x9;
    const BS_REAL x33 = -x3 * x30 + x31 * x32;
    const BS_REAL x34 = x6 * (E_ * x15 + x16);
    const BS_REAL x35 = POW6(wp);
    const BS_REAL x36 = x10 * x35;
    const BS_REAL x37 = x35 * x9;
    const BS_REAL x38 = POW6(w);
    const BS_REAL x39 = x38 * x9;
    const BS_REAL x40 = POW5(wp);
    const BS_REAL x41 = POW5(w);
    const BS_REAL x42 = x1 * x21;
    const BS_REAL x43 = Kokkos::fabs( (w - wp) * POW2((w + wp)) );
    const BS_REAL x44 = 45. * x43;
    const BS_REAL x45 = 15. * x43;
    const BS_REAL x46 = x10 * x38;
    const BS_REAL x47 = x10 * x42;
    const BS_REAL x48 = x0 * x20;
    const BS_REAL x49 = w * x29;
    const BS_REAL x50 = x28 * x41;
    const BS_REAL x51 = Kokkos::fabs(x0 - x1);
    const BS_REAL x52 = 15. * x51;
    const BS_REAL x53 = w * x40;
    const BS_REAL x54 = x10 * x53;
    const BS_REAL x55 = x26 * wp;
    const BS_REAL x56 = x0 * x19;
    const BS_REAL x57 = 45. * x51;
    const BS_REAL x58 = x42 * x9;
    const BS_REAL x59 = x32 * wp;
    const BS_REAL x60 = w * x30;
    const BS_REAL x61 = 40. * x7;
    const BS_REAL x62 = 3. * x9;
    const BS_REAL x63 = 30. * x1;
    const BS_REAL x64 = x10 * x49;
    const BS_REAL x65 = 120. * x7;
    const BS_REAL x66 = x0 * x43 * x63 + 160. * x1 * x13 * x7 + 130. * x10 * x26 * x29 - x13 * x51 * x63 
                      + (- x19 + x20 - x22 + x23 - x26 * x28 - x64) * x61 + 6. * x26 * x30 
                      - x41 * x62 * wp - x53 * x62 + x59 * x65 + x60 * x65;

    const BS_REAL numerator = POW2(kBS_Gf0) * T * (
            POW2(kBS_SinThW2) * (
                4. * x15 * x2 * (
                    -x18 * x45 + x20 * x52 - x21 * x44 + x23 * x57 + x36 - x37
                    - 51. * x39 + x44 * x55 - x45 * x49 + 51. * x46 - 185. * x47
                    - 49. * x48 - 9. * x50 + x52 * x60 + 31. * x54 + 45. * x56
                    - x57 * x59 + 101. * x58 + x66
                )
                - 640. * x17
                - 320. * x34 * (x19 - x20 - 3. * x22 + 3. * x23 - x25 + x27 * x28 + x33)
            )
            - 80. * kBS_SinThW2 * x15 * x2 * x6 * x7 * x8 * (4. * x12 + x14 - x4 * x9)
            + POW2(x8) * (
                x15 * x2 * (
                    -x18 * x44 + x20 * x57 - x21 * x45 + x23 * x52 + 51. * x36
                    - 51. * x37 - x39 + x44 * x49 - x45 * x55 + x46 + 45. * x47
                    + 101. * x48 + 31. * x50 + x52 * x59 - 9. * x54 - 185. * x56
                    - x57 * x60 - 49. * x58 + x66
                )
                - 160. * x17
                - 80. * x34 * (3. * x19 - 3. * x20 - x22 + x23 + x25 + x33 - 4. * x64)
            )
        );

    const BS_REAL denominator = 240. * kBS_Pi * x0 * x1 * x2 * x6 * (1. - SafeExp((-w + wp) / T));

    return numerator / denominator;
}

BS_REAL smoothstep(const BS_REAL wp, const BS_REAL wpzero, const BS_REAL epsilon){

    return (3. * POW2( (wp - (wpzero - epsilon)) / (2. * epsilon) ) 
                - 2. * POW3( (wp - (wpzero - epsilon)) / (2. * epsilon) )
    );
}


// Semi-Analytical expression of the NMS kernel  [MeV^-2]
BS_REAL NMS_SemiAnalyticalKernel(const BS_REAL T, const BS_REAL mu, const BS_REAL w, const BS_REAL wp,
                        const BS_REAL alpha, const BS_REAL beta, const BS_REAL gamma, const BS_REAL delta,
                        const BS_REAL wpbar, const BS_REAL wpzero)
{
    BS_REAL min_val = std::min(wpbar, w);
    BS_REAL max_val = std::max(wpbar, w);
    BS_REAL epsilon = 0.2 * wpzero;

    if (wp <= (wpzero - epsilon)){
        return NMS_analytical_kernel_nu_mu(T, mu, w, wp, -0.4);
    }
    else if (wp > (wpzero - epsilon) && wp < (wpzero + epsilon)){
        return ( NMS_analytical_kernel_nu_mu(T, mu, w, wp, 0.0) * smoothstep(wp, wpzero, epsilon) 
                + NMS_analytical_kernel_nu_mu(T, mu, w, wp, -0.4) * (1 - smoothstep(wp, wpzero, epsilon)) ); 
    }
    else if (wp >= (wpzero + epsilon) && wp <= min_val){
        return NMS_analytical_kernel_nu_mu(T, mu, w, wp, 0.0);
    }
    else if (wp > min_val && wp <= max_val){
        BS_REAL inv_min = 1. / min_val;
        BS_REAL prefactor1 = NMS_analytical_kernel_nu_mu(T, mu, w, min_val, 0.0);
        return prefactor1 * Kokkos::pow(wp * inv_min, gamma) * SafeExp(delta * (min_val - wp));
    }
    else if (wp > max_val){
        BS_REAL inv_min = 1. / min_val;
        BS_REAL inv_max = 1. / max_val;
        BS_REAL prefactor2 = NMS_analytical_kernel_nu_mu(T, mu, w, min_val, 0.0) 
                            * Kokkos::pow(max_val * inv_min, gamma) * SafeExp(delta * (min_val - max_val));
        return prefactor2 * Kokkos::pow(wp * inv_max, alpha) * SafeExp(beta * (max_val - wp));
    }
    else{  //fallback
        return 0.0;
    }
}


// Calculates and saves the neutrino muon scattering in and out kernel
// @TODO: ADD PARAMETER INTERPOLATION FOR OTHER FLAVORS
// @TODO: GENERALIZE TO 6 NEUTRINO SPECIES (nue, nut are equal)
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticNMSKernels_SemiAnalytical(InelasticScattKernelParams* kernel_params,
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
    const NMS_Parameters parameters_values = NMS_parameters_interpolator(T, mu_mu, w);

    // Neutrino-flavor dependent kernel
    // nu_e
    nms_kernel.abs[id_nue] = 0.0;
    
    // anu_e
    nms_kernel.abs[id_anue] = 0.0;
    
    // nu_m
    nms_kernel.abs[id_num] = kBS_NMS_Conv_Const * NMS_SemiAnalyticalKernel(T, mu_mu, w, wp, 
                                parameters_values.alpha, parameters_values.beta, parameters_values.gamma, 
                                parameters_values.delta, parameters_values.wpbar, parameters_values.wpzero);
    
    // anu_m
    nms_kernel.abs[id_anum] = 0.0;

    // nu_t
    nms_kernel.abs[id_nut] = nms_kernel.abs[id_nue];

    // anu_t
    nms_kernel.abs[id_anut] = nms_kernel.abs[id_anue];

    // Absorption and Emission kernels [nm^3 s^-1]
    // nu_e(idx=0), nu_ae(idx=1), nu_x=nu_mu(idx=2), anu_x=anu_mu(idx=3)
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        nms_kernel.em[idx]  = nms_kernel.abs[idx] * det_bal;
    }

    return nms_kernel;
}


#endif // BNS_NURATES_INCLUDE_KERNEL_NMS_SEMI_ANALYTICAL_HPP_
