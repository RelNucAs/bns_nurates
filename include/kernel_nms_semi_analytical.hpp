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


/* 3D linear interpolation of the numerical table of the semi-analytical parameters
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
    constexpr BS_REAL one = 1;

    NMS_Parameters parameters = {0};

    int i0, i1, j0, j1, k0, k1;
    BS_REAL tx, ty, tz;

    if (NMS_find_bracketing_indices(T, NMSParams_T_axis, NMSParams_T_dims, &i0, &i1, &tx) < 0 or
        NMS_find_bracketing_indices(mu, NMSParams_mu_axis, NMSParams_mu_dims, &j0, &j1, &ty) < 0 or
        NMS_find_bracketing_indices(w, NMSParams_w_axis, NMSParams_w_dims, &k0, &k1, &tz) < 0)
    {
        return parameters;
    }

// macro to access the 1D array
#define NMSParams_IDX(i, j, k)                                                        \
    ((i * NMS_mu_dims + j) * NMS_w_dims + k)

    // extract the 8 vertexes of the hypercube for the parameters
    // Alpha
    BS_REAL alpha_numu_c000 = alpha_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL alpha_numu_c001 = alpha_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL alpha_numu_c010 = alpha_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL alpha_numu_c011 = alpha_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL alpha_numu_c100 = alpha_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL alpha_numu_c101 = alpha_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL alpha_numu_c110 = alpha_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL alpha_numu_c111 = alpha_numu_data[NMSParams_IDX(i1, j1, k1)];

    // Beta
    BS_REAL beta_numu_c000 = beta_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL beta_numu_c001 = beta_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL beta_numu_c010 = beta_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL beta_numu_c011 = beta_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL beta_numu_c100 = beta_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL beta_numu_c101 = beta_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL beta_numu_c110 = beta_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL beta_numu_c111 = beta_numu_data[NMSParams_IDX(i1, j1, k1)];

    // Gamma
    BS_REAL gamma_numu_c000 = gamma_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL gamma_numu_c001 = gamma_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL gamma_numu_c010 = gamma_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL gamma_numu_c011 = gamma_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL gamma_numu_c100 = gamma_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL gamma_numu_c101 = gamma_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL gamma_numu_c110 = gamma_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL gamma_numu_c111 = gamma_numu_data[NMSParams_IDX(i1, j1, k1)];

    // Delta
    BS_REAL delta_numu_c000 = delta_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL delta_numu_c001 = delta_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL delta_numu_c010 = delta_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL delta_numu_c011 = delta_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL delta_numu_c100 = delta_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL delta_numu_c101 = delta_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL delta_numu_c110 = delta_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL delta_numu_c111 = delta_numu_data[NMSParams_IDX(i1, j1, k1)];

    // Wp_bar
    BS_REAL wpbar_numu_c000 = wpbar_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL wpbar_numu_c001 = wpbar_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL wpbar_numu_c010 = wpbar_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL wpbar_numu_c011 = wpbar_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL wpbar_numu_c100 = wpbar_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL wpbar_numu_c101 = wpbar_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL wpbar_numu_c110 = wpbar_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL wpbar_numu_c111 = wpbar_numu_data[NMSParams_IDX(i1, j1, k1)];

    // Wp_zero
    BS_REAL wpzero_numu_c000 = wpzero_numu_data[NMSParams_IDX(i0, j0, k0)];
    BS_REAL wpzero_numu_c001 = wpzero_numu_data[NMSParams_IDX(i0, j0, k1)];
    BS_REAL wpzero_numu_c010 = wpzero_numu_data[NMSParams_IDX(i0, j1, k0)];
    BS_REAL wpzero_numu_c011 = wpzero_numu_data[NMSParams_IDX(i0, j1, k1)];
    BS_REAL wpzero_numu_c100 = wpzero_numu_data[NMSParams_IDX(i1, j0, k0)];
    BS_REAL wpzero_numu_c101 = wpzero_numu_data[NMSParams_IDX(i1, j0, k1)];
    BS_REAL wpzero_numu_c110 = wpzero_numu_data[NMSParams_IDX(i1, j1, k0)];
    BS_REAL wpzero_numu_c111 = wpzero_numu_data[NMSParams_IDX(i1, j1, k1)];

    BS_REAL a0 = (one - tx), a1 = tx;
    BS_REAL b0 = (one - ty), b1 = ty;
    BS_REAL c0 = (one - tz), c1 = tz;

    // Alpha
    BS_REAL alpha_numu = alpha_numu_c000 * a0 * b0 * c0 + alpha_numu_c001 * a0 * b0 * c1 +
                        alpha_numu_c010 * a0 * b1 * c0 + alpha_numu_c011 * a0 * b1 * c1 +
                        alpha_numu_c100 * a1 * b0 * c0 + alpha_numu_c101 * a1 * b0 * c1 +
                        alpha_numu_c110 * a1 * b1 * c0 + alpha_numu_c111 * a1 * b1 * c1;

    // Beta
    BS_REAL beta_numu = beta_numu_c000 * a0 * b0 * c0 + beta_numu_c001 * a0 * b0 * c1 +
                        beta_numu_c010 * a0 * b1 * c0 + beta_numu_c011 * a0 * b1 * c1 +
                        beta_numu_c100 * a1 * b0 * c0 + beta_numu_c101 * a1 * b0 * c1 +
                        beta_numu_c110 * a1 * b1 * c0 + beta_numu_c111 * a1 * b1 * c1;

    // Gamma
    BS_REAL gamma_numu = gamma_numu_c000 * a0 * b0 * c0 + gamma_numu_c001 * a0 * b0 * c1 +
                        gamma_numu_c010 * a0 * b1 * c0 + gamma_numu_c011 * a0 * b1 * c1 +
                        gamma_numu_c100 * a1 * b0 * c0 + gamma_numu_c101 * a1 * b0 * c1 +
                        gamma_numu_c110 * a1 * b1 * c0 + gamma_numu_c111 * a1 * b1 * c1;

    // Delta
    BS_REAL delta_numu = delta_numu_c000 * a0 * b0 * c0 + delta_numu_c001 * a0 * b0 * c1 +
                        delta_numu_c010 * a0 * b1 * c0 + delta_numu_c011 * a0 * b1 * c1 +
                        delta_numu_c100 * a1 * b0 * c0 + delta_numu_c101 * a1 * b0 * c1 +
                        delta_numu_c110 * a1 * b1 * c0 + delta_numu_c111 * a1 * b1 * c1;

    // Wp_bar
    BS_REAL wpbar_numu = wpbar_numu_c000 * a0 * b0 * c0 + wpbar_numu_c001 * a0 * b0 * c1 +
                        wpbar_numu_c010 * a0 * b1 * c0 + wpbar_numu_c011 * a0 * b1 * c1 +
                        wpbar_numu_c100 * a1 * b0 * c0 + wpbar_numu_c101 * a1 * b0 * c1 +
                        wpbar_numu_c110 * a1 * b1 * c0 + wpbar_numu_c111 * a1 * b1 * c1;

    // Wp_zero
    BS_REAL wpzero_numu = wpzero_numu_c000 * a0 * b0 * c0 + wpzero_numu_c001 * a0 * b0 * c1 +
                        wpzero_numu_c010 * a0 * b1 * c0 + wpzero_numu_c011 * a0 * b1 * c1 +
                        wpzero_numu_c100 * a1 * b0 * c0 + wpzero_numu_c101 * a1 * b0 * c1 +
                        wpzero_numu_c110 * a1 * b1 * c0 + wpzero_numu_c111 * a1 * b1 * c1;

    parameters.alpha = alpha_numu;
    parameters.beta = beta_numu;
    parameters.gamma = gamma_numu;
    parameters.delta = delta_numu;
    parameters.wpbar = wpbar_numu;
    parameters.wpzero = wpzero_numu;

#undef NMSParams_IDX

    return parameters;
}


//=========================================================//
// --- Inelastic Neutrino Scattering off muons Kernel  --- //
// ---           via SEMI-ANALYTICAL KERNEL            --- //
//=========================================================//
/* Interpolate the parameters grid via 3D linear interpolation
 *
 * NOTE: this function is just for points inside the ranges of the 3D table!
 *
 *      Inputs:
 *          T:          temperature                [MeV]
 *          mu:         muon chemical potential    [MeV]
 *          w:          incoming neutrino energy   [MeV]
 *          wp:         outgoing neutrino energy   [MeV]
 * 
 *      Output:
 *           kernel nu_mu [MeV^-2]
 *           kernel anu_mu [MeV^-2]
 *           kernel nu_e [MeV^-2]
 *           kernel anu_e [MeV^-2]
 */

/*
struct NMS_Kernel_DiffFlavors{
    BS_REAL R[total_num_species];
};

KOKKOS_INLINE_FUNCTION
NMS_Kernel_DiffFlavors NMS_SemiAnalytical(const BS_REAL T, const BS_REAL mu,
                                  const BS_REAL w, const BS_REAL wp)
{

    NMS_Kernel_DiffFlavors kernels_diff_flavors = {};

    // Interpolation of the 4D table
    const NMS_Parameters parameters_values = NMS_parameters_interpolator(T, mu, w);

    // Neutrino-flavor dependent kernel
    //nu_e, nu_t
    kernels_diff_flavors.R[0] = NMS_semi_analytical_kernel(T, mu, w, wp, parameters_values.alpha, 
                                parameters_values.beta, parameters_values.gamma, parameters_values.delta,
                                parameters_values.wpbar, parameters_values.wpzero);
    
    //anu_e, anu_t
    kernels_diff_flavors.R[1] = NMS_semi_analytical_kernel(T, mu, w, wp, parameters_values.alpha, 
                                parameters_values.beta, parameters_values.gamma, parameters_values.delta,
                                parameters_values.wpbar, parameters_values.wpzero);
    
    //nu_mu==nu_x
    kernels_diff_flavors.R[2] = NMS_semi_analytical_kernel(T, mu, w, wp, parameters_values.alpha, 
                                parameters_values.beta, parameters_values.gamma, parameters_values.delta,
                                parameters_values.wpbar, parameters_values.wpzero);
    
    //anu_mu==anu_x
    kernels_diff_flavors.R[3] = NMS_semi_analytical_kernel(T, mu, w, wp, parameters_values.alpha, 
                                parameters_values.beta, parameters_values.gamma, parameters_values.delta,
                                parameters_values.wpbar, parameters_values.wpzero);
    
    return kernels_diff_flavors;
}
*/

// Calculates and saves the neutrino muon scattering in and out kernel
// @TODO: ADD PARAMETER INTERPOLATION FOR OTHER FLAVORS
// @TODO: GENERALIZE TO 6 NEUTRINO SPECIES (nue, nut are equal)
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticMuonScattKernels_SemiAnalytical(InelasticScattKernelParams* kernel_params,
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
    //nu_e, nu_t
    nms_kernel.abs[id_nue] = 0.0;
    
    //anu_e, anu_t
    nms_kernel.abs[id_anue] = 0.0;
    
    //nu_mu==nu_x
    nms_kernel.abs[id_nux] = NMS_semi_analytical_kernel(T, mu_mu, w, wp, parameters_values.alpha, 
                                parameters_values.beta, parameters_values.gamma, parameters_values.delta,
                                parameters_values.wpbar, parameters_values.wpzero);
    
    //anu_mu==anu_x
    nms_kernel.abs[id_anux] = 0.0;

    // Absorption and Emission kernels [nm^3 s^-1]
    // nu_e(idx=0), nu_ae(idx=1), nu_x=nu_mu(idx=2), anu_x=anu_mu(idx=3)
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        nms_kernel.em[idx]  = nms_kernel.abs[idx] * det_bal;
    }

    return nms_kernel;
}


#endif // BNS_NURATES_INCLUDE_KERNEL_NMS_SEMI_ANALYTICAL_HPP_
