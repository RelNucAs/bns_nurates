//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_brem_GP19.hpp
//  \brief contains bremsstrahlung kernels and associated helper functions
//
// Computation of nucleon-nucleon bremsstrahlung kernel using the approach of
// Gang Guo and Gabriel Martínez-Pinedo 2019 ApJ 887 58

#ifndef BNS_NURATES_INCLUDE_KERNEL_BREM_GP19_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_BREM_GP19_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"

#include "GP19Table.hpp"

//===================================//
// --- 4D interpolator structure --- //
//===================================//
/* Finds the idx such that the chosen value is bracketed
 *
 *      Inputs:
 *           value:     value to interpolate    [any]
 *           axis:      1D axis of the table    [any]
 *           size:      size of the 1D axis     [/]
 *           i0:        index 0                 [/]
 *           i1:        index 1                 [/]
 *           t:         relative distance       [/]
 */
KOKKOS_INLINE_FUNCTION
int gp19_find_bracketing_indices(const BS_REAL value, const BS_REAL* axis,
                                 const int size, int* i0, int* i1, BS_REAL* t)
{
    // Out of range
    if (value < axis[0] || value > axis[size - 1])
    {
        return -1;
    }

    for (int i = 0; i < size - 1; ++i)
    {
        if (value >= axis[i] and value <= axis[i + 1])
        {
            *i0 = i;
            *i1 = i + 1;
            *t  = (value - axis[i]) / (axis[i + 1] - axis[i]);

            return 0;
        }
    }

    // Note found
    return -1;
}

/* 4D linear interpolation of the numerical table of the response function
 *
 *      Inputs:
 *          nb:         baryon number density       [fm-3]
 *          Ye:         electron fraction           [/]
 *          T :         temperature                 [MeV]
 *          w :         pair neutrino energy        [MeV]
 *          *table:     GP19Table-type struct
 *
 *      Output:
 *          c:          Interpolated value          [MeV-1]
 */
KOKKOS_INLINE_FUNCTION
BS_REAL gp19_interpolator(const BS_REAL nb, const BS_REAL Ye, const BS_REAL T,
                          const BS_REAL w)
{
    int i0, i1, j0, j1, k0, k1, l0, l1;
    BS_REAL tx, ty, tz, tw;

    // log-linear for nb e w
    BS_REAL log_nb = log(nb);
    BS_REAL log_w  = log(w);

    if (gp19_find_bracketing_indices(log_nb, GP19_lognb_axis, GP19_nb_dims, &i0,
                                     &i1, &tx) < 0 or
        gp19_find_bracketing_indices(Ye, GP19_Ye_axis, GP19_Ye_dims, &j0, &j1,
                                     &ty) < 0 or
        gp19_find_bracketing_indices(T, GP19_T_axis, GP19_T_dims, &k0, &k1,
                                     &tz) < 0 or
        gp19_find_bracketing_indices(log_w, GP19_logw_axis, GP19_w_dims, &l0,
                                     &l1, &tw) < 0)
    {
        return std::numeric_limits<BS_REAL>::quiet_NaN();
    }

// macro to access the 1D array
#define IDX(i, j, k, l)                                                        \
    (((i * GP19_Ye_dims + j) * GP19_T_dims + k) * GP19_w_dims + l)

    // extract the 16 vertexes of the hypercube
    BS_REAL c0000 = GP19_data[IDX(i0, j0, k0, l0)];
    BS_REAL c0001 = GP19_data[IDX(i0, j0, k0, l1)];
    BS_REAL c0010 = GP19_data[IDX(i0, j0, k1, l0)];
    BS_REAL c0011 = GP19_data[IDX(i0, j0, k1, l1)];
    BS_REAL c0100 = GP19_data[IDX(i0, j1, k0, l0)];
    BS_REAL c0101 = GP19_data[IDX(i0, j1, k0, l1)];
    BS_REAL c0110 = GP19_data[IDX(i0, j1, k1, l0)];
    BS_REAL c0111 = GP19_data[IDX(i0, j1, k1, l1)];
    BS_REAL c1000 = GP19_data[IDX(i1, j0, k0, l0)];
    BS_REAL c1001 = GP19_data[IDX(i1, j0, k0, l1)];
    BS_REAL c1010 = GP19_data[IDX(i1, j0, k1, l0)];
    BS_REAL c1011 = GP19_data[IDX(i1, j0, k1, l1)];
    BS_REAL c1100 = GP19_data[IDX(i1, j1, k0, l0)];
    BS_REAL c1101 = GP19_data[IDX(i1, j1, k0, l1)];
    BS_REAL c1110 = GP19_data[IDX(i1, j1, k1, l0)];
    BS_REAL c1111 = GP19_data[IDX(i1, j1, k1, l1)];

    BS_REAL a0 = (1 - tx), a1 = tx;
    BS_REAL b0 = (1 - ty), b1 = ty;
    BS_REAL c0 = (1 - tz), c1 = tz;
    BS_REAL d0 = (1 - tw), d1 = tw;

    BS_REAL c = c0000 * a0 * b0 * c0 * d0 + c0001 * a0 * b0 * c0 * d1 +
                c0010 * a0 * b0 * c1 * d0 + c0011 * a0 * b0 * c1 * d1 +
                c0100 * a0 * b1 * c0 * d0 + c0101 * a0 * b1 * c0 * d1 +
                c0110 * a0 * b1 * c1 * d0 + c0111 * a0 * b1 * c1 * d1 +
                c1000 * a1 * b0 * c0 * d0 + c1001 * a1 * b0 * c0 * d1 +
                c1010 * a1 * b0 * c1 * d0 + c1011 * a1 * b0 * c1 * d1 +
                c1100 * a1 * b1 * c0 * d0 + c1101 * a1 * b1 * c0 * d1 +
                c1110 * a1 * b1 * c1 * d0 + c1111 * a1 * b1 * c1 * d1;

#undef IDX

    return c;
}

//=============================================//
// --- Kernel GP19 inside the table ranges --- //
//=============================================//
/* Compute the Bremsstrahlung absorption reaction kernel of Guo-Martinez Pinedo
 * (2019) via 4D linear interpolation
 *
 * NOTE: this function is just for points inside the ranges of the 4D table!
 *
 *      Inputs:
 *          nb:         baryon number density   [fm-3]
 *          Ye:         electron fraction       [/]
 *          T:          temperature             [MeV]
 *          w1:         neutrino 1 energy       [MeV]
 *          w2:         neutrino 2 energy       [MeV]
 *
 *      Output:
 *           kernel:    Bremsstrahlung reaction kernel      [nm3 s-1]
 */
KOKKOS_INLINE_FUNCTION
BS_REAL gp19_BremKernel_in_ranges(const BS_REAL nb, const BS_REAL Ye,
                                  const BS_REAL T, const BS_REAL w)
{

    // Interpolation of the 4D table
    const BS_REAL S_value = gp19_interpolator(nb, Ye, T, w);

    const BS_REAL kernel =
        kBS_Brem_Const * (nb * 1e+18) * S_value; // [nm^3 s^-1]

    return kernel;
}


//============================================================//
// Definition of the different kernel extrapolating functions //
//============================================================//
KOKKOS_INLINE_FUNCTION
BS_REAL gp19_LowDensityExtrapolation(const BS_REAL nb_below, const BS_REAL Ye,
                                     const BS_REAL T, const BS_REAL w)
{
    constexpr BS_REAL nb_min = GP19_nb_axis[0]; // [fm^-3]

    const BS_REAL kernel_at_nb_min =
        gp19_BremKernel_in_ranges(nb_min, Ye, T, w); // [nm3 s-1]

    // quadratic extrapolation
    const BS_REAL kernel_extra = kernel_at_nb_min * POW2(nb_below / nb_min);

    return kernel_extra;
}


KOKKOS_INLINE_FUNCTION
BS_REAL gp19_LowTemperatureExtrapolation(const BS_REAL nb, const BS_REAL Ye,
                                         const BS_REAL T_below, const BS_REAL w)
{

    constexpr BS_REAL Thalf    = 0.5;                // [MeV]
    constexpr BS_REAL T2       = 2.0;                // [MeV]
    constexpr BS_REAL T4       = 4.0;                // [MeV]
    constexpr BS_REAL ilogT4T2 = 1.4426950408889634; // 1 / log(4 / 2)

    const BS_REAL kernel_T2 =
        gp19_BremKernel_in_ranges(nb, Ye, T2, w); // [nm3 s-1]
    const BS_REAL kernel_T4 =
        gp19_BremKernel_in_ranges(nb, Ye, T4, w); // [nm3 s-1]

    const BS_REAL m = log(kernel_T4 / kernel_T2) * ilogT4T2;

    // Power law
    const BS_REAL kernel_extra = kernel_T2 * pow(T_below * Thalf, m);

    return kernel_extra;
}


KOKKOS_INLINE_FUNCTION
BS_REAL gp19_LowDensity_and_LowTemperatureExtrapolation(const BS_REAL nb_below,
                                                        const BS_REAL Ye,
                                                        const BS_REAL T_below,
                                                        const BS_REAL w)
{
    constexpr BS_REAL nb_min   = GP19_nb_axis[0];    // [fm-3]
    constexpr BS_REAL Thalf    = 0.5;                // [MeV]
    constexpr BS_REAL T2       = 2.0;                // [MeV]
    constexpr BS_REAL T4       = 4.0;                // [MeV]
    constexpr BS_REAL ilogT4T2 = 1.4426950408889634; // 1 / log(4 / 2)

    BS_REAL kernel_LowDensity_T2 =
        gp19_LowDensityExtrapolation(nb_below, Ye, T2, w); // [nm3 s-1]
    BS_REAL kernel_LowDensity_T4 =
        gp19_LowDensityExtrapolation(nb_below, Ye, T4, w); // [nm3 s-1]

    // Small floor value to avoid log(0)
    constexpr BS_REAL eps = 1e-30;
    if (kernel_LowDensity_T2 < eps)
    {
        kernel_LowDensity_T2 = eps;
    }
    if (kernel_LowDensity_T4 < eps)
    {
        kernel_LowDensity_T4 = eps;
    }

    const BS_REAL m =
        log(kernel_LowDensity_T4 / kernel_LowDensity_T2) * ilogT4T2;

    BS_REAL kernel_extra = kernel_LowDensity_T2 * pow(T_below * Thalf, m);

    // Guarantees a minimum kernel
    if (kernel_extra < eps)
    {
        kernel_extra = eps;
    }

    return kernel_extra;
}


KOKKOS_INLINE_FUNCTION
BS_REAL gp19_HighTemperatureExtrapolation(const BS_REAL nb, const BS_REAL Ye,
                                          const BS_REAL T_above,
                                          const BS_REAL w)
{
    constexpr BS_REAL one          = 1.0;
    constexpr BS_REAL two          = 2.0;
    constexpr BS_REAL four         = 4.0;
    constexpr BS_REAL eight        = 8.0;
    const BS_REAL nb_nm            = nb * 1e+18;         // [nm^3]
    const BS_REAL nb_nm_squared    = POW2(nb_nm);        // [nm^6]
    const BS_REAL omegaPairSquared = POW2(w);            // [MeV^2]
    constexpr BS_REAL Temp48       = 48.0;               // [MeV]
    constexpr BS_REAL Temp50       = 50.0;               // [MeV]
    constexpr BS_REAL sqrtTemp50   = 7.0710678118654755; // [MeV^1/2]

    // Kernel @ T = 48, 50 MeV
    const BS_REAL R48 =
        gp19_BremKernel_in_ranges(nb, Ye, Temp48, w); // [nm^3 s^-1]
    const BS_REAL R50 =
        gp19_BremKernel_in_ranges(nb, Ye, Temp50, w); // [nm^3 s^-1]

    const BS_REAL K = four * kBS_Brem_Const * nb_nm_squared; // [MeV s^-1]
    const BS_REAL A = one / R50;                             // [nm^3 s]
    const BS_REAL Asquared = POW2(A);
    const BS_REAL B        = one / (K * sqrtTemp50); // [MeV^-3/2 s]
    const BS_REAL Bsquared = POW2(B);

    const BS_REAL m = (R50 - R48) / (Temp50 - Temp48); // [nm^3 s^-1 MeV^-1]

    // Beta and Gamma params
    const BS_REAL beta_num  = eight * omegaPairSquared * K * Bsquared;
    const BS_REAL beta_den1 = two * m * Asquared * sqrtTemp50;
    const BS_REAL beta_den2 = A * K * B;
    const BS_REAL beta = beta_num / (beta_den1 + beta_den2); // [nm^3 MeV^1/2]

    const BS_REAL gamma_num = (A * beta) - (four * B * omegaPairSquared);
    const BS_REAL gamma_den = B * nb_nm_squared * Temp50;
    const BS_REAL gamma     = gamma_num / gamma_den; // [nm^6 MeV]

    // Lorentzian function for a general T > 50 MeV
    const BS_REAL Lor_den1 = (w / T_above) * (w / T_above);
    const BS_REAL Lor_den2 = (nb_nm_squared / (four * T_above)) * gamma;
    const BS_REAL inv_den_Lorentzian = one / (Lor_den1 + Lor_den2);

    const BS_REAL Lor_scaling = (beta * nb_nm) / (T_above * sqrt(T_above));

    const BS_REAL S_Lorentzian = Lor_scaling * inv_den_Lorentzian; // [MeV^-1]

    const BS_REAL kernel_extra =
        kBS_Brem_Const * nb_nm * S_Lorentzian; // [nm^3 s^-1]

    return kernel_extra;
}


KOKKOS_INLINE_FUNCTION
BS_REAL gp19_HighTemperature_and_LowDensityExtrapolation(BS_REAL nb_below,
                                                         const BS_REAL Ye,
                                                         const BS_REAL T_above,
                                                         const BS_REAL w)
{
    // NOTE: verified that extrapolating in T first and in density second the
    // result changes at the 16th figure

    constexpr BS_REAL nb_min = GP19_nb_axis[0]; // [fm-3]

    const BS_REAL kernel_at_nb_min = gp19_HighTemperatureExtrapolation(
        nb_min, Ye, T_above, w); // [nm^3 s^-1]
                                 //
    const BS_REAL kernel_extra =
        kernel_at_nb_min *
        ((nb_below / nb_min) * (nb_below / nb_min)); // [nm^3 s^-1]

    return kernel_extra;
}


//==============================================//
// --- Kernel GP19 over the whole 4D domain --- //
//==============================================//
/* Compute the Bremsstrahlung absorption reaction kernel of Guo-Martinez Pinedo
 * (2019) over the full 4D domain
 *
 *      Inputs:
 *          *eos:               MyEOSParams-type struct
 *              nb_nm:          baryon number density           [nm^-3]
 *              Ye:             electron fraction               [/]
 *              T:              temperature                     [MeV]
 *          *bremParams:        BremKernelParams-type struct
 *              omega:          neutrino 1 energy               [MeV]
 *              omega_prime:    neutrino 2 energy               [MeV]
 *
 *      Output:
 *           kernel:    Bremsstrahlung reaction kernel      [nm^3 s-1]
 */
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelAbsGP19(const BremKernelParams* bremParams,
                                 const MyEOSParams* eos)
{
    constexpr BS_REAL three = 3;
    constexpr BS_REAL nm2fm_3 = 1e-18;

    BS_REAL nb_nm = eos->nb;                                     // [nm^-3]
    BS_REAL Ye    = eos->ye;                                     // [/]
    BS_REAL T     = eos->temp;                                   // [MeV]
    BS_REAL w     = bremParams->omega + bremParams->omega_prime; // [MeV]

    BS_REAL nb = nb_nm * nm2fm_3; // [fm^-3]

    // Min & Max definition
    const BS_REAL nb_max = GP19_nb_axis[GP19_nb_dims - 1]; // [1     fm^-3]
    const BS_REAL nb_min = GP19_nb_axis[0];                // [1e-4  fm^-3]
    const BS_REAL Ye_max = GP19_Ye_axis[GP19_Ye_dims - 1]; // [0.5]
    const BS_REAL T_max  = GP19_T_axis[GP19_T_dims - 1];   // [50.0  MeV]
    const BS_REAL T_min  = GP19_T_axis[0];                 // [2.0   MeV]
    const BS_REAL w_max  = GP19_w_axis[GP19_w_dims - 1];   // [~398  MeV]
    const BS_REAL w_min  = GP19_w_axis[0];                 // [~0.05 MeV]

    //==================================//
    // --- Manual shifts / clipping --- //
    //==================================//

    // High Ye
    if (Ye > Ye_max)
    {
        Ye = Ye_max;
    }

    // High nb
    if (nb > nb_max)
    {
        nb = nb_max;
    }

    // High w
    if (w > w_max)
    {
        w = w_max;
    }

    // Low w
    if (w < w_min)
    {
        w = w_min;
    }

    BS_REAL s_abs;

    //================================//
    // --- Standard interpolation --- //
    //================================//
    if (nb >= nb_min && T >= T_min && T <= T_max)
    {
        s_abs = gp19_BremKernel_in_ranges(nb, Ye, T, w);
    }

    //==========================================//
    // --- Extrapolation in special domains --- //
    //==========================================//

    // Case 1: low-density (< 1e-4 fm-3)
    if (nb < nb_min)
    {
        // Low-temperature (< 2 MeV)
        if (T < T_min)
        {
            s_abs =
                gp19_LowDensity_and_LowTemperatureExtrapolation(nb, Ye, T, w);
        }
        // High-temperature (> 50 MeV)
        else if (T > T_max)
        {
            s_abs =
                gp19_HighTemperature_and_LowDensityExtrapolation(nb, Ye, T, w);
        }
        // Tempearature in range
        else
        {
            s_abs = gp19_LowDensityExtrapolation(nb, Ye, T, w);
        }
    }
    // Case 2: nb in range but T out of range
    // Low-temperature (< 2 MeV)
    else if (T < T_min)
    {
        s_abs = gp19_LowTemperatureExtrapolation(nb, Ye, T, w);
    }
    // High-temperature (> 50 MeV)
    else if (T > T_max)
    {
        s_abs = gp19_HighTemperatureExtrapolation(nb, Ye, T, w);
    }

    // The factor of three comes from the fact that this is the 0th order term
    // of the Legendre expansion, e.g. the angular-independent one. Confront
    // with the analogue factors present in the HR98 verion of the
    // Bremsstrahlung kernel.
    s_abs *= three;

    // Production kernel from detailed balance
    BS_REAL s_em = s_abs * SafeExp(-w / T);

    MyKernelOutput brem_kernel;
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_BREM_GP19_HPP_
