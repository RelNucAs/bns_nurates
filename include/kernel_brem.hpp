//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_brem.c
//  \brief contains bremsstrahlung kernels and associated helper functions
//
// Computation of nucleon-nucleon bremsstrahlung kernel using the analytic
// fitting formula in Hannestad & Raffelt 1998, Apj, 507, 339
// (https://iopscience.iop.org/article/10.1086/306303/pdf)

#ifndef BNS_NURATES_INCLUDE_KERNEL_BREM_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_BREM_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"


/* Compute the analytical fit for the s-component of the kernel for
 * neutrino bremsstrahlung and inelastic scattering in a nucleon field
 *
 * Note: Does not support negative x!
 *
 * Inputs:
 *      x:        rescaled total neutrino energy (w+wp/T)
 *      y:        pion mass parameter defined in Eqn. (38)
 *      eta_star: nucleon degeneracy parameter
 *
 * Output:
 *      s:  a dimensionless quantity as defined in Eqn. (49)
 */
KOKKOS_INLINE_FUNCTION
BS_REAL BremKernelS(BS_REAL x, BS_REAL y, BS_REAL eta_star)
{
    constexpr BS_REAL zero   = 0;
    constexpr BS_REAL one    = 1;
    constexpr BS_REAL two    = 2;
    constexpr BS_REAL three  = 3;
    constexpr BS_REAL four   = 4;
    constexpr BS_REAL five   = 5;
    constexpr BS_REAL six    = 6;
    constexpr BS_REAL twelve = 12;
    constexpr BS_REAL thirty = 30;

    constexpr BS_REAL three_halves         = 1.5;
    constexpr BS_REAL five_halves          = 2.5;
    constexpr BS_REAL eleven_tenth         = 1.1;
    constexpr BS_REAL twelve_tenth         = 1.2;
    constexpr BS_REAL one_fifth            = 0.2;
    constexpr BS_REAL two_fifth            = 0.4;
    constexpr BS_REAL four_fifth           = 0.8;
    constexpr BS_REAL fourteen_fifth       = 2.8;
    constexpr BS_REAL five_thousands       = 0.005;
    constexpr BS_REAL six_hundreds         = 0.06;
    constexpr BS_REAL one_tenth            = 0.1;
    constexpr BS_REAL twentythree_tenth    = 2.3;
    constexpr BS_REAL ninetythree_hundreds = 0.93;
    constexpr BS_REAL sixtyseven_hundreds  = 0.67;
    constexpr BS_REAL eighteen_hundreds    = 0.18;

    constexpr BS_REAL ten_minus_fourteen = 1e-14;
    constexpr BS_REAL ten_minus_ten      = 1.e-10;

    constexpr BS_REAL c1 = 0.0044;
    constexpr BS_REAL c2 = 1.05;
    constexpr BS_REAL c3 = 0.0001;
    constexpr BS_REAL c4 = 2.39;

    BS_ASSERT(x >= zero);
    BS_ASSERT(y >= zero);
    BS_ASSERT(isfinite(eta_star));
    BS_ASSERT(eta_star >= zero);

    // Prevent singular behavior
    x        = (x > kBS_Brem_Xmin) ? x : kBS_Brem_Xmin;
    y        = (y > kBS_Brem_Ymin) ? y : kBS_Brem_Ymin;
    eta_star = (eta_star > kBS_Brem_Etamin) ? eta_star : kBS_Brem_Etamin;

    // Compute non-degenerate approximation, s_nd in Eqn. (45)
    const BS_REAL s_nd_numerator =
        two * kBS_SqrtPi * pow(x + two - SafeExp(-y / twelve), three_halves) *
        (POW2(x) + two * x * y + five * POW2(y) / three + one);
    const BS_REAL s_nd_denominator =
        kBS_SqrtPi + POW4(kBS_Pi2OneEighth + x + y);
    const BS_REAL s_nd = s_nd_numerator / s_nd_denominator;

    // Compute degenerate approximation, s_d in Eqn. (46)
    const BS_REAL u     = sqrt(y / (two * eta_star)) + ten_minus_ten;
    const BS_REAL u2    = POW2(u);
    const BS_REAL u_arg = u2 / (two * sqrt(two * u2 + four));
    BS_REAL f_u =
        (one - five * u * atan(two / u) / six + u2 / (three * (u2 + four)) +
         atan(one / u_arg) * u_arg / three);

    // @TODO: Leonardo check this! Doing this to prevent s_d from being a large
    // negative number
    f_u = (fabs(f_u) < ten_minus_fourteen) ? ten_minus_fourteen : f_u;

    // f_u = fabs(f_u);

    const BS_REAL s_d = three * kBS_PiHalfToFiveHalves *
                        pow(eta_star, -five_halves) *
                        (POW2(x) + kBS_FourPiSquared) * x * f_u /
                        (kBS_FourPiSquared * (one - SafeExp(-x)));

    const BS_REAL pow_x_1_1 = pow(x, eleven_tenth);

    // F, Eqn. (50)
    const BS_REAL f_denominator =
        (three + POW2(x - twelve_tenth) + pow(x, -four)) *
        (one + POW2(eta_star)) * (one + POW4(y));
    // Eq.(50)
    const BS_REAL f_brem = one + one / f_denominator;

    // G, Eqn. (50)
    const BS_REAL g_brem = one - c1 * pow_x_1_1 * y /
                                     (four_fifth + six_hundreds * pow(y, c2)) *
                                     sqrt(eta_star) / (eta_star + one_fifth);

    // h and C, Eqn. (50)
    const BS_REAL h_brem =
        one_tenth * eta_star / (c4 + one_tenth * pow(eta_star, eleven_tenth));
    const BS_REAL c_brem =
        eleven_tenth * pow_x_1_1 * h_brem /
        (twentythree_tenth + h_brem * pow(x, ninetythree_hundreds) +
         c3 * pow(x, twelve_tenth)) *
        thirty / (thirty + five_thousands * pow(x, fourteen_fifth));

    // p, Eqn. (50)
    const BS_REAL p_brem =
        sixtyseven_hundreds + eighteen_hundreds * pow(y, two_fifth);

    // interpolated formula for s in Eqn. (49)
    const BS_REAL s_brem =
        pow(pow(s_nd, -p_brem) + pow(s_d, -p_brem), -one / p_brem) * f_brem *
        (one + c_brem * g_brem);

    BS_ASSERT(s_brem >= zero);

    return s_brem;
}

/* Compute the analytic fit of the g-component of the kernel for neutrino
 * bremsstrahlung and inelastic scattering in a nucleon field. Implements
 * Eqn. (52) of Hannestad & Raffelt (1998).
 *
 * Inputs:
 *    y:        pion mass parameter defined in Eqn. (38)
 *    eta_star: nucleon degeneracy parameter
 *
 * Output:
 *    g: a dimensionless quantity as defined in Eqn. (52)
 */
KOKKOS_INLINE_FUNCTION
BS_REAL BremKernelG(BS_REAL y, BS_REAL eta_star)
{
    constexpr BS_REAL zero       = 0;
    constexpr BS_REAL one        = 1;
    constexpr BS_REAL two        = 2;
    constexpr BS_REAL half       = 0.5;
    constexpr BS_REAL twentyfive = 25;

    constexpr BS_REAL three_halves       = 1.5;
    constexpr BS_REAL five_halves        = 2.5;
    constexpr BS_REAL two_fifth          = 0.4;
    constexpr BS_REAL two_hundreds       = 0.02;
    constexpr BS_REAL four_hundreds      = 0.04;
    constexpr BS_REAL five_hundreds      = 0.05;
    constexpr BS_REAL twelve_tenth       = 1.2;
    constexpr BS_REAL six_tenth          = 0.6;
    constexpr BS_REAL eighteen_tenth     = 1.8;
    constexpr BS_REAL twentythree_tenth  = 2.3;
    constexpr BS_REAL fortyfive_hundreds = 0.45;
    constexpr BS_REAL fifteen_hundreds   = 0.15;

    constexpr BS_REAL c1 = 15.6;
    constexpr BS_REAL c2 = 0.63;
    constexpr BS_REAL c3 = 1.45;
    constexpr BS_REAL c4 = 0.025;
    constexpr BS_REAL c5 = 13.75;

    BS_ASSERT(y >= zero);
    BS_ASSERT(eta_star >= zero);

    // prevent singular behavior
    y        = (y > kBS_Brem_Ymin) ? y : kBS_Brem_Ymin;
    eta_star = (eta_star > kBS_Brem_Etamin) ? eta_star : kBS_Brem_Etamin;

    // alpha_1, Eqn. (53)
    const BS_REAL y2            = POW2(y);
    const BS_REAL eta_star_inv  = one / eta_star;
    const BS_REAL alpha_1_denom = twentyfive * y2 + one;

    const BS_REAL alpha_1 =
        (half + eta_star_inv) / (one + eta_star_inv) * (one / alpha_1_denom) +
        (half + eta_star / c1) * twentyfive * y2 / alpha_1_denom;

    // alpha_2, Eqn. (53)
    const BS_REAL alpha_2 = (c2 + four_hundreds * pow(eta_star, c3)) /
                            (one + two_hundreds * pow(eta_star, five_halves));

    const BS_REAL pow_eta_star_1_5 = pow(eta_star, three_halves);

    // alpha_3, Eqn. (53)
    const BS_REAL alpha_3 =
        twelve_tenth *
        SafeExp(six_tenth * eta_star - two_fifth * pow_eta_star_1_5);

    // p_1, Eqn. (53)
    const BS_REAL p_1 = (eighteen_tenth + fortyfive_hundreds * eta_star) /
                        (one + fifteen_hundreds * pow_eta_star_1_5);

    // p_2, Eqn. (53)
    const BS_REAL p_2 =
        twentythree_tenth - five_hundreds * eta_star / (one + c4 * eta_star);

    // g, Eqn. (52)
    const BS_REAL g =
        (alpha_1 + alpha_2 * pow(y, p_1)) /
        (one + alpha_3 * pow(y, p_2) + alpha_2 * pow(y, p_1 + two) / c5);

    BS_ASSERT(g >= zero);

    return g;
}

/* Compute the absorption kernels for a given NN Bremsstrahlung channel */
KOKKOS_INLINE_FUNCTION
BS_REAL BremSingleChannelAbsKernel(const BS_REAL n_nuc, const BS_REAL m_nuc,
                                   BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
    constexpr BS_REAL two          = 2;
    constexpr BS_REAL three        = 3;
    constexpr BS_REAL half         = 0.5;
    constexpr BS_REAL one_tenth    = 0.1;
    constexpr BS_REAL two_thirds   = 2. / 3.;
    constexpr BS_REAL three_halves = 1.5;

    constexpr BS_REAL c1 = 1.63;
    constexpr BS_REAL c2 = 1.94;

    // EOS parameters
    // Temperature
    const BS_REAL T = eos_params->temp; // [MeV]

    // kernel parameters
    // Neutrino energy
    const BS_REAL omega = kernel_params->omega; // [MeV]
    // Primed neutrino energy
    const BS_REAL omega_prime = kernel_params->omega_prime; //  [MeV]

    // Dimensionless neutrino energy sum
    const BS_REAL x = (omega + omega_prime) / T;

    // Temperature in units of 10 MeV
    const BS_REAL T_10 = T * one_tenth;

    // Nucleon effective degeneracy parameter, Eqn. (36) using baryon number
    // density instead of matter density
    const BS_REAL eta_star = pow(three * kBS_PiSquared * n_nuc, two_thirds) *
                             kBS_Brem_Aux1 / (two * m_nuc * T);

    // Spin-fluctuation rate gamma, Eqn. (37)
    const BS_REAL gamma = c1 * pow(eta_star, three_halves) * T_10;

    // Pion mass parameter y, Eqn. (38)
    const BS_REAL y = c2 / T_10;

    // Dimensionless fitting parameter s
    const BS_REAL sb = BremKernelS(x, y, eta_star);

    // Dimensionless fitting parameter g
    const BS_REAL gb = BremKernelG(y, eta_star);

    // Differential absorption kernel, Eqn. (35)
    return gamma / (POW2(x) + POW2(half * gamma * gb)) * sb / T;
}

/* Compute the angular independent part of the absorption kernels for the
 Bremsstrahlung reactions by summing the contributions of all NN channels */
KOKKOS_INLINE_FUNCTION
BS_REAL BremAllChannelsAbsKernel(BremKernelParams* kernel_params,
                                 MyEOSParams* eos_params)
{
    constexpr BS_REAL twentyeight_thirds = 28. / 3.;
    constexpr BS_REAL one                = 1;
    constexpr BS_REAL three              = 3;

    // EOS parameters
    const BS_REAL nb = eos_params->nb; // baryon number density [nm^-3]
    const BS_REAL xn = eos_params->yn; // neutron abundance/mass fraction
    const BS_REAL xp = eos_params->yp; // proton abundance/mass fraction

    const BS_REAL x_mean =
        sqrt(xn * xp); // geometric mean of nucleon abundances/mass fractions

    const BS_REAL nn = nb * xn; // neutron number density [nm^-3]
    const BS_REAL np = nb * xp; // protron number density [nm^-3]
    const BS_REAL n_mean =
        nb * x_mean; // geometric mean of nucleon number densities [nm^-3]

    // compute single channel kernels
    // Neutron-neutron
    const BS_REAL s_abs_nn =
        BremSingleChannelAbsKernel(nn, kBS_MnGrams, kernel_params, eos_params);
    // Neutron-neutron
    const BS_REAL s_abs_pp =
        BremSingleChannelAbsKernel(np, kBS_MpGrams, kernel_params, eos_params);
    // Neutron-proton
    const BS_REAL s_abs_np = BremSingleChannelAbsKernel(
        n_mean, kBS_MAvgGrams, kernel_params, eos_params);

    // Total absorption kernel
    BS_REAL s_abs_tot =
        kBS_Brem_Const * (nn * s_abs_nn + np * s_abs_pp +
                          twentyeight_thirds * n_mean * s_abs_np);

    // kernel correction due to medium dependence as in Fischer2016
    if (kernel_params->use_NN_medium_corr == true)
    {
        s_abs_tot =
            s_abs_tot / POW6((one + cbrt(nb / kBS_Saturation_n) / three));
    }

    return s_abs_tot;
}

/* Compute a specific Legendre coefficient in the expansion of production and
 * absorption kernels for the Bremsstrahlung reactions */
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsLegCoeff(BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
    constexpr BS_REAL one   = 1;
    constexpr BS_REAL three = 3;

    // kernel parameters
    const int l         = kernel_params->l;     // order of Legendre coefficient
    const BS_REAL omega = kernel_params->omega; // neutrino energy [MeV]
    const BS_REAL omega_prime =
        kernel_params->omega_prime; // primed neutrino energy [MeV]

    BS_ASSERT(l >= 0 && l <= 1);

    // EOS parameters
    const BS_REAL temp = eos_params->temp; // temperature [MeV]

    // dimensionless neutrino energy sum
    const BS_REAL x = (omega + omega_prime) / temp;

    // angular independent part of absorption kernel
    BS_REAL s_abs = BremAllChannelsAbsKernel(kernel_params, eos_params);

    switch (l)
    {
    case 0:
        s_abs = three * s_abs; // zeroth Legedre coefficient
        break;
    case 1:
        s_abs = -one * s_abs; // first Legedre coefficient
        break;
    default:
        printf("BremKernelsLegCoeff (kernel_brem.c): l = %d must be either 0 "
               "or 1\n",
               l);
    }

    // production kernel from detailed balance
    BS_REAL s_em = s_abs * SafeExp(-x);

    MyKernelOutput brem_kernel;

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}

KOKKOS_INLINE_FUNCTION
void BremKernelsTable(const int n, BS_REAL* nu_array,
                      GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput brem_ker;

    grey_pars->kernel_pars.brem_kernel_params.l = 0;
    grey_pars->kernel_pars.brem_kernel_params.use_NN_medium_corr =
        grey_pars->opacity_pars.use_NN_medium_corr;

    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            // compute the brem kernels
            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[j];

            brem_ker =
                BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params,
                                    &grey_pars->eos_pars);

            out->m1_mat_em[0][i][j] = brem_ker.em[0];
            out->m1_mat_em[0][j][i] = brem_ker.em[0];

            out->m1_mat_ab[0][i][j] = brem_ker.abs[0];
            out->m1_mat_ab[0][j][i] = brem_ker.abs[0];
        }
    }

    return;
}


/* NN bremsstrahlung rates from BRT06 */


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

KOKKOS_INLINE_FUNCTION
void BremKernelsTableBRT06(const int n, BS_REAL* nu_array,
                           GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput brem_ker;

    for (int i = 0; i < n; ++i)
    {

        for (int j = i; j < n; ++j)
        {

            // compute the brem kernels
            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[j];

            brem_ker =
                BremKernelsBRT06(&grey_pars->kernel_pars.brem_kernel_params,
                                 &grey_pars->eos_pars);

            out->m1_mat_em[0][i][j] = brem_ker.em[0];
            out->m1_mat_em[0][j][i] = brem_ker.em[0];

            out->m1_mat_ab[0][i][j] = brem_ker.abs[0];
            out->m1_mat_ab[0][j][i] = brem_ker.abs[0];
        }
    }

    return;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_BREM_HPP_
