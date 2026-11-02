//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_brem_HR98.hpp
//  \brief contains bremsstrahlung kernels and associated helper functions
//
// Computation of nucleon-nucleon bremsstrahlung kernel using the analytic
// fitting formula in Hannestad & Raffelt 1998, Apj, 507, 339
// (https://iopscience.iop.org/article/10.1086/306303/pdf)

#ifndef BNS_NURATES_INCLUDE_KERNEL_BREM_HR98_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_BREM_HR98_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"


// ------------------------------------------------------------------------ //
// ------------------------- OPTIMIZED HR98 CODE -------------------------- //
// ------------------------------------------------------------------------ //

/* The HR98 kernel is split in two stages, so that the same code is used for
 * spectral (1D integral over nu_bar) and grey (2D integral over nu, nu_bar)
 * rates:
 *
 *   1. BremHR98PrecomputeEOS: everything that depends only on the EOS
 *      (y, eta_star of each NN channel, f(u), the whole g-component, gamma,
 *      F/G/h factors in eta_star, Fischer16 medium correction). To be called
 *      once per EOS point, outside of the loops over neutrino energies.
 *   2. BremKernelsLegCoeffFromEOS: the part that depends on
 *      x = (omega + omega_prime) / T. To be called for each energy pair.
 *
 * BremKernelsLegCoeff = 1. + 2., for standalone use.
 *
 * Non-integer powers are derived from a single log per variable
 * (e.g. x^1.2 = x (x^0.1)^2, y^0.4 = (y^0.05)^8, eta^1.45 = eta (eta^0.05)^9),
 * powers with integer or half-integer exponents from products and sqrt.
 */


/* Compute the EOS-dependent quantities of the HR98 kernel for all the NN
 * channels.
 *
 * Inputs:
 *      eos_params:         EOS parameters (nb, yn, yp, temp)
 *      use_NN_medium_corr: include the medium correction as in Fischer16
 *
 * Output:
 *      BremHR98EOSTerms struct
 */
KOKKOS_INLINE_FUNCTION
BremHR98EOSTerms BremHR98PrecomputeEOS(const MyEOSParams* eos_params,
                                       const bool use_NN_medium_corr)
{
    [[maybe_unused]] constexpr BS_REAL zero = 0;
    constexpr BS_REAL one    = 1;
    constexpr BS_REAL two    = 2;
    constexpr BS_REAL three  = 3;
    constexpr BS_REAL four   = 4;
    constexpr BS_REAL five   = 5;
    constexpr BS_REAL six    = 6;
    constexpr BS_REAL twelve = 12;
    constexpr BS_REAL half   = 0.5;

    constexpr BS_REAL twentyeight_thirds = 28. / 3.;
    constexpr BS_REAL one_twentieth      = 0.05;
    constexpr BS_REAL one_tenth          = 0.1;
    constexpr BS_REAL one_fifth          = 0.2;
    constexpr BS_REAL two_fifth          = 0.4;
    constexpr BS_REAL four_fifth         = 0.8;
    constexpr BS_REAL six_tenth          = 0.6;
    constexpr BS_REAL twelve_tenth       = 1.2;
    constexpr BS_REAL eighteen_tenth     = 1.8;
    constexpr BS_REAL twentythree_tenth  = 2.3;
    constexpr BS_REAL twentyfive         = 25;
    constexpr BS_REAL two_hundreds       = 0.02;
    constexpr BS_REAL four_hundreds      = 0.04;
    constexpr BS_REAL five_hundreds      = 0.05;
    constexpr BS_REAL six_hundreds       = 0.06;
    constexpr BS_REAL fifteen_hundreds   = 0.15;
    constexpr BS_REAL fortyfive_hundreds = 0.45;
    constexpr BS_REAL sixtyseven_hundreds = 0.67;
    constexpr BS_REAL eighteen_hundreds  = 0.18;

    constexpr BS_REAL ten_minus_fourteen = 1.e-14;
    constexpr BS_REAL ten_minus_ten      = 1.e-10;
    constexpr BS_REAL ten_minus_seven    = 1.e-7;

    // gamma and y constants, Eqns. (37) and (38)
    constexpr BS_REAL c1 = 1.63;
    constexpr BS_REAL c2 = 1.94;

    // s-component constants, Eqn. (50)
    constexpr BS_REAL cs1 = 0.0044;
    constexpr BS_REAL cs4 = 2.39;

    // g-component constants, Eqns. (52) and (53)
    constexpr BS_REAL cg1 = 15.6;
    constexpr BS_REAL cg2 = 0.63;
    constexpr BS_REAL cg4 = 0.025;
    constexpr BS_REAL cg5 = 13.75;

    BS_REAL fu_threshold;

    if constexpr (std::is_same_v<BS_REAL, float>)
    {
        fu_threshold = ten_minus_seven;
    }
    else
    {
        fu_threshold = ten_minus_fourteen;
    }

    // EOS parameters
    const BS_REAL nb = eos_params->nb;   // baryon number density [nm^-3]
    const BS_REAL xn = eos_params->yn;   // neutron abundance/mass fraction
    const BS_REAL xp = eos_params->yp;   // proton abundance/mass fraction
    const BS_REAL T  = eos_params->temp; // temperature [MeV]

    // Temperature in units of 10 MeV
    const BS_REAL T_10 = T * one_tenth;

    BremHR98EOSTerms HR98eos_terms;

    HR98eos_terms.temp = T;

    /* y-dependent quantities */

    // Pion mass parameter y, Eqn. (38)
    BS_REAL y = c2 / T_10;

    BS_ASSERT(y >= zero);

    // Prevent singular behavior
    y = (y > kBS_Brem_Ymin) ? y : kBS_Brem_Ymin;

    // Powers of y
    const BS_REAL y2     = POW2(y);
    const BS_REAL log_y  = Kokkos::log(y);
    const BS_REAL y_0_05 = Kokkos::exp(one_twentieth * log_y);
    const BS_REAL y_0_4  = POW4(POW2(y_0_05));
    const BS_REAL y_1_05 = y * y_0_05;

    HR98eos_terms.y           = y;
    HR98eos_terms.s_nd_yexp   = two - SafeExp(-y / twelve);
    HR98eos_terms.s_nd_ypoly  = five * y2 / three + one;
    HR98eos_terms.s_nd_yshift = kBS_Pi2OneEighth + y;
    HR98eos_terms.f_yfac      = one + POW4(y);
    HR98eos_terms.g_yfac      = cs1 * y / (four_fifth + six_hundreds * y_1_05);
    HR98eos_terms.p           = sixtyseven_hundreds + eighteen_hundreds * y_0_4;
    HR98eos_terms.minus_inv_p = -one / HR98eos_terms.p;

    // y-dependent factors of alpha_1, Eqn. (53)
    const BS_REAL alpha_1_denom = twentyfive * y2 + one;
    const BS_REAL alpha_1_a     = one / alpha_1_denom;
    const BS_REAL alpha_1_b     = twentyfive * y2 / alpha_1_denom;
    const BS_REAL y2_over_c5    = y2 / cg5;

    /* Channel-dependent quantities */

    const BS_REAL x_mean =
        Kokkos::sqrt(xn * xp); // geometric mean of nucleon abundances/mass fractions

    // nucleon number densities [nm^-3]: neutron, proton, geometric mean
    const BS_REAL n_nuc[3] = {nb * xn, nb * xp, nb * x_mean};
    // nucleon masses [g]
    const BS_REAL m_nuc[3] = {kBS_MnGrams, kBS_MpGrams, kBS_MAvgGrams};
    // channel multiplicity
    const BS_REAL mult[3] = {one, one, twentyeight_thirds};

    // kernel correction due to medium dependence as in Fischer2016
    const BS_REAL medium_corr =
        use_NN_medium_corr
            ? one / POW6((one + Kokkos::cbrt(nb / kBS_Saturation_n) / three))
            : one;

    for (int ich = 0; ich < 3; ++ich)
    {
        BremHR98Channel* ch = &HR98eos_terms.channel[ich];

        // Nucleon effective degeneracy parameter, Eqn. (36) using baryon number
        // density instead of matter density, with (.)^(2/3) = cbrt(.)^2
        const BS_REAL cbrt_n = Kokkos::cbrt(three * kBS_PiSquared * n_nuc[ich]);
        BS_REAL eta_star =
            POW2(cbrt_n) * kBS_Brem_Aux1 / (two * m_nuc[ich] * T);

        BS_ASSERT(isfinite(eta_star));
        BS_ASSERT(eta_star >= zero);

        // Spin-fluctuation rate gamma, Eqn. (37) (with the unclamped eta_star)
        const BS_REAL gamma = c1 * eta_star * Kokkos::sqrt(eta_star) * T_10;

        // Prevent singular behavior
        eta_star = (eta_star > kBS_Brem_Etamin) ? eta_star : kBS_Brem_Etamin;

        // Powers of eta_star
        const BS_REAL eta2     = POW2(eta_star);
        const BS_REAL sqrt_eta = Kokkos::sqrt(eta_star);
        const BS_REAL eta_1_5  = eta_star * sqrt_eta;
        const BS_REAL eta_2_5  = eta2 * sqrt_eta;
        const BS_REAL eta_0_05 = Kokkos::exp(one_twentieth * Kokkos::log(eta_star));
        const BS_REAL eta_0_4  = POW4(POW2(eta_0_05));
        const BS_REAL eta_1_1  = eta_star * POW2(eta_0_05);
        const BS_REAL eta_1_45 = eta_star * eta_0_4 * eta_0_05;

        /* s-component */

        // f(u) of the degenerate approximation s_d, Eqn. (46)
        const BS_REAL u     = Kokkos::sqrt(y / (two * eta_star)) + ten_minus_ten;
        const BS_REAL u2    = POW2(u);
        const BS_REAL u_arg = u2 / (two * Kokkos::sqrt(two * u2 + four));
        BS_REAL f_u =
            (one - five * u * Kokkos::atan(two / u) / six + u2 / (three * (u2 + four)) +
             Kokkos::atan(one / u_arg) * u_arg / three);

        // @TODO: check this fix! Doing this to prevent s_d from being a large
        // negative number
        f_u = (Kokkos::fabs(f_u) < fu_threshold) ? fu_threshold : f_u;

        ch->s_d_etafac     = f_u / eta_2_5;
        ch->log_s_d_etafac = Kokkos::log(ch->s_d_etafac);

        // eta_star-dependent factors of F, G and h, Eqn. (50)
        ch->f_etafac = one / (one + eta2);
        ch->g_etafac = sqrt_eta / (eta_star + one_fifth);
        ch->h        = one_tenth * eta_star / (cs4 + one_tenth * eta_1_1);

        /* g-component */

        // alpha_1, Eqn. (53)
        const BS_REAL eta_star_inv = one / eta_star;
        const BS_REAL alpha_1 = (half + eta_star_inv) / (one + eta_star_inv) * alpha_1_a +
                                (half + eta_star / cg1) * alpha_1_b;

        // alpha_2, Eqn. (53)
        const BS_REAL alpha_2 =
            (cg2 + four_hundreds * eta_1_45) / (one + two_hundreds * eta_2_5);

        // alpha_3, Eqn. (53)
        const BS_REAL alpha_3 =
            twelve_tenth * SafeExp(six_tenth * eta_star - two_fifth * eta_1_5);

        // p_1, Eqn. (53)
        const BS_REAL p_1 =
            (eighteen_tenth + fortyfive_hundreds * eta_star) / (one + fifteen_hundreds * eta_1_5);

        // p_2, Eqn. (53)
        const BS_REAL p_2 = twentythree_tenth - five_hundreds * eta_star / (one + cg4 * eta_star);

        // g, Eqn. (52), with y^(p_1+2) = y^p_1 * y^2
        const BS_REAL y_p1 = Kokkos::exp(p_1 * log_y);
        const BS_REAL y_p2 = Kokkos::exp(p_2 * log_y);

        const BS_REAL g_brem = (alpha_1 + alpha_2 * y_p1) /
                               (one + alpha_3 * y_p2 + alpha_2 * y_p1 * y2_over_c5);

        BS_ASSERT(g_brem >= zero);

        /* Prefactors of the differential absorption kernel, Eqn. (35) */

        ch->weight = kBS_Brem_Const * mult[ich] * n_nuc[ich] * gamma / T * medium_corr;
        ch->half_gamma_g_sq = POW2(half * gamma * g_brem);
    }

    return HR98eos_terms;
}

/* Compute the angular independent part of the absorption kernel for the
 * Bremsstrahlung reactions, summing the contributions of all NN channels,
 * from the precomputed EOS-dependent quantities.
 *
 * Note: Does not support negative x!
 *
 * Inputs:
 *      HR98eos_terms: EOS-dependent quantities (BremHR98PrecomputeEOS)
 *      x:             rescaled total neutrino energy (w+wp)/T
 *      exp_minus_x:   exp(-x)
 */
KOKKOS_INLINE_FUNCTION
BS_REAL BremAllChannelsAbsKernelFromEOS(const BremHR98EOSTerms* HR98eos_terms, BS_REAL x,
                                        const BS_REAL exp_minus_x)
{
    [[maybe_unused]] constexpr BS_REAL zero = 0;
    constexpr BS_REAL one    = 1;
    constexpr BS_REAL two    = 2;
    constexpr BS_REAL three  = 3;
    constexpr BS_REAL thirty = 30;

    constexpr BS_REAL one_tenth            = 0.1;
    constexpr BS_REAL eleven_tenth         = 1.1;
    constexpr BS_REAL twelve_tenth         = 1.2;
    constexpr BS_REAL twentythree_tenth    = 2.3;
    constexpr BS_REAL five_thousands       = 0.005;
    constexpr BS_REAL ninetythree_hundreds = 0.93;

    // s-component constant, Eqn. (50)
    constexpr BS_REAL cs3 = 0.0001;

    BS_ASSERT(x >= zero);

    const BS_REAL x2_noclamp = POW2(x);

    // Prevent singular behavior (1 - exp(-x) ~ x for x -> 0)
    const BS_REAL one_minus_exp_minus_x =
        (x > kBS_Brem_Xmin) ? one - exp_minus_x : kBS_Brem_Xmin;
    x = (x > kBS_Brem_Xmin) ? x : kBS_Brem_Xmin;

    // Powers of x
    const BS_REAL x2       = POW2(x);
    const BS_REAL log_x    = Kokkos::log(x);
    const BS_REAL x_0_1    = Kokkos::exp(one_tenth * log_x);
    const BS_REAL x_0_2    = POW2(x_0_1);
    const BS_REAL x_0_8    = POW4(x_0_2);
    const BS_REAL x_1_1    = x * x_0_1;
    const BS_REAL x_1_2    = x * x_0_2;
    const BS_REAL x_2_8    = x2 * x_0_8;
    const BS_REAL x_0_93   = Kokkos::exp(ninetythree_hundreds * log_x);
    const BS_REAL x_minus4 = one / POW2(x2);

    // Non-degenerate approximation, s_nd in Eqn. (45)
    const BS_REAL s_nd_base = x + HR98eos_terms->s_nd_yexp;
    const BS_REAL s_nd_numerator =
        two * kBS_SqrtPi * s_nd_base * Kokkos::sqrt(s_nd_base) *
        (x * (x + two * HR98eos_terms->y) + HR98eos_terms->s_nd_ypoly);
    const BS_REAL s_nd_denominator = kBS_SqrtPi + POW4(HR98eos_terms->s_nd_yshift + x);
    const BS_REAL s_nd = s_nd_numerator / s_nd_denominator;

    // Degenerate approximation s_d in Eqn. (46) is s_d_xfac * ch.s_d_etafac
    const BS_REAL s_d_xfac = three * kBS_PiHalfToFiveHalves * (x2 + kBS_FourPiSquared) * x /
                             (kBS_FourPiSquared * one_minus_exp_minus_x);

    // log(s_nd / s_d) = log_s_nd_over_xfac - ch.log_s_d_etafac
    const BS_REAL log_s_nd_over_xfac = Kokkos::log(s_nd / s_d_xfac);

    // x-dependent factors of F, G and C, Eqn. (50)
    const BS_REAL f_xy_denom = (three + POW2(x - twelve_tenth) + x_minus4) * HR98eos_terms->f_yfac;
    const BS_REAL g_xy_fac   = x_1_1 * HR98eos_terms->g_yfac;
    const BS_REAL c_xfac     = eleven_tenth * x_1_1 * thirty / (thirty + five_thousands * x_2_8);
    const BS_REAL c_xdenom   = twentythree_tenth + cs3 * x_1_2;

    BS_REAL s_abs_tot = zero;

    for (int ich = 0; ich < 3; ++ich)
    {
        const BremHR98Channel* ch = &HR98eos_terms->channel[ich];

        const BS_REAL s_d = s_d_xfac * ch->s_d_etafac;

        // F, G, h and C, Eqn. (50)
        const BS_REAL f_brem = one + ch->f_etafac / f_xy_denom;
        const BS_REAL g_brem = one - g_xy_fac * ch->g_etafac;
        const BS_REAL c_brem = c_xfac * ch->h / (c_xdenom + ch->h * x_0_93);

        // interpolated formula for s in Eqn. (49), written as
        // s_min * (1 + (s_min/s_max)^p)^(-1/p), with
        // (s_min/s_max)^p = exp(-p |log(s_nd/s_d)|)
        const BS_REAL log_s_nd_over_s_d = log_s_nd_over_xfac - ch->log_s_d_etafac;
        const BS_REAL s_min   = (log_s_nd_over_s_d < zero) ? s_nd : s_d;
        const BS_REAL ratio_p = Kokkos::exp(-HR98eos_terms->p * Kokkos::fabs(log_s_nd_over_s_d));
        // (log1p is not needed: for ratio_p << 1 the factor is ~ 1 - ratio_p / p and the
        // absolute rounding error of log(1 + ratio_p) is harmless)
        const BS_REAL s_brem = s_min * Kokkos::exp(HR98eos_terms->minus_inv_p * Kokkos::log(one + ratio_p)) *
                               f_brem * (one + c_brem * g_brem);

        BS_ASSERT(s_brem >= zero);

        // Differential absorption kernel, Eqn. (35)
        s_abs_tot += ch->weight * s_brem / (x2_noclamp + ch->half_gamma_g_sq);
    }

    return s_abs_tot;
}

/* Compute a specific Legendre coefficient in the expansion of production and
 * absorption kernels for the Bremsstrahlung reactions, from the precomputed
 * EOS-dependent quantities.
 *
 * Inputs:
 *      kernel_params: omega, omega_prime [MeV] and l
 *      HR98eos_terms: EOS-dependent quantities (BremHR98PrecomputeEOS)
 */
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsLegCoeffFromEOS(const BremKernelParams* kernel_params,
                                          const BremHR98EOSTerms* HR98eos_terms)
{
    constexpr BS_REAL one   = 1;
    constexpr BS_REAL three = 3;

    // kernel parameters
    const int l         = kernel_params->l;     // order of Legendre coefficient
    const BS_REAL omega = kernel_params->omega; // neutrino energy [MeV]
    const BS_REAL omega_prime =
        kernel_params->omega_prime; // primed neutrino energy [MeV]

    BS_ASSERT(l >= 0 && l <= 1);

    // dimensionless neutrino energy sum
    const BS_REAL x = (omega + omega_prime) / HR98eos_terms->temp;

    // shared by the s-component and the detailed balance
    const BS_REAL exp_minus_x = SafeExp(-x);

    // angular independent part of absorption kernel
    BS_REAL s_abs = BremAllChannelsAbsKernelFromEOS(HR98eos_terms, x, exp_minus_x);

    switch (l)
    {
    case 0:
        s_abs = three * s_abs; // zeroth Legedre coefficient
        break;
    case 1:
        s_abs = -one * s_abs; // first Legedre coefficient
        break;
    default:
        BS_PRINTF(
            "BremKernelsLegCoeff (kernel_brem.c): l = %d must be either 0 "
            "or 1\n",
            l);
    }

    // production kernel from detailed balance
    BS_REAL s_em = s_abs * exp_minus_x;

    MyKernelOutput brem_kernel;

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}

/* Compute a specific Legendre coefficient in the expansion of production and
 * absorption kernels for the Bremsstrahlung reactions.
 *
 * Standalone version: inside loops over neutrino energies call
 * BremHR98PrecomputeEOS once and BremKernelsLegCoeffFromEOS for each energy
 * pair instead. */
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsLegCoeff(BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
    const BremHR98EOSTerms HR98eos_terms =
        BremHR98PrecomputeEOS(eos_params, kernel_params->use_NN_medium_corr);
    return BremKernelsLegCoeffFromEOS(kernel_params, &HR98eos_terms);
}


// ------------------------------------------------------------------------ //
// ----------------------- NON-OPTIMIZED HR98 CODE ------------------------ //
// ------------------------------------------------------------------------ //

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
/*
KOKKOS_INLINE_FUNCTION
BS_REAL BremKernelS(BS_REAL x, BS_REAL y, BS_REAL eta_star)
{
    [[maybe_unused]] constexpr BS_REAL zero   = 0;
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

    constexpr BS_REAL ten_minus_fourteen = 1.e-14;
    constexpr BS_REAL ten_minus_ten      = 1.e-10;
    constexpr BS_REAL ten_minus_seven    = 1.e-7;

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
        two * kBS_SqrtPi * Kokkos::pow(x + two - SafeExp(-y / twelve), three_halves) *
        (POW2(x) + two * x * y + five * POW2(y) / three + one);
    const BS_REAL s_nd_denominator =
        kBS_SqrtPi + POW4(kBS_Pi2OneEighth + x + y);
    const BS_REAL s_nd = s_nd_numerator / s_nd_denominator;

    // Compute degenerate approximation, s_d in Eqn. (46)
    const BS_REAL u     = Kokkos::sqrt(y / (two * eta_star)) + ten_minus_ten;
    const BS_REAL u2    = POW2(u);
    const BS_REAL u_arg = u2 / (two * Kokkos::sqrt(two * u2 + four));
    BS_REAL f_u =
        (one - five * u * Kokkos::atan(two / u) / six + u2 / (three * (u2 + four)) +
         Kokkos::atan(one / u_arg) * u_arg / three);

    BS_REAL fu_threshold;

    if constexpr (std::is_same_v<BS_REAL, float>)
    {
        fu_threshold = ten_minus_seven;
    }
    else
    {
        fu_threshold = ten_minus_fourteen;
    }

    // @TODO: check this fix! Doing this to prevent s_d from being a large
    // negative number
    f_u = (Kokkos::fabs(f_u) < fu_threshold) ? fu_threshold : f_u;

    const BS_REAL s_d = three * kBS_PiHalfToFiveHalves *
                        Kokkos::pow(eta_star, -five_halves) *
                        (POW2(x) + kBS_FourPiSquared) * x * f_u /
                        (kBS_FourPiSquared * (one - SafeExp(-x)));

    const BS_REAL pow_x_1_1 = Kokkos::pow(x, eleven_tenth);

    // F, Eqn. (50)
    const BS_REAL f_denominator =
        (three + POW2(x - twelve_tenth) + Kokkos::pow(x, -four)) *
        (one + POW2(eta_star)) * (one + POW4(y));
    // Eq.(50)
    const BS_REAL f_brem = one + one / f_denominator;

    // G, Eqn. (50)
    const BS_REAL g_brem = one - c1 * pow_x_1_1 * y /
                                     (four_fifth + six_hundreds * Kokkos::pow(y, c2)) *
                                     Kokkos::sqrt(eta_star) / (eta_star + one_fifth);

    // h and C, Eqn. (50)
    const BS_REAL h_brem =
        one_tenth * eta_star / (c4 + one_tenth * Kokkos::pow(eta_star, eleven_tenth));
    const BS_REAL c_brem =
        eleven_tenth * pow_x_1_1 * h_brem /
        (twentythree_tenth + h_brem * Kokkos::pow(x, ninetythree_hundreds) +
         c3 * Kokkos::pow(x, twelve_tenth)) *
        thirty / (thirty + five_thousands * Kokkos::pow(x, fourteen_fifth));

    // p, Eqn. (50)
    const BS_REAL p_brem =
        sixtyseven_hundreds + eighteen_hundreds * Kokkos::pow(y, two_fifth);

    // interpolated formula for s in Eqn. (49)
    const BS_REAL s_brem =
        Kokkos::pow(Kokkos::pow(s_nd, -p_brem) + Kokkos::pow(s_d, -p_brem), -one / p_brem) * f_brem *
        (one + c_brem * g_brem);

    BS_ASSERT(s_brem >= zero);

    return s_brem;
}
*/

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
/*
KOKKOS_INLINE_FUNCTION
BS_REAL BremKernelG(BS_REAL y, BS_REAL eta_star)
{
    [[maybe_unused]] constexpr BS_REAL zero       = 0;
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
    const BS_REAL alpha_2 = (c2 + four_hundreds * Kokkos::pow(eta_star, c3)) /
                            (one + two_hundreds * Kokkos::pow(eta_star, five_halves));

    const BS_REAL pow_eta_star_1_5 = Kokkos::pow(eta_star, three_halves);

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
        (alpha_1 + alpha_2 * Kokkos::pow(y, p_1)) /
        (one + alpha_3 * Kokkos::pow(y, p_2) + alpha_2 * Kokkos::pow(y, p_1 + two) / c5);

    BS_ASSERT(g >= zero);

    return g;
}
*/

/* Compute the absorption kernels for a given NN Bremsstrahlung channel */
/*
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
    const BS_REAL eta_star = Kokkos::pow(three * kBS_PiSquared * n_nuc, two_thirds) *
                             kBS_Brem_Aux1 / (two * m_nuc * T);

    // Spin-fluctuation rate gamma, Eqn. (37)
    const BS_REAL gamma = c1 * Kokkos::pow(eta_star, three_halves) * T_10;

    // Pion mass parameter y, Eqn. (38)
    const BS_REAL y = c2 / T_10;

    // Dimensionless fitting parameter s
    const BS_REAL sb = BremKernelS(x, y, eta_star);

    // Dimensionless fitting parameter g
    const BS_REAL gb = BremKernelG(y, eta_star);

    // Differential absorption kernel, Eqn. (35)
    return gamma / (POW2(x) + POW2(half * gamma * gb)) * sb / T;
}
*/

/* Compute the angular independent part of the absorption kernels for the
 Bremsstrahlung reactions by summing the contributions of all NN channels */
/*
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
        Kokkos::sqrt(xn * xp); // geometric mean of nucleon abundances/mass fractions

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
            s_abs_tot / POW6((one + Kokkos::cbrt(nb / kBS_Saturation_n) / three));
    }

    return s_abs_tot;
}
*/


#endif // BNS_NURATES_INCLUDE_KERNEL_BREM_HR98_HPP_
