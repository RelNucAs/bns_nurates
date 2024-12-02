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
    BS_ASSERT(x >= 0.);
    BS_ASSERT(y >= 0.);
    BS_ASSERT(isfinite(eta_star));
    BS_ASSERT(eta_star >= 0.);

    // Prevent singular behavior
    x        = (x > kBS_Brem_Xmin) ? x : kBS_Brem_Xmin;
    y        = (y > kBS_Brem_Ymin) ? y : kBS_Brem_Ymin;
    eta_star = (eta_star > kBS_Brem_Etamin) ? eta_star : kBS_Brem_Etamin;

    // Compute non-degenerate approximation, s_nd in Eqn. (45)
    const BS_REAL s_nd_numerator =
        2. * kBS_SqrtPi * pow(x + 2. - SafeExp(-y / 12.), 1.5) *
        (POW2(x) + 2. * x * y + 5. * POW2(y) / 3. + 1.);
    const BS_REAL s_nd_denominator =
        kBS_SqrtPi + POW4(kBS_Pi2OneEighth + x + y);
    const BS_REAL s_nd = s_nd_numerator / s_nd_denominator;

    // Compute degenerate approximation, s_d in Eqn. (46)
    const BS_REAL u     = sqrt(y / (2. * eta_star)) + 1.e-10;
    const BS_REAL u2    = POW2(u);
    const BS_REAL u_arg = u2 / (2. * sqrt(2. * u2 + 4.));
    BS_REAL f_u = (1. - 5. * u * atan(2. / u) / 6. + u2 / (3. * (u2 + 4.)) +
                   atan(1. / u_arg) * u_arg / 3.);

    // @TODO: Leonardo check this! Doing this to prevent s_d from being a large
    // negative number
    f_u = (fabs(f_u) < 1e-14) ? 1e-14 : f_u;

    const BS_REAL s_d = 3. * kBS_PiHalfToFiveHalves * pow(eta_star, -2.5) *
                        (POW2(x) + kBS_FourPiSquared) * x * f_u /
                        (kBS_FourPiSquared * (1. - SafeExp(-x)));

    const BS_REAL pow_x_1_1 = pow(x, 1.1);

    // F, Eqn. (50)
    const BS_REAL f_denominator = (3. + POW2(x - 1.2) + pow(x, -4.)) *
                                  (1. + POW2(eta_star)) * (1. + POW4(y));
    // Eq.(50)
    const BS_REAL f_brem = 1. + 1. / f_denominator;

    // G, Eqn. (50)
    const BS_REAL g_brem = 1. - 0.0044 * pow_x_1_1 * y /
                                    (0.8 + 0.06 * pow(y, 1.05)) *
                                    sqrt(eta_star) / (eta_star + 0.2);

    // h and C, Eqn. (50)
    const BS_REAL h_brem = 0.1 * eta_star / (2.39 + 0.1 * pow(eta_star, 1.1));
    const BS_REAL c_brem =
        1.1 * pow_x_1_1 * h_brem /
        (2.3 + h_brem * pow(x, 0.93) + 0.0001 * pow(x, 1.2)) * 30. /
        (30. + 0.005 * pow(x, 2.8));

    // p, Eqn. (50)
    const BS_REAL p_brem = 0.67 + 0.18 * pow(y, 0.4);

    // interpolated formula for s in Eqn. (49)
    const BS_REAL s_brem =
        pow(pow(s_nd, -p_brem) + pow(s_d, -p_brem), -1. / p_brem) * f_brem *
        (1. + c_brem * g_brem);

    BS_ASSERT(s_brem >= 0.);

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
    BS_ASSERT(y >= 0.);
    BS_ASSERT(eta_star >= 0.);

    // prevent singular behavior
    y        = (y > kBS_Brem_Ymin) ? y : kBS_Brem_Ymin;
    eta_star = (eta_star > kBS_Brem_Etamin) ? eta_star : kBS_Brem_Etamin;

    // alpha_1, Eqn. (53)
    const BS_REAL y2            = POW2(y);
    const BS_REAL eta_star_inv  = 1. / eta_star;
    const BS_REAL alpha_1_denom = 25. * y2 + 1.;

    const BS_REAL alpha_1 =
        (0.5 + eta_star_inv) / (1. + eta_star_inv) * (1. / alpha_1_denom) +
        (0.5 + eta_star / 15.6) * 25. * y2 / alpha_1_denom;

    // alpha_2, Eqn. (53)
    const BS_REAL alpha_2 =
        (0.63 + 0.04 * pow(eta_star, 1.45)) / (1. + 0.02 * pow(eta_star, 2.5));

    const BS_REAL pow_eta_star_1_5 = pow(eta_star, 1.5);

    // alpha_3, Eqn. (53)
    const BS_REAL alpha_3 =
        1.2 * SafeExp(0.6 * eta_star - 0.4 * pow_eta_star_1_5);

    // p_1, Eqn. (53)
    const BS_REAL p_1 =
        (1.8 + 0.45 * eta_star) / (1. + 0.15 * pow_eta_star_1_5);

    // p_2, Eqn. (53)
    const BS_REAL p_2 = 2.3 - 0.05 * eta_star / (1. + 0.025 * eta_star);

    // g, Eqn. (52)
    const BS_REAL g =
        (alpha_1 + alpha_2 * pow(y, p_1)) /
        (1. + alpha_3 * pow(y, p_2) + alpha_2 * pow(y, p_1 + 2.) / 13.75);

    BS_ASSERT(g >= 0.);

    return g;
}

/* Compute the absorption kernels for a given NN Bremsstrahlung channel */
KOKKOS_INLINE_FUNCTION
BS_REAL BremSingleChannelAbsKernel(const BS_REAL n_nuc, const BS_REAL m_nuc,
                                   BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
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
    const BS_REAL T_10 = T * 0.1;

    // Nucleon effective degeneracy parameter, Eqn. (36) using baryon number
    // density instead of matter density
    const BS_REAL eta_star = pow(3. * kBS_PiSquared * n_nuc, 2. / 3.) *
                             kBS_Brem_Aux1 / (2. * m_nuc * T);

    // Spin-fluctuation rate gamma, Eqn. (37)
    const BS_REAL gamma = 1.63 * pow(eta_star, 3. / 2.) * T_10;

    // Pion mass parameter y, Eqn. (38)
    const BS_REAL y = 1.94 / T_10;

    // Dimensionless fitting parameter s
    const BS_REAL sb = BremKernelS(x, y, eta_star);

    // Dimensionless fitting parameter g
    const BS_REAL gb = BremKernelG(y, eta_star);

    // Differential absorption kernel, Eqn. (35)
    return gamma / (POW2(x) + POW2(0.5 * gamma * gb)) * sb / T;
}

/* Compute the angular independent part of the absorption kernels for the
 Bremsstrahlung reactions by summing the contributions of all NN channels */
KOKKOS_INLINE_FUNCTION
BS_REAL BremAllChannelsAbsKernel(BremKernelParams* kernel_params,
                                 MyEOSParams* eos_params)
{
    static const BS_REAL kTwentyeightThirds = 28. / 3.;

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
    BS_REAL s_abs_tot = kBS_Brem_Const * (nn * s_abs_nn + np * s_abs_pp +
                                          28. * n_mean * s_abs_np / 3.);

    // kernel correction due to medium dependence as in Fischer2016
    if (kernel_params->use_NN_medium_corr == true)
    {
        s_abs_tot = s_abs_tot / POW6((1. + cbrt(nb / kBS_Saturation_n) / 3.));
    }

    return s_abs_tot;
}

/* Compute a specific Legendre coefficient in the expansion of production and
 * absorption kernels for the Bremsstrahlung reactions */
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsLegCoeff(BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
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
        s_abs = 3. * s_abs; // zeroth Legedre coefficient
        break;
    case 1:
        s_abs = -1. * s_abs; // first Legedre coefficient
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
    const BS_REAL rho = nb * kBS_Mb; // mass density [g nm-3]
    return 2.0778e+02 * 0.5 * kBS_MeV *
           (POW2(xn) + POW2(xp) + 28. * xn * xp / 3.) * POW2(rho) *
           pow(T, 5.5); // [MeV nm-3 s-1]
}

// Bremsstrahlung kernel from BRT06 Eq.(143) rewritten consistently
// to fit within the framework of the present library
KOKKOS_INLINE_FUNCTION
MyKernelOutput BremKernelsBRT06(BremKernelParams* kernel_params,
                                MyEOSParams* eos_pars)
{
    const BS_REAL omega       = kernel_params->omega;
    const BS_REAL omega_prime = kernel_params->omega_prime;
    const BS_REAL temp        = eos_pars->temp;

    const BS_REAL x = 0.5 * (omega + omega_prime) / temp;
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
