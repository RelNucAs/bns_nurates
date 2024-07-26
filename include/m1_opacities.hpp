#ifndef BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  m1_opacities.hpp
//  \brief header file for all integration routines

#include <Kokkos_Core.hpp>

#include "bns_nurates.hpp"
#include "opacities.hpp"
#include "kernels.hpp"
#include "integration.hpp"
#include "distribution.hpp"

/* Compute the 2d integrands for all reactions from Leonardo's notes [Eqns. (51)
 * & (52)] There are a total of two expressions for 'e' and 'x' neutrinos, so 4
 * integrands in total
 *
 * 1. Contribution to emissivity: (4 pi)^2/(hc)^6 nu^3 nubar^2 [R^prod(pair) +
 * R^prod(brem)][1 - g_nubar]
 * 2. Contribution to absorption coefficient: (1/(c J)) *(4 pi)^2/(hc)^6 * (nu^3
 * nubar^2 [R_pro(Pair) + R_pro(Brem)][1 - g_nubar] g_nu
 *                                                                        + nu^3
 * nubar^2 [R_abs(Pair) + R_abs(Brem)]g_nubar g_nu)
 *
 * Note that there are no double integrals for the computation of the scattering
 * coefficient.
 */
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand M1CoeffsDoubleIntegrand(double* var, void* p)
{

    // energies and parameters
    double nu     = var[0]; // [MeV]
    double nu_bar = var[1]; // [MeV]

    GreyOpacityParams* my_grey_opacity_params = (GreyOpacityParams*)p;
    OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species];
    double g_nu_bar[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
        g_nu_bar[idx] =
            TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
    }

    // compute the pair kernels
    MyKernelOutput pair_kernels_m1 = {
        0}; //{.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
    if (opacity_flags.use_pair)
    {
        my_grey_opacity_params->kernel_pars.pair_kernel_params.omega = nu;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.omega_prime =
            nu_bar;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = 1.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.filter    = 0.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax      = 0;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu        = 1.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime  = 1.;
        // pair_kernels_m1 = PairKernelsM1(&my_grey_opacity_params->eos_pars,
        // &my_grey_opacity_params->kernel_pars.pair_kernel_params);
        pair_kernels_m1 = PairKernelsOptimized(
            &my_grey_opacity_params->eos_pars,
            &my_grey_opacity_params->kernel_pars.pair_kernel_params);
    }

    // compute the bremsstrahlung kernels
    MyKernelOutput brem_kernels_m1 = {0};
    if (opacity_flags.use_brem)
    {
        my_grey_opacity_params->kernel_pars.brem_kernel_params.omega = nu;
        my_grey_opacity_params->kernel_pars.brem_kernel_params.omega_prime =
            nu_bar;
        my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
        brem_kernels_m1 = BremKernelsLegCoeff(
            &my_grey_opacity_params->kernel_pars.brem_kernel_params,
            &my_grey_opacity_params->eos_pars);
    }

    // compute the inelastic NES/NPS kernels
    MyKernelOutput inelastic_kernels_m1 = {0};
    if (opacity_flags.use_inelastic_scatt)
    {
        my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;
        my_grey_opacity_params->kernel_pars.inelastic_kernel_params
            .omega_prime     = nu_bar;
        inelastic_kernels_m1 = InelasticScattKernels(
            &my_grey_opacity_params->kernel_pars.inelastic_kernel_params,
            &my_grey_opacity_params->eos_pars);
    }

    double pro_term[total_num_species];

    pro_term[id_nue] =
        (pair_kernels_m1.em[id_nue] + brem_kernels_m1.em[id_nue]) *
            (1. - g_nu_bar[id_anue]) +
        inelastic_kernels_m1.em[id_nue] * g_nu_bar[id_nue];
    pro_term[id_anue] =
        (pair_kernels_m1.em[id_anue] + brem_kernels_m1.em[id_anue]) *
            (1. - g_nu_bar[id_nue]) +
        inelastic_kernels_m1.em[id_anue] * g_nu_bar[id_anue];
    pro_term[id_nux] =
        (pair_kernels_m1.em[id_nux] + brem_kernels_m1.em[id_nux]) *
            (1. - g_nu_bar[id_anux]) +
        inelastic_kernels_m1.em[id_nux] * g_nu_bar[id_nux];
    pro_term[id_anux] =
        (pair_kernels_m1.em[id_anux] + brem_kernels_m1.em[id_anux]) *
            (1. - g_nu_bar[id_nux]) +
        inelastic_kernels_m1.em[id_anux] * g_nu_bar[id_anux];

    // pro_term[id_nue] = (pair_kernels_m1.em_e + brem_kernels_m1.em[id_nue]) +
    // inelastic_kernels_m1.em[id_nue] * g_nu_bar[id_nue]; pro_term[id_anue] =
    // (pair_kernels_m1.em_e + brem_kernels_m1.em[id_anue]) +
    // inelastic_kernels_m1.em[id_anue] * g_nu_bar[id_anue]; pro_term[id_nux] =
    // (pair_kernels_m1.em_x + brem_kernels_m1.em[id_nux]) +
    // inelastic_kernels_m1.em[id_nux] * g_nu_bar[id_nux]; pro_term[id_anux] =
    // (pair_kernels_m1.em_x + brem_kernels_m1.em[id_anux]) +
    // inelastic_kernels_m1.em[id_anux] * g_nu_bar[id_anux];

    double ann_term[total_num_species];

    ann_term[id_nue] =
        (pair_kernels_m1.abs[id_nue] + brem_kernels_m1.abs[id_nue]) *
            g_nu_bar[id_anue] +
        inelastic_kernels_m1.abs[id_nue] * (1. - g_nu_bar[id_nue]);
    ann_term[id_anue] =
        (pair_kernels_m1.abs[id_anue] + brem_kernels_m1.abs[id_anue]) *
            g_nu_bar[id_nue] +
        inelastic_kernels_m1.abs[id_anue] * (1. - g_nu_bar[id_anue]);
    ann_term[id_nux] =
        (pair_kernels_m1.abs[id_nux] + brem_kernels_m1.abs[id_nux]) *
            g_nu_bar[id_anux] +
        inelastic_kernels_m1.abs[id_nux] * (1. - g_nu_bar[id_nux]);
    ann_term[id_anux] =
        (pair_kernels_m1.abs[id_anux] + brem_kernels_m1.abs[id_anux]) *
            g_nu_bar[id_nux] +
        inelastic_kernels_m1.abs[id_anux] * (1. - g_nu_bar[id_anux]);

    double n_integrand_1[total_num_species], e_integrand_1[total_num_species];
    double n_integrand_2[total_num_species], e_integrand_2[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        n_integrand_1[idx] = nu * nu * nu_bar * nu_bar * pro_term[idx];
        e_integrand_1[idx] = nu * n_integrand_1[idx];

        n_integrand_2[idx] = g_nu[idx] * nu * nu * nu_bar * nu_bar *
                             (pro_term[idx] + ann_term[idx]);
        // n_integrand_2[idx] = g_nu[idx] * nu * nu * nu_bar * nu_bar *
        // ann_term[idx];
        e_integrand_2[idx] = nu * n_integrand_2[idx];
    }

    MyQuadratureIntegrand result = {.n = 16};

    result.integrand[0]  = n_integrand_1[id_nue];
    result.integrand[1]  = n_integrand_1[id_anue];
    result.integrand[2]  = n_integrand_1[id_nux];
    result.integrand[3]  = n_integrand_1[id_anux];
    result.integrand[4]  = n_integrand_2[id_nue];
    result.integrand[5]  = n_integrand_2[id_anue];
    result.integrand[6]  = n_integrand_2[id_nux];
    result.integrand[7]  = n_integrand_2[id_anux];
    result.integrand[8]  = e_integrand_1[id_nue];
    result.integrand[9]  = e_integrand_1[id_anue];
    result.integrand[10] = e_integrand_1[id_nux];
    result.integrand[11] = e_integrand_1[id_anux];
    result.integrand[12] = e_integrand_2[id_nue];
    result.integrand[13] = e_integrand_2[id_anue];
    result.integrand[14] = e_integrand_2[id_nux];
    result.integrand[15] = e_integrand_2[id_anux];

    return result;
}

/* Computes the integrand for all single integrals from Leonardo's notes
 *
 * There are a total of 3 expressions for electron-type neutrinos, electron-type
 * antineutrinos and 'x' neutrinos, so a total of 9 integrands should be
 * computed
 *
 * However, two of them (those in Eq.(51) and Eq.(52) for 'x' neutrinos) are
 * trivially equal to zero)
 *
 * 1. Contribution to emissivity: (4 pi /(h c)^3) nu^3 j_x
 * 2. Contribution to absorption coefficient: (1/(c J)) (4 pi /(h c)^3) nu^3
 * g_nu (j_x + 1/lambda_x)
 * 3. Contribution to scattering coefficient: (1/(c J)) (4 pi)^2 nu^5 g_nu
 * (R_iso(1)/3 - R_iso(0))
 */

/* Computes the integrand for all single integrals from Leonardo's notes
 *
 * There are a total of 3 expressions for electron-type neutrinos, electron-type
 * antineutrinos and 'x' neutrinos, so a total of 9 integrands should be
 * computed
 *
 * However, two of them (those in Eq.(51) and Eq.(52) for 'x' neutrinos) are
 * trivially equal to zero)
 *
 * 1. Contribution to emissivity: (4 pi /(h c)^3) nu^3 j_x
 * 2. Contribution to absorption coefficient: (1/(c J)) (4 pi /(h c)^3) nu^3
 * g_nu (j_x + 1/lambda_x)
 * 3. Contribution to scattering coefficient: (1/(c J)) (4 pi)^2 nu^5 g_nu
 * (R_iso(1)/3 - R_iso(0))
 */
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand M1CoeffsSingleIntegrandNotStimulated(double* var, void* p)
{

    // energy and parameters
    double nu = var[0]; // [MeV]

    GreyOpacityParams* my_grey_opacity_params = (GreyOpacityParams*)p;
    OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    }

    double n_integrand_1[total_num_species] = {0.},
           e_integrand_1[total_num_species] = {0.};
    double n_integrand_2[total_num_species] = {0.},
           e_integrand_2[total_num_species] = {0.};
    double e_integrand_3[total_num_species] = {0.};

    if (opacity_flags.use_abs_em == 1)
    {
        MyOpacity abs_em_beta =
            AbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                       &my_grey_opacity_params->eos_pars); // [s^-1]

        for (int idx = 0; idx < total_num_species; idx++)
        {
            n_integrand_1[idx] =
                nu * nu * abs_em_beta.em[idx];            // [MeV^-1 cm^-3 s^-1]
            e_integrand_1[idx] = nu * n_integrand_1[idx]; // [cm^-3 s^-1]

            n_integrand_2[idx] =
                nu * nu * g_nu[idx] *
                abs_em_beta.abs[idx]; // ab = ab (NOT stimulated absorption)
            if (idx == 0)
            {
                double tmp = my_grey_opacity_params->eos_pars.mu_n -
                             my_grey_opacity_params->eos_pars.mu_p -
                             my_grey_opacity_params->eos_pars.mu_e;
            }
            e_integrand_2[idx] = nu * n_integrand_2[idx];
        }
    }

    if (opacity_flags.use_iso == 1)
    {
        const double iso_scatt =
            IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars,
                          &my_grey_opacity_params->eos_pars);

        for (int idx = 0; idx < total_num_species; idx++)
        {
            e_integrand_3[idx] =
                4. * kPi * nu * nu * nu * nu * nu * g_nu[idx] * iso_scatt;
        }
    }

    MyQuadratureIntegrand result = {.n = 12};

    result.integrand[0]  = n_integrand_1[id_nue];
    result.integrand[1]  = n_integrand_1[id_anue];
    result.integrand[2]  = n_integrand_2[id_nue];
    result.integrand[3]  = n_integrand_2[id_anue];
    result.integrand[4]  = e_integrand_1[id_nue];
    result.integrand[5]  = e_integrand_1[id_anue];
    result.integrand[6]  = e_integrand_2[id_nue];
    result.integrand[7]  = e_integrand_2[id_anue];
    result.integrand[8]  = e_integrand_3[id_nue];
    result.integrand[9]  = e_integrand_3[id_anue];
    result.integrand[10] = e_integrand_3[id_nux];
    result.integrand[11] = e_integrand_3[id_anux];

    return result;
}

KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand M1CoeffsSingleIntegrand(double* var, void* p)
{

    // energy and parameters
    double nu = var[0]; // [MeV]

    GreyOpacityParams* my_grey_opacity_params = (GreyOpacityParams*)p;
    OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    }

    double n_integrand_1[total_num_species] = {0.},
           e_integrand_1[total_num_species] = {0.};
    double n_integrand_2[total_num_species] = {0.},
           e_integrand_2[total_num_species] = {0.};
    double e_integrand_3[total_num_species] = {0.};

    if (opacity_flags.use_abs_em == 1)
    {
        MyOpacity abs_em_beta =
            StimAbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                           &my_grey_opacity_params->eos_pars); // [s^-1]

        for (int idx = 0; idx < total_num_species; idx++)
        {
            n_integrand_1[idx] =
                nu * nu * abs_em_beta.em[idx];            // [MeV^-1 cm^-3 s^-1]
            e_integrand_1[idx] = nu * n_integrand_1[idx]; // [cm^-3 s^-1]

            n_integrand_2[idx] =
                nu * nu * g_nu[idx] *
                abs_em_beta.abs[idx]; // ab = em + ab (stimulated absorption)
            e_integrand_2[idx] = nu * n_integrand_2[idx];
        }
    }

    if (opacity_flags.use_iso == 1)
    {
        const double iso_scatt =
            IsoScattTotal(nu, &my_grey_opacity_params->opacity_pars,
                          &my_grey_opacity_params->eos_pars);

        for (int idx = 0; idx < total_num_species; idx++)
        {
            e_integrand_3[idx] =
                4. * kPi * nu * nu * nu * nu * nu * g_nu[idx] * iso_scatt;
        }
    }

    MyQuadratureIntegrand result = {.n = 12};

    result.integrand[0]  = n_integrand_1[id_nue];
    result.integrand[1]  = n_integrand_1[id_anue];
    result.integrand[2]  = n_integrand_2[id_nue];
    result.integrand[3]  = n_integrand_2[id_anue];
    result.integrand[4]  = e_integrand_1[id_nue];
    result.integrand[5]  = e_integrand_1[id_anue];
    result.integrand[6]  = e_integrand_2[id_nue];
    result.integrand[7]  = e_integrand_2[id_anue];
    result.integrand[8]  = e_integrand_3[id_nue];
    result.integrand[9]  = e_integrand_3[id_anue];
    result.integrand[10] = e_integrand_3[id_nux];
    result.integrand[11] = e_integrand_3[id_anux];

    return result;
}

KOKKOS_INLINE_FUNCTION
void ComputeM1DoubleIntegrandNotStimulated(MyQuadrature* quad_1d,
                                           GreyOpacityParams* grey_pars,
                                           double t, M1MatrixKokkos* out_n,
                                           M1MatrixKokkos* out_j)
{
    const int n = quad_1d->nx;

    double nu, nu_bar;

    double* nu_array;

    for (int i = 0; i < n; i++)
    {
        nu_array[i]     = t * quad_1d->points[i];
        nu_array[n + i] = t / quad_1d->points[i];
    }

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species], g_nu_bar[total_num_species];
    double block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    double pro_term_ij[total_num_species] = {0},
           pro_term_ji[total_num_species] = {0};
    double ann_term_ij[total_num_species] = {0},
           ann_term_ji[total_num_species] = {0};

    OpacityParams opacity_pars = grey_pars->opacity_pars;

    // M1Matrix beta;
    M1MatrixKokkos pair, brem, inel;

    // BetaOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars,
    // &grey_pars->opacity_pars, t, &beta);

    if (grey_pars->opacity_flags.use_pair == 1)
    {
        PairKernelsTable(2 * n, nu_array, grey_pars, &pair);
    }

    if (grey_pars->opacity_flags.use_brem == 1)
    {
        if (grey_pars->opacity_pars.use_BRT_brem == true)
        {
            BremKernelsTableBRT06(2 * n, nu_array, grey_pars, &brem);
        }
        else
        {
            BremKernelsTable(2 * n, nu_array, grey_pars, &brem);
        }
    }

    if (grey_pars->opacity_flags.use_inelastic_scatt == 1)
    {
        InelasticKernelsTable(2 * n, nu_array, grey_pars, &inel);
    }

    for (int i = 0; i < 2 * n; i++)
    {
        nu = nu_array[i];

        for (int j = i; j < 2 * n; j++)
        {
            nu_bar = nu_array[j];

            for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx]     = 1. - g_nu[idx];
                    block_factor_nu_bar[idx] = 1. - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu[idx]     = 1.;
                    block_factor_nu_bar[idx] = 1.;
                }
            }

            pro_term_ij[id_nue] =
                (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) *
                (1. - g_nu[id_nue]) * block_factor_nu_bar[id_anue];
            pro_term_ij[id_anue] =
                (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) *
                (1. - g_nu[id_anue]) * block_factor_nu_bar[id_nue];
            pro_term_ij[id_nux] =
                (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) *
                (1. - g_nu[id_nux]) * block_factor_nu_bar[id_anux];
            pro_term_ij[id_anux] =
                (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) *
                (1. - g_nu[id_anux]) * block_factor_nu_bar[id_nux];

            pro_term_ji[id_nue] =
                (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) *
                (1. - g_nu_bar[id_nue]) * block_factor_nu[id_anue];
            pro_term_ji[id_anue] =
                (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) *
                (1. - g_nu_bar[id_anue]) * block_factor_nu[id_nue];
            pro_term_ji[id_nux] =
                (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) *
                (1. - g_nu_bar[id_nux]) * block_factor_nu[id_anux];
            pro_term_ji[id_anux] =
                (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) *
                (1. - g_nu_bar[id_anux]) * block_factor_nu[id_nux];

            ann_term_ij[id_nue] =
                (pair.m1_mat_ab[id_nue][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_anue];
            ann_term_ij[id_anue] =
                (pair.m1_mat_ab[id_anue][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_nue];
            ann_term_ij[id_nux] =
                (pair.m1_mat_ab[id_nux][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_anux];
            ann_term_ij[id_anux] =
                (pair.m1_mat_ab[id_anux][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_nux];

            ann_term_ji[id_nue] =
                (pair.m1_mat_ab[id_nue][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_anue];
            ann_term_ji[id_anue] =
                (pair.m1_mat_ab[id_anue][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_nue];
            ann_term_ji[id_nux] =
                (pair.m1_mat_ab[id_nux][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_anux];
            ann_term_ji[id_anux] =
                (pair.m1_mat_ab[id_anux][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_nux];

            //////////////////////////////////////////////
            ////// ONLY FOR COMPARISON WITH NULIB ////////
            //////////////////////////////////////////////
            if (kirchoff_flag)
            {
                ann_term_ij[id_nue] +=
                    (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_anue] / g_nu[id_nue];
                ann_term_ij[id_anue] +=
                    (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_nue] / g_nu[id_anue];
                ann_term_ij[id_nux] +=
                    (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_anux] / g_nu[id_nux];
                ann_term_ij[id_anux] +=
                    (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_nux] / g_nu[id_anux];

                ann_term_ji[id_nue] +=
                    (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_anue] / g_nu_bar[id_nue];
                ann_term_ji[id_anue] +=
                    (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_nue] / g_nu_bar[id_anue];
                ann_term_ji[id_nux] +=
                    (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_anux] / g_nu_bar[id_nux];
                ann_term_ji[id_anux] +=
                    (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_nux] / g_nu_bar[id_anux];
            }
            //////////////////////////////////////////////

            for (int idx = 0; idx < total_num_species; idx++)
            {
                pro_term_ij[idx] += inel.m1_mat_em[idx][i][j] * g_nu_bar[idx] *
                                    (1. - g_nu[idx]); // @TODO: check this!
                pro_term_ji[idx] += inel.m1_mat_em[idx][j][i] * g_nu[idx] *
                                    (1. - g_nu_bar[idx]);

                ann_term_ij[idx] +=
                    inel.m1_mat_ab[idx][i][j] * block_factor_nu_bar[idx];
                ann_term_ji[idx] +=
                    inel.m1_mat_ab[idx][j][i] * block_factor_nu[idx];

                out_n->m1_mat_em[idx][i][j] =
                    nu * nu * nu_bar * nu_bar * pro_term_ij[idx];
                out_n->m1_mat_em[idx][j][i] =
                    nu_bar * nu_bar * nu * nu * pro_term_ji[idx];

                out_n->m1_mat_ab[idx][i][j] =
                    nu * nu * g_nu[idx] * nu_bar * nu_bar * ann_term_ij[idx];
                out_n->m1_mat_ab[idx][j][i] = nu_bar * nu_bar * g_nu_bar[idx] *
                                              nu * nu * ann_term_ji[idx];

                out_j->m1_mat_em[idx][i][j] = nu * out_n->m1_mat_em[idx][i][j];
                out_j->m1_mat_em[idx][j][i] =
                    nu_bar * out_n->m1_mat_em[idx][j][i];

                out_j->m1_mat_ab[idx][i][j] = nu * out_n->m1_mat_ab[idx][i][j];
                out_j->m1_mat_ab[idx][j][i] =
                    nu_bar * out_n->m1_mat_ab[idx][j][i];
            }
        }
    }

    //FreeM1Matrix(&pair, n);
    //FreeM1Matrix(&inel, n);

    //FreeM1MatrixSingleFlavor(&brem, n, 0);

    return;
}

KOKKOS_INLINE_FUNCTION
void ComputeM1DoubleIntegrand(MyQuadrature* quad_1d,
                              GreyOpacityParams* grey_pars, double t,
                              M1MatrixKokkos* out_n, M1MatrixKokkos* out_j)
{
    const int n = quad_1d->nx;

    double nu, nu_bar;

    double nu_array[n_max];

    for (int i = 0; i < n; i++)
    {
        nu_array[i]     = t * quad_1d->points[i];
        nu_array[n + i] = t / quad_1d->points[i];
    }


    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species], g_nu_bar[total_num_species];
    double block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    double pro_term_ij[total_num_species] = {0},
           pro_term_ji[total_num_species] = {0};
    double ann_term_ij[total_num_species] = {0},
           ann_term_ji[total_num_species] = {0};

    OpacityParams opacity_pars = grey_pars->opacity_pars;


    // M1Matrix beta;
    M1MatrixKokkos pair, brem, inel;

    // BetaOpacitiesTable(quad_1d, &my_grey_opacity_params->eos_pars,
    // &grey_pars->opacity_pars, t, &beta);

    //InitializeM1Matrix(&pair, n);
    if (grey_pars->opacity_flags.use_pair == 1)
    {
        PairKernelsTable(2 * n, nu_array, grey_pars, &pair);
    }

    //InitializeM1MatrixSingleFlavor(&brem, n, 0);
    if (grey_pars->opacity_flags.use_brem == 1)
    {
        if (grey_pars->opacity_pars.use_BRT_brem == true)
        {
            BremKernelsTableBRT06(2 * n, nu_array, grey_pars, &brem);
        }
        else
        {
            BremKernelsTable(2 * n, nu_array, grey_pars, &brem);
        }
    }

    //InitializeM1Matrix(&inel, n);
    if (grey_pars->opacity_flags.use_inelastic_scatt == 1)
    {
        InelasticKernelsTable(2 * n, nu_array, grey_pars, &inel);
    }

    //InitializeM1Matrix(out_n, n);
    //InitializeM1Matrix(out_j, n);

    for (int i = 0; i < 2 * n; i++)
    {
        nu = nu_array[i];

        for (int j = i; j < 2 * n; j++)
        {
            nu_bar = nu_array[j];

            for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx]     = 1. - g_nu[idx];
                    block_factor_nu_bar[idx] = 1. - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu[idx]     = 1.;
                    block_factor_nu_bar[idx] = 1.;
                }
            }

            pro_term_ij[id_nue] =
                (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) *
                block_factor_nu_bar[id_anue];
            pro_term_ij[id_anue] =
                (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) *
                block_factor_nu_bar[id_nue];
            pro_term_ij[id_nux] =
                (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) *
                block_factor_nu_bar[id_anux];
            pro_term_ij[id_anux] =
                (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) *
                block_factor_nu_bar[id_nux];

            pro_term_ji[id_nue] =
                (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) *
                block_factor_nu[id_anue];
            pro_term_ji[id_anue] =
                (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) *
                block_factor_nu[id_nue];
            pro_term_ji[id_nux] =
                (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) *
                block_factor_nu[id_anux];
            pro_term_ji[id_anux] =
                (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) *
                block_factor_nu[id_nux];

            ann_term_ij[id_nue] =
                (pair.m1_mat_ab[id_nue][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_anue];
            ann_term_ij[id_anue] =
                (pair.m1_mat_ab[id_anue][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_nue];
            ann_term_ij[id_nux] =
                (pair.m1_mat_ab[id_nux][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_anux];
            ann_term_ij[id_anux] =
                (pair.m1_mat_ab[id_anux][i][j] + brem.m1_mat_ab[0][i][j]) *
                g_nu_bar[id_nux];

            ann_term_ji[id_nue] =
                (pair.m1_mat_ab[id_nue][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_anue];
            ann_term_ji[id_anue] =
                (pair.m1_mat_ab[id_anue][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_nue];
            ann_term_ji[id_nux] =
                (pair.m1_mat_ab[id_nux][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_anux];
            ann_term_ji[id_anux] =
                (pair.m1_mat_ab[id_anux][j][i] + brem.m1_mat_ab[0][j][i]) *
                g_nu[id_nux];

            //////////////////////////////////////////////
            ////// ONLY FOR COMPARISON WITH NULIB ////////
            //////////////////////////////////////////////
            if (kirchoff_flag)
            {
                ann_term_ij[id_nue] +=
                    (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_anue] / g_nu[id_nue];
                ann_term_ij[id_anue] +=
                    (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_nue] / g_nu[id_anue];
                ann_term_ij[id_nux] +=
                    (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_anux] / g_nu[id_nux];
                ann_term_ij[id_anux] +=
                    (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) *
                    g_nu_bar[id_nux] / g_nu[id_anux];

                ann_term_ji[id_nue] +=
                    (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_anue] / g_nu_bar[id_nue];
                ann_term_ji[id_anue] +=
                    (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_nue] / g_nu_bar[id_anue];
                ann_term_ji[id_nux] +=
                    (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_anux] / g_nu_bar[id_nux];
                ann_term_ji[id_anux] +=
                    (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) *
                    g_nu[id_nux] / g_nu_bar[id_anux];
            }
            //////////////////////////////////////////////

            for (int idx = 0; idx < total_num_species; idx++)
            {
                pro_term_ij[idx] += inel.m1_mat_em[idx][i][j] * g_nu_bar[idx];
                pro_term_ji[idx] += inel.m1_mat_em[idx][j][i] * g_nu[idx];

                ann_term_ij[idx] +=
                    inel.m1_mat_ab[idx][i][j] * block_factor_nu_bar[idx];
                ann_term_ji[idx] +=
                    inel.m1_mat_ab[idx][j][i] * block_factor_nu[idx];

                out_n->m1_mat_em[idx][i][j] =
                    nu * nu * nu_bar * nu_bar * pro_term_ij[idx];
                out_n->m1_mat_em[idx][j][i] =
                    nu_bar * nu_bar * nu * nu * pro_term_ji[idx];

                out_n->m1_mat_ab[idx][i][j] =
                    nu * nu * g_nu[idx] * nu_bar * nu_bar *
                    (pro_term_ij[idx] + ann_term_ij[idx]);
                out_n->m1_mat_ab[idx][j][i] =
                    nu_bar * nu_bar * g_nu_bar[idx] * nu * nu *
                    (pro_term_ji[idx] + ann_term_ji[idx]);

                out_j->m1_mat_em[idx][i][j] = nu * out_n->m1_mat_em[idx][i][j];
                out_j->m1_mat_em[idx][j][i] =
                    nu_bar * out_n->m1_mat_em[idx][j][i];

                out_j->m1_mat_ab[idx][i][j] = nu * out_n->m1_mat_ab[idx][i][j];
                out_j->m1_mat_ab[idx][j][i] =
                    nu_bar * out_n->m1_mat_ab[idx][j][i];
            }
        }
    }

    //FreeM1Matrix(&pair, n);
    //FreeM1Matrix(&inel, n);

    //FreeM1MatrixSingleFlavor(&brem, n, 0);

    return;
}


/* Computes the opacities for the M1 code
 *
 */
KOKKOS_INLINE_FUNCTION
M1Opacities
ComputeM1OpacitiesNotStimulated(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params)
{

    // compute some constants
    static const double four_pi_hc3 =
        (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
    static const double four_pi_hc3_sqr =
        four_pi_hc3 * four_pi_hc3; // [MeV^-6 cm^-6]

    // set up 1d integration
    MyFunctionMultiD integrand_m1_1d;
    MyQuadratureIntegrand integrand_m1_1d_info = {.n = 12};
    integrand_m1_1d.function = &M1CoeffsSingleIntegrandNotStimulated;
    integrand_m1_1d.dim      = 1;
    integrand_m1_1d.params   = my_grey_opacity_params;
    integrand_m1_1d.my_quadrature_integrand = integrand_m1_1d_info;

    double n[total_num_species];
    double J[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        // m1_pars.n and m1_pars.J are assumed to be parsed in cgs
        n[idx] = my_grey_opacity_params->m1_pars.n[idx];
        J[idx] = my_grey_opacity_params->m1_pars.J[idx];

        J[idx] = J[idx] / kMeV; // erg cm-3 -> MeV cm-3 conversion
    }

    const double temp  = my_grey_opacity_params->eos_pars.temp;
    const double eta_e = my_grey_opacity_params->eos_pars.mu_e / temp;

    const double s     = 1.5 * temp;
    double s_array[12] = {0};

    s_array[0] = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
    s_array[1] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    s_array[2] = s;
    s_array[3] = s;
    s_array[4] = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
    s_array[5] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    s_array[6] = s;
    s_array[7] = s;

    for (int idx = 0; idx < total_num_species; idx++)
    {
        s_array[idx + 8] = J[idx] / n[idx];
    }

    // MyQuadratureIntegrand integrals_1d =
    // GaussLegendreIntegrateFixedSplit1D(quad_1d, &integrand_m1_1d, s);
    MyQuadratureIntegrand integrals_1d =
        GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_array);

    M1MatrixKokkos out_n, out_j;
    ComputeM1DoubleIntegrandNotStimulated(quad_1d, my_grey_opacity_params, s,
                                          &out_n, &out_j);
    MyQuadratureIntegrand n_integrals_2d =
        GaussLegendreIntegrate2DMatrix(quad_1d, &out_n, s);
    MyQuadratureIntegrand e_integrals_2d =
        GaussLegendreIntegrate2DMatrix(quad_1d, &out_j, s);

    M1Opacities m1_opacities;

    m1_opacities.eta_0[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[0] + integrals_1d.integrand[0]);
    m1_opacities.eta_0[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[1] + integrals_1d.integrand[1]);
    m1_opacities.eta_0[id_nux]  = four_pi_hc3_sqr * n_integrals_2d.integrand[2];
    m1_opacities.eta_0[id_anux] = four_pi_hc3_sqr * n_integrals_2d.integrand[3];

    m1_opacities.kappa_0_a[id_nue] =
        four_pi_hc3 / (kClight * n[id_nue]) *
        (four_pi_hc3 * n_integrals_2d.integrand[4] + integrals_1d.integrand[2]);
    m1_opacities.kappa_0_a[id_anue] =
        four_pi_hc3 / (kClight * n[id_anue]) *
        (four_pi_hc3 * n_integrals_2d.integrand[5] + integrals_1d.integrand[3]);
    m1_opacities.kappa_0_a[id_nux] =
        four_pi_hc3_sqr / (kClight * n[id_nux]) * n_integrals_2d.integrand[6];
    m1_opacities.kappa_0_a[id_anux] =
        four_pi_hc3_sqr / (kClight * n[id_anux]) * n_integrals_2d.integrand[7];

    m1_opacities.eta[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[0] + integrals_1d.integrand[4]);
    m1_opacities.eta[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[1] + integrals_1d.integrand[5]);
    m1_opacities.eta[id_nux]  = four_pi_hc3_sqr * e_integrals_2d.integrand[2];
    m1_opacities.eta[id_anux] = four_pi_hc3_sqr * e_integrals_2d.integrand[3];

    m1_opacities.kappa_a[id_nue] =
        four_pi_hc3 / (kClight * J[id_nue]) *
        (four_pi_hc3 * e_integrals_2d.integrand[4] + integrals_1d.integrand[6]);
    m1_opacities.kappa_a[id_anue] =
        four_pi_hc3 / (kClight * J[id_anue]) *
        (four_pi_hc3 * e_integrals_2d.integrand[5] + integrals_1d.integrand[7]);
    m1_opacities.kappa_a[id_nux] =
        four_pi_hc3_sqr / (kClight * J[id_nux]) * e_integrals_2d.integrand[6];
    m1_opacities.kappa_a[id_anux] =
        four_pi_hc3_sqr / (kClight * J[id_anux]) * e_integrals_2d.integrand[7];

    m1_opacities.kappa_s[id_nue] =
        four_pi_hc3 / (kClight * J[id_nue]) * integrals_1d.integrand[8];
    m1_opacities.kappa_s[id_anue] =
        four_pi_hc3 / (kClight * J[id_anue]) * integrals_1d.integrand[9];
    m1_opacities.kappa_s[id_nux] =
        four_pi_hc3 / (kClight * J[id_nux]) * integrals_1d.integrand[10];
    m1_opacities.kappa_s[id_anux] =
        four_pi_hc3 / (kClight * J[id_anux]) * integrals_1d.integrand[11];

    return m1_opacities;
}

KOKKOS_INLINE_FUNCTION
M1Opacities ComputeM1Opacities(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                               GreyOpacityParams* my_grey_opacity_params)
{

    // compute some constants
    static const double four_pi_hc3 =
        (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
    static const double four_pi_hc3_sqr =
        four_pi_hc3 * four_pi_hc3; // [MeV^-6 cm^-6]

    const double temp = my_grey_opacity_params->eos_pars.temp;

    // set up 1d integration
    MyFunctionMultiD integrand_m1_1d;
    MyQuadratureIntegrand integrand_m1_1d_info = {.n = 12};
    integrand_m1_1d.function                   = &M1CoeffsSingleIntegrand;
    integrand_m1_1d.dim                        = 1;
    integrand_m1_1d.params                     = my_grey_opacity_params;
    integrand_m1_1d.my_quadrature_integrand    = integrand_m1_1d_info;

    // set up 2d integration
    MyFunctionMultiD integrand_m1_2d;
    MyQuadratureIntegrand integrand_m1_2d_info = {.n = 16};
    integrand_m1_2d.function                   = &M1CoeffsDoubleIntegrand;
    integrand_m1_2d.dim                        = 2;
    integrand_m1_2d.params                     = my_grey_opacity_params;
    integrand_m1_2d.my_quadrature_integrand    = integrand_m1_2d_info;

    double n[total_num_species];
    double J[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        // m1_pars.n and m1_pars.J are assumed to be parsed in cgs
        n[idx] = my_grey_opacity_params->m1_pars.n[idx];
        J[idx] = my_grey_opacity_params->m1_pars.J[idx];

        J[idx] = J[idx] / kMeV; // erg cm-3 -> MeV cm-3 conversion
    }

    // double s_2d_x[12];
    // double s_2d_y[12];

    // for (int i=0; i<11; i++) {
    // s_2d_x[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
    // s_2d_y[i] = s_2d_x[i];
    //}
    // @TODO: choose this appropriately

    // double s_1d[11];
    // for (int i=0; i<11; i++) {
    // s_1d[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
    //}

    const double eta_e = my_grey_opacity_params->eos_pars.mu_e / temp;

    // const double s = 1.5 * temp;
    const double s = 0.5 * 4.364 * temp;
    // const double s = 0.5 * temp * (FDI_p4(eta_e) / FDI_p3(eta_e) +
    // FDI_p4(-eta_e) / FDI_p3(-eta_e));

    double s_array[12] = {0};
    for (int j = 0; j < 12; j++)
    {
        s_array[j] = 1.5 * temp;
    }

    for (int idx = 0; idx < total_num_species; idx++)
    {
        s_array[idx + 8] = J[idx] / n[idx];
    }

    /*
    s_array[0] = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
    s_array[1] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    s_array[4] = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
    s_array[5] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    */

    /*
    s_array[0] = temp * fmax(eta_e, 4.) - kQ;
    s_array[1] = temp * fmax(-eta_e, 4.) + kQ;
    s_array[4] = temp * fmax(eta_e, 5.) - kQ;
    s_array[5] = temp * fmax(-eta_e, 5.) + kQ;
    */

    // MyQuadratureIntegrand integrals_1d =
    // GaussLegendreIntegrateFixedSplit1D(quad_1d, &integrand_m1_1d, s);
    //MyQuadratureIntegrand integrals_1d =
        GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_array);
    MyQuadratureIntegrand integrals_1d = {0};

    // MyQuadratureIntegrand integrals_2d =
    // GaussLegendreIntegrateFixedSplit2D(quad_2d, &integrand_m1_2d, s);
    M1MatrixKokkos out_n, out_j;
    ComputeM1DoubleIntegrand(quad_1d, my_grey_opacity_params, s, &out_n,
                             &out_j);

    MyQuadratureIntegrand n_integrals_2d =
        GaussLegendreIntegrate2DMatrix(quad_1d, &out_n, s);
    MyQuadratureIntegrand e_integrals_2d =
        GaussLegendreIntegrate2DMatrix(quad_1d, &out_j, s);

    M1Opacities m1_opacities;

    m1_opacities.eta_0[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[0] + integrals_1d.integrand[0]);
    m1_opacities.eta_0[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[1] + integrals_1d.integrand[1]);
    m1_opacities.eta_0[id_nux]  = four_pi_hc3_sqr * n_integrals_2d.integrand[2];
    m1_opacities.eta_0[id_anux] = four_pi_hc3_sqr * n_integrals_2d.integrand[3];

    m1_opacities.kappa_0_a[id_nue] =
        n[id_nue] == 0. ? 0. :
                          four_pi_hc3 / (kClight * n[id_nue]) *
                              (four_pi_hc3 * n_integrals_2d.integrand[4] +
                               integrals_1d.integrand[2]);
    m1_opacities.kappa_0_a[id_anue] =
        n[id_anue] == 0. ? 0. :
                           four_pi_hc3 / (kClight * n[id_anue]) *
                               (four_pi_hc3 * n_integrals_2d.integrand[5] +
                                integrals_1d.integrand[3]);
    m1_opacities.kappa_0_a[id_nux] =
        n[id_nux] == 0. ? 0. :
                          four_pi_hc3_sqr / (kClight * n[id_nux]) *
                              n_integrals_2d.integrand[6];
    m1_opacities.kappa_0_a[id_anux] =
        n[id_anux] == 0. ? 0. :
                           four_pi_hc3_sqr / (kClight * n[id_anux]) *
                               n_integrals_2d.integrand[7];

    m1_opacities.eta[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[0] + integrals_1d.integrand[4]);
    m1_opacities.eta[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[1] + integrals_1d.integrand[5]);
    m1_opacities.eta[id_nux]  = four_pi_hc3_sqr * e_integrals_2d.integrand[2];
    m1_opacities.eta[id_anux] = four_pi_hc3_sqr * e_integrals_2d.integrand[3];

    m1_opacities.kappa_a[id_nue] =
        n[id_nue] == 0. ? 0. :
                          four_pi_hc3 / (kClight * J[id_nue]) *
                              (four_pi_hc3 * e_integrals_2d.integrand[4] +
                               integrals_1d.integrand[6]);
    m1_opacities.kappa_a[id_anue] =
        n[id_anue] == 0. ? 0. :
                           four_pi_hc3 / (kClight * J[id_anue]) *
                               (four_pi_hc3 * e_integrals_2d.integrand[5] +
                                integrals_1d.integrand[7]);
    m1_opacities.kappa_a[id_nux] = n[id_nux] == 0. ?
                                       0. :
                                       four_pi_hc3_sqr / (kClight * J[id_nux]) *
                                           e_integrals_2d.integrand[6];
    m1_opacities.kappa_a[id_anux] =
        n[id_anux] == 0. ? 0. :
                           four_pi_hc3_sqr / (kClight * J[id_anux]) *
                               e_integrals_2d.integrand[7];

    m1_opacities.kappa_s[id_nue] =
        n[id_nue] == 0. ?
            0. :
            four_pi_hc3 / (kClight * J[id_nue]) * integrals_1d.integrand[8];
    m1_opacities.kappa_s[id_anue] =
        n[id_anue] == 0. ?
            0. :
            four_pi_hc3 / (kClight * J[id_anue]) * integrals_1d.integrand[9];
    m1_opacities.kappa_s[id_nux] =
        n[id_nux] == 0. ?
            0. :
            four_pi_hc3 / (kClight * J[id_nux]) * integrals_1d.integrand[10];
    m1_opacities.kappa_s[id_anux] =
        n[id_anux] == 0. ?
            0. :
            four_pi_hc3 / (kClight * J[id_anux]) * integrals_1d.integrand[11];

    return m1_opacities;
}

/* Compute the integrands for the computation of the spectral emissivity and
 * inverse mean free path */
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand SpectralIntegrand(double* var, void* p)
{
    // energies and parameters
    double nu_bar = var[0]; // [MeV]

    GreyOpacityParams* my_grey_opacity_params = (GreyOpacityParams*)p;
    MyEOSParams my_eos_params  = my_grey_opacity_params->eos_pars;
    OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;
    OpacityParams opacity_pars = my_grey_opacity_params->opacity_pars;

    double nu = my_grey_opacity_params->kernel_pars.pair_kernel_params.omega;

    double block_factor[total_num_species]; // blocking factor

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species], g_nu_bar[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
        g_nu_bar[idx] =
            TotalNuF(nu_bar, &my_grey_opacity_params->distr_pars, idx);
    }

    // compute the pair kernels
    MyKernelOutput pair_kernels_m1 = {
        0}; //{.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.};
    if (opacity_flags.use_pair)
    {
        my_grey_opacity_params->kernel_pars.pair_kernel_params.omega_prime =
            nu_bar;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = 1.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.filter    = 0.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax      = 0;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu        = 1.;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime  = 1.;
        // pair_kernels_m1 = PairKernelsM1(&my_eos_params,
        // &my_grey_opacity_params->kernel_pars.pair_kernel_params);
        pair_kernels_m1 = PairKernelsOptimized(
            &my_eos_params,
            &my_grey_opacity_params->kernel_pars.pair_kernel_params);
    }

    // compute the bremsstrahlung kernels
    MyKernelOutput brem_kernels_m1 = {0};
    if (opacity_flags.use_brem)
    {
        my_grey_opacity_params->kernel_pars.brem_kernel_params.omega_prime =
            nu_bar;
        if (opacity_pars.use_BRT_brem == true)
        {
            brem_kernels_m1 = BremKernelsBRT06(
                &my_grey_opacity_params->kernel_pars.brem_kernel_params,
                &my_eos_params);
        }
        else
        {
            my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
            my_grey_opacity_params->kernel_pars.brem_kernel_params
                .use_NN_medium_corr =
                my_grey_opacity_params->opacity_pars.use_NN_medium_corr;
            brem_kernels_m1 = BremKernelsLegCoeff(
                &my_grey_opacity_params->kernel_pars.brem_kernel_params,
                &my_eos_params);
        }
    }

    // compute the inelastic NES/NPS kernels
    MyKernelOutput inelastic_kernels_m1 = {0};
    if (opacity_flags.use_inelastic_scatt)
    {
        my_grey_opacity_params->kernel_pars.inelastic_kernel_params
            .omega_prime     = nu_bar;
        inelastic_kernels_m1 = InelasticScattKernels(
            &my_grey_opacity_params->kernel_pars.inelastic_kernel_params,
            &my_grey_opacity_params->eos_pars);
    }

    double pro_term[total_num_species] = {0};

    if (opacity_pars.neglect_blocking == false)
    {
        for (int idx = 0; idx < total_num_species; idx++)
        {
            block_factor[idx] = 1. - g_nu_bar[idx];
        }
    }
    else
    {
        for (int idx = 0; idx < total_num_species; idx++)
        {
            block_factor[idx] = 1.;
        }
    }

    pro_term[id_nue] =
        (pair_kernels_m1.em[id_nue] + brem_kernels_m1.em[id_nue]) *
        block_factor[id_anue];
    pro_term[id_anue] =
        (pair_kernels_m1.em[id_anue] + brem_kernels_m1.em[id_anue]) *
        block_factor[id_nue];
    pro_term[id_nux] =
        (pair_kernels_m1.em[id_nux] + brem_kernels_m1.em[id_nux]) *
        block_factor[id_anux];
    pro_term[id_anux] =
        (pair_kernels_m1.em[id_anux] + brem_kernels_m1.em[id_anux]) *
        block_factor[id_nux];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        pro_term[idx] += inelastic_kernels_m1.em[idx] * g_nu_bar[idx];
    }

    double ann_term[total_num_species] = {0};

    ann_term[id_nue] =
        (pair_kernels_m1.abs[id_nue] + brem_kernels_m1.abs[id_nue]) *
        g_nu_bar[id_anue];
    ann_term[id_anue] =
        (pair_kernels_m1.abs[id_anue] + brem_kernels_m1.abs[id_anue]) *
        g_nu_bar[id_nue];
    ann_term[id_nux] =
        (pair_kernels_m1.abs[id_nux] + brem_kernels_m1.abs[id_nux]) *
        g_nu_bar[id_anux];
    ann_term[id_anux] =
        (pair_kernels_m1.abs[id_anux] + brem_kernels_m1.abs[id_anux]) *
        g_nu_bar[id_nux];

    //////////////////////////////////////////////
    ////// ONLY FOR COMPARISON WITH NULIB ////////
    //////////////////////////////////////////////
    if (kirchoff_flag)
    {
        ann_term[id_nue] +=
            pair_kernels_m1.em[id_nue] * g_nu_bar[id_anue] / g_nu[id_nue];
        ann_term[id_anue] +=
            pair_kernels_m1.em[id_anue] * g_nu_bar[id_nue] / g_nu[id_anue];
        ann_term[id_nux] +=
            pair_kernels_m1.em[id_nux] * g_nu_bar[id_anux] / g_nu[id_nux];
        ann_term[id_anux] +=
            pair_kernels_m1.em[id_anux] * g_nu_bar[id_nux] / g_nu[id_anux];
    }
    //////////////////////////////////////////////

    for (int idx = 0; idx < total_num_species; idx++)
    {
        ann_term[idx] += inelastic_kernels_m1.abs[idx] * block_factor[idx];
    }

    double integrand_1[total_num_species], integrand_2[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        integrand_1[idx] = nu_bar * nu_bar * pro_term[idx];
        integrand_2[idx] = nu_bar * nu_bar * ann_term[idx];
    }

    MyQuadratureIntegrand result = {.n = 8};

    result.integrand[0] = integrand_1[id_nue];
    result.integrand[1] = integrand_1[id_anue];
    result.integrand[2] = integrand_1[id_nux];
    result.integrand[3] = integrand_1[id_anux];
    result.integrand[4] = integrand_2[id_nue];
    result.integrand[5] = integrand_2[id_anue];
    result.integrand[6] = integrand_2[id_nux];
    result.integrand[7] = integrand_2[id_anux];

    return result;
}

/* Computes the spectral emissivity and inverse mean free path */

// Version without stimulated absorption
KOKKOS_INLINE_FUNCTION
SpectralOpacities ComputeSpectralOpacitiesNotStimulatedAbs(
    const double nu, MyQuadrature* quad_1d,
    GreyOpacityParams* my_grey_opacity_params)
{
    // compute some constants
    static const double four_pi_hc3 =
        (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]

    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;

    // set up 1d integration
    MyFunctionMultiD integrand_m1_1d;
    MyQuadratureIntegrand integrand_m1_1d_info = {.n = 8};
    integrand_m1_1d.function                   = &SpectralIntegrand;
    integrand_m1_1d.dim                        = 1;
    integrand_m1_1d.params                     = my_grey_opacity_params;
    integrand_m1_1d.my_quadrature_integrand    = integrand_m1_1d_info;

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    }

    const double eta_e = my_grey_opacity_params->eos_pars.mu_e /
                         my_grey_opacity_params->eos_pars.temp;

    double s[8];
    for (int i = 0; i < 8; i++)
    {
        // s[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
        // s[i] = 0.5 * 4.364 * my_grey_opacity_params->eos_pars.temp;
        s[i] =
            0.5 * my_grey_opacity_params->eos_pars.temp *
            (FDI_p4(eta_e) / FDI_p3(eta_e) + FDI_p4(-eta_e) / FDI_p3(-eta_e));
        // s[i] = my_grey_opacity_params->eos_pars.temp;
        // s[i] = 2.425E-03 * my_grey_opacity_params->eos_pars.temp;
    }

    MyQuadratureIntegrand integrals_1d =
        GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s);

    MyOpacity abs_em_beta = {0};
    if (my_grey_opacity_params->opacity_flags.use_abs_em)
    {
        abs_em_beta = AbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                                 &my_grey_opacity_params->eos_pars); // [s^-1]
    }

    double iso_scatt = 0.;
    if (my_grey_opacity_params->opacity_flags.use_iso)
    {
        iso_scatt = IsoScattLegCoeff(nu, &my_grey_opacity_params->opacity_pars,
                                     &my_grey_opacity_params->eos_pars, 0);
    }

    SpectralOpacities sp_opacities;

    sp_opacities.j[id_nue] =
        abs_em_beta.em[id_nue] + four_pi_hc3 * integrals_1d.integrand[0];
    sp_opacities.j[id_anue] =
        abs_em_beta.em[id_anue] + four_pi_hc3 * integrals_1d.integrand[1];
    sp_opacities.j[id_nux] =
        abs_em_beta.em[id_nux] + four_pi_hc3 * integrals_1d.integrand[2];
    sp_opacities.j[id_anux] =
        abs_em_beta.em[id_anux] + four_pi_hc3 * integrals_1d.integrand[3];

    sp_opacities.kappa[id_nue] =
        (abs_em_beta.abs[id_nue] + four_pi_hc3 * integrals_1d.integrand[4]) /
        kClight;
    sp_opacities.kappa[id_anue] =
        (abs_em_beta.abs[id_anue] + four_pi_hc3 * integrals_1d.integrand[5]) /
        kClight;
    sp_opacities.kappa[id_nux] =
        (abs_em_beta.abs[id_nux] + four_pi_hc3 * integrals_1d.integrand[6]) /
        kClight;
    sp_opacities.kappa[id_anux] =
        (abs_em_beta.abs[id_anux] + four_pi_hc3 * integrals_1d.integrand[7]) /
        kClight;

    sp_opacities.j_s[id_nue]  = 4. * kPi * nu * nu * g_nu[id_nue] * iso_scatt;
    sp_opacities.j_s[id_anue] = 4. * kPi * nu * nu * g_nu[id_anue] * iso_scatt;
    sp_opacities.j_s[id_nux]  = 4. * kPi * nu * nu * g_nu[id_nux] * iso_scatt;
    sp_opacities.j_s[id_anux] = 4. * kPi * nu * nu * g_nu[id_anux] * iso_scatt;

    sp_opacities.kappa_s[id_nue] =
        4. * kPi * nu * nu * (1. - g_nu[id_nue]) * iso_scatt / kClight;
    sp_opacities.kappa_s[id_anue] =
        4. * kPi * nu * nu * (1. - g_nu[id_anue]) * iso_scatt / kClight;
    sp_opacities.kappa_s[id_nux] =
        4. * kPi * nu * nu * (1. - g_nu[id_nux]) * iso_scatt / kClight;
    sp_opacities.kappa_s[id_anux] =
        4. * kPi * nu * nu * (1. - g_nu[id_anux]) * iso_scatt / kClight;

    return sp_opacities;
}


// Version with stimulated absorption
KOKKOS_INLINE_FUNCTION
SpectralOpacities
ComputeSpectralOpacitiesStimulatedAbs(const double nu, MyQuadrature* quad_1d,
                                      GreyOpacityParams* my_grey_opacity_params)
{
    SpectralOpacities spec_opacs = ComputeSpectralOpacitiesNotStimulatedAbs(
        nu, quad_1d, my_grey_opacity_params);

    for (int idx = 0; idx < total_num_species; idx++)
    {
        spec_opacs.kappa[idx] =
            spec_opacs.j[idx] / kClight + spec_opacs.kappa[idx];
        spec_opacs.kappa_s[idx] =
            spec_opacs.j_s[idx] / kClight + spec_opacs.kappa_s[idx];
    }

    return spec_opacs;
}

#endif // BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_