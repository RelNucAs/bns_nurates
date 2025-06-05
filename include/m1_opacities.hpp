#ifndef BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  m1_opacities.hpp
//  \brief header file for all integration routines

#include "bns_nurates.hpp"
#include "kernels.hpp"
#include "integration.hpp"
#include "distribution.hpp"

/* Thresholds on the neutrino energy number and energy density. If values are
below the thresholds, absorption opacities or scattering opacities are set to 0.
*/
constexpr BS_REAL THRESHOLD_N = 1e-21;
constexpr BS_REAL THRESHOLD_J = 1e-25;

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
void Scattering1DIntegrand(const MyQuadrature* quad,
                           GreyOpacityParams* grey_pars, const BS_REAL* t,
                           BS_REAL out[][BS_N_MAX])
{
    constexpr BS_REAL four_pi = 4 * kBS_Pi;

    const int n = quad->nx;

    BS_REAL nu, iso_scatt, aux;

    BS_REAL g_nu[total_num_species];

    for (int i = 0; i < n; ++i)
    {
        for (int idx = 0; idx < total_num_species; ++idx)
        {
            nu = t[idx] * quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[%d]=%e, "
                      "quad->points[%d]=%e)",
                      nu, idx, t[idx], i, quad->points[i]);

            // compute the neutrino distribution function
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            iso_scatt = IsoScattTotal(nu, &grey_pars->opacity_pars,
                                      &grey_pars->eos_pars);

            aux = four_pi * POW2(nu) * POW2(nu) * nu * iso_scatt;

            out[idx][i] = g_nu[idx] * aux;

            nu = t[idx] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[%d]=%e, "
                      "quad->points[%d]=%e)",
                      nu, idx, t[idx], i, quad->points[i]);

            // compute the neutrino distribution function
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            iso_scatt = IsoScattTotal(nu, &grey_pars->opacity_pars,
                                      &grey_pars->eos_pars);

            aux = four_pi * POW2(nu) * POW2(nu) * nu * iso_scatt;

            out[idx][n + i] = g_nu[idx] * aux;
        }
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void Beta1DIntegrand(const MyQuadrature* quad, GreyOpacityParams* grey_pars,
                     const BS_REAL* t, BS_REAL out_em[][BS_N_MAX],
                     BS_REAL out_ab[][BS_N_MAX], const int stim_abs)
{
    const int n = quad->nx;
    BS_REAL nu, nu_sqr, g_nu;
    MyOpacity abs_em_beta;

    if (stim_abs == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            nu = t[id_nue] * quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_nue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_nue], i, quad->points[i]);
            nu_sqr = POW2(nu);

            g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

            abs_em_beta = StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]


            out_em[id_nue][i] = nu_sqr * abs_em_beta.em[id_nue];

            // ab = em + ab (stimulated absorption)
            out_ab[id_nue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue];

            nu = t[id_anue] * quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

            abs_em_beta = StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_anue][i] = nu_sqr * abs_em_beta.em[id_anue];

            // ab = em + ab (stimulated absorption)
            out_ab[id_anue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue];

            nu = t[id_nue] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_nue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_nue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

            abs_em_beta = StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_nue][n + i] = nu_sqr * abs_em_beta.em[id_nue];

            // ab = em + ab (stimulated absorption)
            out_ab[id_nue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue];

            nu = t[id_anue] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

            abs_em_beta = StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_anue][n + i] = nu_sqr * abs_em_beta.em[id_anue];

            // ab = em + ab (stimulated absorption)
            out_ab[id_anue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue];
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            nu = t[id_nue] * quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_nue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_nue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

            abs_em_beta = AbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]


            out_em[id_nue][i] = nu_sqr * abs_em_beta.em[id_nue];
            out_ab[id_nue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue];

            nu = t[id_anue] * quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

            abs_em_beta = AbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_anue][i] = nu_sqr * abs_em_beta.em[id_anue];
            out_ab[id_anue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue];

            nu = t[id_nue] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_nue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_nue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

            abs_em_beta = AbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_nue][n + i] = nu_sqr * abs_em_beta.em[id_nue];
            out_ab[id_nue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue];

            nu = t[id_anue] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anue]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anue], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

            abs_em_beta = AbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_anue][n + i] = nu_sqr * abs_em_beta.em[id_anue];
            out_ab[id_anue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue];
        }
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void AddBetaReactionToIntegrand(int n, BS_REAL* nu_array,
                                GreyOpacityParams* grey_pars,
                                M1MatrixKokkos2D* out, const int stim_abs)
{
    BS_REAL nu;
    MyOpacity abs_em_beta;

    if (stim_abs == 1)
    {
        for (int i = 0; i < 2 * n; ++i)
        {
            nu = nu_array[i];

            abs_em_beta = StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            for (int j = 0; i < 2 * n; ++j)
            {
                out->m1_mat_em[id_nue][i][j] += abs_em_beta.em[id_nue];
                out->m1_mat_em[id_anue][i][j] += abs_em_beta.em[id_anue];

                // ab = em + ab (stimulated absorption)
                out->m1_mat_ab[id_nue][i][j] += abs_em_beta.abs[id_nue];
                out->m1_mat_ab[id_anue][i][j] += abs_em_beta.abs[id_anue];
            }
        }
    }
    else
    {
        for (int i = 0; i < 2 * n; ++i)
        {
            nu = nu_array[i];

            abs_em_beta = AbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            for (int j = 0; i < 2 * n; ++j)
            {
                out->m1_mat_em[id_nue][i][j] += abs_em_beta.em[id_nue];
                out->m1_mat_em[id_anue][i][j] += abs_em_beta.em[id_anue];

                out->m1_mat_ab[id_nue][i][j] += abs_em_beta.abs[id_nue];
                out->m1_mat_ab[id_anue][i][j] += abs_em_beta.abs[id_anue];
            }
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void AddPairKernelsToIntegrand(int n, BS_REAL* nu_array,
                               GreyOpacityParams* grey_pars,
                               M1MatrixKokkos2D* out)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    MyKernelOutput pair_1, pair_2;

    grey_pars->kernel_pars.pair_kernel_params.cos_theta = one;
    grey_pars->kernel_pars.pair_kernel_params.filter    = zero;
    grey_pars->kernel_pars.pair_kernel_params.lmax      = zero;
    grey_pars->kernel_pars.pair_kernel_params.mu        = one;
    grey_pars->kernel_pars.pair_kernel_params.mu_prime  = one;

    for (int i = 0; i < 2 * n; ++i)
    {
        grey_pars->kernel_pars.pair_kernel_params.omega       = nu_array[i];
        grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[i];

        PairKernels(&grey_pars->eos_pars,
                    &grey_pars->kernel_pars.pair_kernel_params, &pair_1,
                    &pair_2);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            out->m1_mat_em[idx][i][i] += pair_1.em[idx];
            out->m1_mat_ab[idx][i][i] += pair_1.abs[idx];
        }

        for (int j = i + 1; j < 2 * n; ++j)
        {
            grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[j];

            PairKernels(&grey_pars->eos_pars,
                        &grey_pars->kernel_pars.pair_kernel_params, &pair_1,
                        &pair_2);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] += pair_1.em[idx];
                out->m1_mat_em[idx][j][i] += pair_2.em[idx];

                out->m1_mat_ab[idx][i][j] += pair_1.abs[idx];
                out->m1_mat_ab[idx][j][i] += pair_2.abs[idx];
            }
        }
    }

    return;
}


KOKKOS_INLINE_FUNCTION
void AddBremKernelsToIntegrand(int n, BS_REAL* nu_array,
                               GreyOpacityParams* grey_pars,
                               M1MatrixKokkos2D* out)
{
    MyKernelOutput brem_ker;

    if (grey_pars->opacity_pars.use_BRT_brem == true)
    {
        for (int i = 0; i < 2 * n; ++i)
        {

            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[i];

            // compute the brem kernels
            brem_ker =
                BremKernelsBRT06(&grey_pars->kernel_pars.brem_kernel_params,
                                 &grey_pars->eos_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][i] += brem_ker.em[0];
                out->m1_mat_ab[idx][i][i] += brem_ker.abs[0];
            }

            for (int j = i + 1; j < 2 * n; ++j)
            {
                grey_pars->kernel_pars.brem_kernel_params.omega_prime =
                    nu_array[j];

                // compute the brem kernels
                brem_ker =
                    BremKernelsBRT06(&grey_pars->kernel_pars.brem_kernel_params,
                                     &grey_pars->eos_pars);

                for (int idx = 0; idx < total_num_species; ++idx)
                {
                    out->m1_mat_em[idx][i][j] += brem_ker.em[0];
                    out->m1_mat_em[idx][j][i] += brem_ker.em[0];

                    out->m1_mat_ab[idx][i][j] += brem_ker.abs[0];
                    out->m1_mat_ab[idx][j][i] += brem_ker.abs[0];
                }
            }
        }
    }
    else
    {
        grey_pars->kernel_pars.brem_kernel_params.l = 0;
        grey_pars->kernel_pars.brem_kernel_params.use_NN_medium_corr =
            grey_pars->opacity_pars.use_NN_medium_corr;

        for (int i = 0; i < 2 * n; ++i)
        {

            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[i];

            // compute the brem kernels
            brem_ker =
                BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params,
                                    &grey_pars->eos_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][i] += brem_ker.em[0];
                out->m1_mat_ab[idx][i][i] += brem_ker.abs[0];
            }

            for (int j = i + 1; j < 2 * n; ++j)
            {
                // compute the brem kernels
                grey_pars->kernel_pars.brem_kernel_params.omega_prime =
                    nu_array[j];
                brem_ker = BremKernelsLegCoeff(
                    &grey_pars->kernel_pars.brem_kernel_params,
                    &grey_pars->eos_pars);

                for (int idx = 0; idx < total_num_species; ++idx)
                {
                    out->m1_mat_em[idx][i][j] += brem_ker.em[0];
                    out->m1_mat_em[idx][j][i] += brem_ker.em[0];

                    out->m1_mat_ab[idx][i][j] += brem_ker.abs[0];
                    out->m1_mat_ab[idx][j][i] += brem_ker.abs[0];
                }
            }
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void AddInelKernelsToIntegrand(int n, BS_REAL* nu_array,
                               GreyOpacityParams* grey_pars,
                               M1MatrixKokkos2D* out)
{
    constexpr BS_REAL one = 1;

    BS_REAL nu, nu_bar;
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    BS_REAL block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];
    MyKernelOutput inel_1, inel_2;

    for (int i = 0; i < 2 * n; ++i)
    {
        nu = nu_array[i];
        BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            if (grey_pars->opacity_pars.neglect_blocking == false)
            {
                block_factor_nu[idx] = one - g_nu[idx];
            }
            else
            {
                block_factor_nu[idx] = one;
            }
        }

        // compute the pair kernels
        grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
        grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

        inel_1 = InelasticScattKernels(
            &grey_pars->kernel_pars.inelastic_kernel_params,
            &grey_pars->eos_pars);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            out->m1_mat_em[idx][i][i] += inel_1.em[idx] * g_nu[idx];
            out->m1_mat_ab[idx][i][i] += inel_1.abs[idx] * block_factor_nu[idx];
        }


        for (int j = i + 1; j < 2 * n; ++j)
        {
            nu_bar = nu_array[j];
            BS_ASSERT(nu_bar >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      nu_bar);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu_bar[idx] = one - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu_bar[idx] = one;
                }
            }

            // compute the pair kernels
            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

            inel_1 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] += inel_1.em[idx] * g_nu_bar[idx];
                out->m1_mat_em[idx][j][i] += inel_2.em[idx] * g_nu[idx];

                out->m1_mat_ab[idx][i][j] +=
                    inel_1.abs[idx] * block_factor_nu_bar[idx];
                out->m1_mat_ab[idx][j][i] +=
                    inel_2.abs[idx] * block_factor_nu[idx];
            }
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void WeightNuNuBarReactionsWithDistr(int n, BS_REAL* nu_array,
                                     GreyOpacityParams* grey_pars,
                                     M1MatrixKokkos2D* out)
{
    constexpr BS_REAL one = 1;

    BS_REAL nu, nu_bar;
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    BS_REAL block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    for (int i = 0; i < 2 * n; ++i)
    {
        nu = nu_array[i];
        BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            if (grey_pars->opacity_pars.neglect_blocking == false)
            {
                block_factor_nu[idx] = one - g_nu[idx];
            }
            else
            {
                block_factor_nu[idx] = one;
            }
        }

        out->m1_mat_em[id_nue][i][i] *= block_factor_nu[id_anue];
        out->m1_mat_em[id_anue][i][i] *= block_factor_nu[id_nue];
        out->m1_mat_em[id_nux][i][i] *= block_factor_nu[id_anux];
        out->m1_mat_em[id_anux][i][i] *= block_factor_nu[id_nux];

        out->m1_mat_ab[id_nue][i][i] *= g_nu[id_anue];
        out->m1_mat_ab[id_anue][i][i] *= g_nu[id_nue];
        out->m1_mat_ab[id_nux][i][i] *= g_nu[id_anux];
        out->m1_mat_ab[id_anux][i][i] *= g_nu[id_nux];

        for (int j = i + 1; j < 2 * n; ++j)
        {
            nu_bar = nu_array[j];
            BS_ASSERT(nu_bar >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      nu_bar);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu_bar[idx] = one - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu_bar[idx] = one;
                }
            }

            out->m1_mat_em[id_nue][i][j] *= block_factor_nu_bar[id_anue];
            out->m1_mat_em[id_anue][i][j] *= block_factor_nu_bar[id_nue];
            out->m1_mat_em[id_nux][i][j] *= block_factor_nu_bar[id_anux];
            out->m1_mat_em[id_anux][i][j] *= block_factor_nu_bar[id_nux];

            out->m1_mat_em[id_nue][j][i] *= block_factor_nu[id_anue];
            out->m1_mat_em[id_anue][j][i] *= block_factor_nu[id_nue];
            out->m1_mat_em[id_nux][j][i] *= block_factor_nu[id_anux];
            out->m1_mat_em[id_anux][j][i] *= block_factor_nu[id_nux];

            out->m1_mat_ab[id_nue][i][j] *= g_nu_bar[id_anue];
            out->m1_mat_ab[id_anue][i][j] *= g_nu_bar[id_nue];
            out->m1_mat_ab[id_nux][i][j] *= g_nu_bar[id_anux];
            out->m1_mat_ab[id_anux][i][j] *= g_nu_bar[id_nux];

            out->m1_mat_ab[id_nue][j][i] *= g_nu[id_anue];
            out->m1_mat_ab[id_anue][j][i] *= g_nu[id_nue];
            out->m1_mat_ab[id_nux][j][i] *= g_nu[id_anux];
            out->m1_mat_ab[id_anux][j][i] *= g_nu[id_nux];
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void AddCommonWeightsToIntegrand(int n, BS_REAL* nu_array,
                                 GreyOpacityParams* grey_pars,
                                 M1MatrixKokkos2D* out, int stim_abs)
{
    constexpr BS_REAL one = 1;

    BS_REAL nu, nu_bar, nu_squared, nu_fourth;
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];

    BS_ASSERT((stim_abs == 0) || (stim_abs == 1));

    if (stim_abs == 1)
    {
        for (int i = 0; i < 2 * n; ++i)
        {
            nu = nu_array[i];
            BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
            nu_squared = POW2(nu);
            nu_fourth  = POW2(nu_squared);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

                out->m1_mat_ab[idx][i][i] =
                    nu_fourth * g_nu[idx] *
                    (out->m1_mat_em[idx][i][i] + out->m1_mat_ab[idx][i][i]);
                out->m1_mat_em[idx][i][i] *= nu_fourth;
            }

            for (int j = i + 1; j < 2 * n; ++j)
            {
                nu_bar = nu_array[j];
                BS_ASSERT(nu_bar >= 0,
                          "Neutrino energy is negative (nu_bar=%e)", nu_bar);
                nu_fourth = nu_squared * POW2(nu_bar);

                for (int idx = 0; idx < total_num_species; ++idx)
                {
                    g_nu_bar[idx] =
                        TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                    out->m1_mat_ab[idx][i][j] =
                        nu_fourth * g_nu[idx] *
                        (out->m1_mat_em[idx][i][j] + out->m1_mat_ab[idx][i][j]);
                    out->m1_mat_ab[idx][j][i] =
                        nu_fourth * g_nu_bar[idx] *
                        (out->m1_mat_em[idx][j][i] + out->m1_mat_ab[idx][j][i]);

                    out->m1_mat_em[idx][i][j] *= nu_fourth;
                    out->m1_mat_em[idx][j][i] *= nu_fourth;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < 2 * n; ++i)
        {
            nu = nu_array[i];
            BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
            nu_squared = POW2(nu);
            nu_fourth  = POW2(nu_squared);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

                out->m1_mat_ab[idx][i][i] *= nu_fourth * g_nu[idx];
                out->m1_mat_em[idx][i][i] *= nu_fourth * (one - g_nu[idx]);
            }

            for (int j = i + 1; j < 2 * n; ++j)
            {
                nu_bar = nu_array[j];
                BS_ASSERT(nu_bar >= 0,
                          "Neutrino energy is negative (nu_bar=%e)", nu_bar);
                nu_fourth = nu_squared * POW2(nu_bar);

                for (int idx = 0; idx < total_num_species; ++idx)
                {
                    g_nu_bar[idx] =
                        TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                    out->m1_mat_ab[idx][i][j] *= nu_fourth * g_nu[idx];
                    out->m1_mat_ab[idx][j][i] *= nu_fourth * g_nu_bar[idx];

                    out->m1_mat_em[idx][i][j] *= nu_fourth * (one - g_nu[idx]);
                    out->m1_mat_em[idx][j][i] *=
                        nu_fourth * (one - g_nu_bar[idx]);
                }
            }
        }
    }
    return;
}

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
 * Note that there are no BS_REAL integrals for the computation of the
 * scattering coefficient.
 */
KOKKOS_INLINE_FUNCTION
M1MatrixKokkos2D ComputeDoubleIntegrand(const MyQuadrature* quad, BS_REAL t,
                                        GreyOpacityParams* grey_pars,
                                        const int stim_abs)
{
    const int n = quad->nx;
    BS_REAL nu, nu_bar;
    BS_REAL nu_array[BS_N_MAX];
    M1MatrixKokkos2D out = {0};

    for (int i = 0; i < n; ++i)
    {
        nu_array[i]     = t * quad->points[i];
        nu_array[n + i] = t / quad->points[i];
    }

    // compute the neutrino & anti-neutrino distribution function
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    BS_REAL block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    if (grey_pars->opacity_flags.use_pair == 1)
    {
        AddPairKernelsToIntegrand(n, nu_array, grey_pars, &out);
    }


    if (grey_pars->opacity_flags.use_brem == 1)
    {
        AddBremKernelsToIntegrand(n, nu_array, grey_pars, &out);
    }

    if ((grey_pars->opacity_flags.use_pair == 1) ||
        (grey_pars->opacity_flags.use_brem == 1))
    {
        WeightNuNuBarReactionsWithDistr(n, nu_array, grey_pars, &out);
    }

    // if (grey_pars->opacity_flags.use_inelastic_scatt == 1)
    // {
    //     AddInelKernelsToIntegrand(n, nu_array, grey_pars, &out);
    // }

    // /*
    // if (grey_pars->opacity_flags.use_abs_em == 1)
    // {
    //   AddBetaReactionToIntegrand(n, nu_array, grey_pars, &out, stim_abs);
    // }
    // */

    // /*
    // //////////////////////////////////////////////
    // ////// ONLY FOR COMPARISON WITH NULIB ////////
    // //////////////////////////////////////////////
    // if (kirchoff_flag)
    // {
    //     ann_term_ij[id_nue] +=
    //         (pair.m1_mat_em[id_nue][i][j] + brem.m1_mat_em[0][i][j]) *
    //         g_nu_bar[id_anue] / g_nu[id_nue];
    //     ann_term_ij[id_anue] +=
    //         (pair.m1_mat_em[id_anue][i][j] + brem.m1_mat_em[0][i][j]) *
    //         g_nu_bar[id_nue] / g_nu[id_anue];
    //     ann_term_ij[id_nux] +=
    //         (pair.m1_mat_em[id_nux][i][j] + brem.m1_mat_em[0][i][j]) *
    //         g_nu_bar[id_anux] / g_nu[id_nux];
    //     ann_term_ij[id_anux] +=
    //         (pair.m1_mat_em[id_anux][i][j] + brem.m1_mat_em[0][i][j]) *
    //         g_nu_bar[id_nux] / g_nu[id_anux];

    //     ann_term_ji[id_nue] +=
    //         (pair.m1_mat_em[id_nue][j][i] + brem.m1_mat_em[0][j][i]) *
    //         g_nu[id_anue] / g_nu_bar[id_nue];
    //     ann_term_ji[id_anue] +=
    //         (pair.m1_mat_em[id_anue][j][i] + brem.m1_mat_em[0][j][i]) *
    //         g_nu[id_nue] / g_nu_bar[id_anue];
    //     ann_term_ji[id_nux] +=
    //         (pair.m1_mat_em[id_nux][j][i] + brem.m1_mat_em[0][j][i]) *
    //         g_nu[id_anux] / g_nu_bar[id_nux];
    //     ann_term_ji[id_anux] +=
    //         (pair.m1_mat_em[id_anux][j][i] + brem.m1_mat_em[0][j][i]) *
    //         g_nu[id_nux] / g_nu_bar[id_anux];
    // }
    // */

    // //////////////////////////////////////////////

    AddCommonWeightsToIntegrand(n, nu_array, grey_pars, &out, stim_abs);

    return out;
}

KOKKOS_INLINE_FUNCTION
M1MatrixKokkos2D ComputeNEPSIntegrand(const MyQuadrature* quad, BS_REAL t,
                                      GreyOpacityParams* grey_pars,
                                      const int stim_abs)
{
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL one  = 1;

    BS_ASSERT((stim_abs == 0) || (stim_abs == 1));

    const int n = quad->nx;
    BS_REAL x_i, x_j;
    BS_REAL nu, nu_bar, nu_fourth;

    BS_REAL tmp_em_1, tmp_em_2;
    BS_REAL tmp_abs_1, tmp_abs_2;

    // compute the neutrino & anti-neutrino distribution function
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    BS_REAL block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    MyKernelOutput inel_1, inel_2;

    M1MatrixKokkos2D out = {0};

    // for (int i = 0; i < n; ++i)
    // {
    //     x = quad->points[i];

    //     nu_array_1[i]     = t * x;
    //     nu_array_1[n + i] = t / x;

    //     nu_array_2[i] = x / (1. - x) - (1. - x) / x;
    // }

    for (int i = 0; i < n; ++i)
    {

        x_i = quad->points[i];

        for (int j = 0; j < n; ++j)
        {

            x_j = quad->points[j];

            nu = half * t * x_i * (one - x_j);
            BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
            nu_bar = half * t * x_i * (one + x_j);
            BS_ASSERT(nu_bar >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      nu_bar);

            nu_fourth = POW2(nu) * POW2(nu_bar);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx]     = one - g_nu[idx];
                    block_factor_nu_bar[idx] = one - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu[idx]     = one;
                    block_factor_nu_bar[idx] = one;
                }
            }

            // compute the pair kernels
            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

            inel_1 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                tmp_em_1  = inel_1.em[idx] * g_nu_bar[idx];
                tmp_abs_1 = inel_1.abs[idx] * block_factor_nu_bar[idx];

                tmp_em_2  = inel_2.em[idx] * g_nu[idx];
                tmp_abs_2 = inel_2.abs[idx] * block_factor_nu[idx];

                if (stim_abs == 1)
                {
                    out.m1_mat_ab[idx][i][j] =
                        nu_fourth * g_nu[idx] * (tmp_em_1 + tmp_abs_1);
                    out.m1_mat_em[idx][i][j] = nu_fourth * tmp_em_1;

                    out.m1_mat_ab[idx][i][n + j] =
                        nu_fourth * g_nu_bar[idx] * (tmp_em_2 + tmp_abs_2);
                    out.m1_mat_em[idx][i][n + j] = nu_fourth * tmp_em_2;
                }
                else
                {
                    out.m1_mat_ab[idx][i][j] =
                        nu_fourth * g_nu[idx] * tmp_abs_1;
                    out.m1_mat_em[idx][i][j] =
                        nu_fourth * (one - g_nu[idx]) * tmp_em_1;

                    out.m1_mat_ab[idx][i][n + j] =
                        nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                    out.m1_mat_em[idx][i][n + j] =
                        nu_fourth * (one - g_nu_bar[idx]) * tmp_em_2;
                }
            }

            nu = half * t * (one - x_j) / x_i;
            BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
            nu_bar = half * t * (one + x_j) / x_i;
            BS_ASSERT(nu_bar >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      nu_bar);

            nu_fourth = POW2(nu) * POW2(nu_bar);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx]     = one - g_nu[idx];
                    block_factor_nu_bar[idx] = one - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu[idx]     = one;
                    block_factor_nu_bar[idx] = one;
                }
            }

            // compute the pair kernels
            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

            inel_1 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticScattKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                tmp_em_1  = inel_1.em[idx] * g_nu_bar[idx];
                tmp_abs_1 = inel_1.abs[idx] * block_factor_nu_bar[idx];

                tmp_em_2  = inel_2.em[idx] * g_nu[idx];
                tmp_abs_2 = inel_2.abs[idx] * block_factor_nu[idx];

                if (stim_abs == 1)
                {
                    out.m1_mat_ab[idx][n + i][j] =
                        nu_fourth * g_nu[idx] * (tmp_em_1 + tmp_abs_1);
                    out.m1_mat_em[idx][n + i][j] = nu_fourth * tmp_em_1;

                    out.m1_mat_ab[idx][n + i][n + j] =
                        nu_fourth * g_nu_bar[idx] * (tmp_em_2 + tmp_abs_2);
                    out.m1_mat_em[idx][n + i][n + j] = nu_fourth * tmp_em_2;
                }
                else
                {
                    out.m1_mat_ab[idx][n + i][j] =
                        nu_fourth * g_nu[idx] * tmp_abs_1;
                    out.m1_mat_em[idx][n + i][j] =
                        nu_fourth * (one - g_nu[idx]) * tmp_em_1;

                    out.m1_mat_ab[idx][n + i][n + j] =
                        nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                    out.m1_mat_em[idx][n + i][n + j] =
                        nu_fourth * (one - g_nu_bar[idx]) * tmp_em_2;
                }
            }
        }
    }

    return out;
}

/* Computes the opacities for the M1 code
 *
 */
KOKKOS_INLINE_FUNCTION
M1Opacities ComputeM1OpacitiesGenericFormalism(
    const MyQuadrature* quad_1d, const MyQuadrature* quad_2d,
    GreyOpacityParams* my_grey_opacity_params, const int stim_abs)
{
    constexpr BS_REAL four    = 4;
    constexpr BS_REAL c_light = kBS_Clight;

    BS_REAL n[total_num_species];
    BS_REAL J[total_num_species];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        // m1_pars.n and m1_pars.J are assumed to be parsed in [nm^-3]
        // and [g nm^-1 s^-2]
        n[idx] = my_grey_opacity_params->m1_pars.n[idx];
        J[idx] = my_grey_opacity_params->m1_pars.J[idx];

        // g nm-1 s^-1 -> MeV nm-3 conversion
        J[idx] = J[idx] / kBS_MeV;
    }

    const BS_REAL temp  = my_grey_opacity_params->eos_pars.temp;
    const BS_REAL eta_e = my_grey_opacity_params->eos_pars.mu_e / temp;

    constexpr BS_REAL three_halves  = 1.5;
    constexpr BS_REAL five_sixths   = 0.8333333333333333;
    constexpr BS_REAL five          = 5;
    constexpr BS_REAL temp_multiple = 0.5 * 4.364;
    const BS_REAL s_pair            = temp_multiple * temp;
    // const BS_REAL s_pair              = temp * (FDI_p4(eta_e) / FDI_p3(eta_e)
    // + FDI_p4(-eta_e) / FDI_p3(-eta_e));
    const BS_REAL s_nux  = three_halves * temp;
    const BS_REAL s_neps = temp_multiple * temp;

    BS_REAL s_beta[total_num_species] = {0}, s_iso[total_num_species] = {0};

    if (eta_e > -30. and eta_e < 30.)
    {
        s_beta[id_nue]  = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
        s_beta[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    }
    else if (eta_e > 30.)
    {
        s_beta[id_nue]  = temp * eta_e * five_sixths;
        s_beta[id_anue] = temp * five;
    }
    else
    {
        s_beta[id_nue]  = temp * five;
        s_beta[id_anue] = temp * eta_e * five_sixths;
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        s_iso[idx] = (n[idx] > THRESHOLD_N) ? (J[idx] / n[idx]) : s_nux;
    }



    MyQuadratureIntegrand iso_integrals = {0};
    if (my_grey_opacity_params->opacity_flags.use_iso == 1)
    {
        BS_REAL out_iso[total_num_species][BS_N_MAX];
        Scattering1DIntegrand(quad_1d, my_grey_opacity_params, s_iso, out_iso);
        iso_integrals = GaussLegendreIntegrate1DMatrix(
            quad_1d, total_num_species, out_iso, s_iso);
    }

    MyQuadratureIntegrand beta_n_em_integrals  = {0};
    MyQuadratureIntegrand beta_j_em_integrals  = {0};
    MyQuadratureIntegrand beta_n_abs_integrals = {0};
    MyQuadratureIntegrand beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_abs_em == 1)
    {
        BS_REAL out_beta_em[total_num_species][BS_N_MAX];
        BS_REAL out_beta_ab[total_num_species][BS_N_MAX];

        Beta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta, out_beta_em,
                        out_beta_ab, stim_abs);

        GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_beta_em,
                                                 s_beta, &beta_n_em_integrals,
                                                 &beta_j_em_integrals);
        GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_beta_ab,
                                                 s_beta, &beta_n_abs_integrals,
                                                 &beta_j_abs_integrals);
    }


    MyQuadratureIntegrand n_integrals_2d = {0};
    MyQuadratureIntegrand e_integrals_2d = {0};
    
    if ((my_grey_opacity_params->opacity_flags.use_pair == 1) ||
        (my_grey_opacity_params->opacity_flags.use_brem == 1))
    {
        M1MatrixKokkos2D out_pair =
           ComputeDoubleIntegrand(quad_2d, s_pair, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForM1Coeffs(quad_2d, &out_pair, s_pair,
                                                  &n_integrals_2d, &e_integrals_2d);
    }

    MyQuadratureIntegrand n_neps_2d = {0};
    MyQuadratureIntegrand e_neps_2d = {0};
   
    if (my_grey_opacity_params->opacity_flags.use_inelastic_scatt == 1)
    {
        M1MatrixKokkos2D out_inel = ComputeNEPSIntegrand(
            quad_2d, four * s_neps, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForNEPS(quad_2d, &out_inel, four * s_neps,
                                              &n_neps_2d, &e_neps_2d);
    }

    M1Opacities m1_opacities = {0};

    /* Set all opacities to zero. They'll be left as 0 if the neutrino
    number/energy density is too low (to avoid inf/nan, since the number/energy
    density appears in the denominator). Note that emissivities do no need this
    precaution (and it also make sense physically: you can produce neutrinos
    even if there are none to start with). */
    constexpr BS_REAL zero = 0;
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        m1_opacities.kappa_0_a[idx] = zero;
        m1_opacities.kappa_a[idx]   = zero;
        m1_opacities.kappa_s[idx]   = zero;
    }

    /* Electron neutrinos */
    m1_opacities.eta_0[id_nue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[0] +
                                            n_neps_2d.integrand[0]) +
                          beta_n_em_integrals.integrand[id_nue]);
    m1_opacities.eta[id_nue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[0] +
                                            e_neps_2d.integrand[0]) +
                          beta_j_em_integrals.integrand[id_nue]);
    if (n[id_nue] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_nue] =
            kBS_FourPi_hc3 / (c_light * n[id_nue]) *
            (kBS_FourPi_hc3 *
                 (n_integrals_2d.integrand[4] + n_neps_2d.integrand[4]) +
             beta_n_abs_integrals.integrand[id_nue]);
    }
    if (J[id_nue] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_nue] =
            n[id_nue] == zero ?
                zero :
                kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                    (kBS_FourPi_hc3 * (e_integrals_2d.integrand[4] +
                                       e_neps_2d.integrand[4]) +
                     beta_j_abs_integrals.integrand[id_nue]);
        m1_opacities.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                       iso_integrals.integrand[id_nue];
    }

    /* Electron anti-neutrinos */
    m1_opacities.eta_0[id_anue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[1] +
                                            n_neps_2d.integrand[1]) +
                          beta_n_em_integrals.integrand[id_anue]);
    m1_opacities.eta[id_anue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[1] +
                                            e_neps_2d.integrand[1]) +
                          beta_j_em_integrals.integrand[id_anue]);
    if (n[id_anue] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_anue] =
            kBS_FourPi_hc3 / (c_light * n[id_anue]) *
            (kBS_FourPi_hc3 *
                 (n_integrals_2d.integrand[5] + n_neps_2d.integrand[5]) +
             beta_n_abs_integrals.integrand[id_anue]);
    }
    if (J[id_anue] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_anue] =
            kBS_FourPi_hc3 / (c_light * J[id_anue]) *
            (kBS_FourPi_hc3 *
                 (e_integrals_2d.integrand[5] + e_neps_2d.integrand[5]) +
             beta_j_abs_integrals.integrand[id_anue]);
        m1_opacities.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                        (c_light * J[id_anue]) *
                                        iso_integrals.integrand[id_anue];
    }

    /* Heavy neutrinos */
    m1_opacities.eta_0[id_nux] =
        kBS_FourPi_hc3_sqr *
        (n_integrals_2d.integrand[2] + n_neps_2d.integrand[2]);
    m1_opacities.eta[id_nux] =
        kBS_FourPi_hc3_sqr *
        (e_integrals_2d.integrand[2] + e_neps_2d.integrand[2]);
    if (n[id_nux] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_nux] =
            kBS_FourPi_hc3_sqr / (c_light * n[id_nux]) *
            (n_integrals_2d.integrand[6] + n_neps_2d.integrand[6]);
    }
    if (J[id_nux] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_nux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
            (e_integrals_2d.integrand[6] + e_neps_2d.integrand[6]);
        m1_opacities.kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J[id_nux]) *
                                       iso_integrals.integrand[id_nux];
    }

    /* Heavy anti-neutrinos */
    m1_opacities.eta_0[id_anux] =
        kBS_FourPi_hc3_sqr *
        (n_integrals_2d.integrand[3] + n_neps_2d.integrand[3]);
    m1_opacities.eta[id_anux] =
        kBS_FourPi_hc3_sqr *
        (e_integrals_2d.integrand[3] + e_neps_2d.integrand[3]);
    if (n[id_anux] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_anux] =
            n[id_anux] == zero ?
                zero :
                kBS_FourPi_hc3_sqr / (c_light * n[id_anux]) *
                    (n_integrals_2d.integrand[7] + n_neps_2d.integrand[7]);
    }
    if (J[id_anux] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_anux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
            (e_integrals_2d.integrand[7] + e_neps_2d.integrand[7]);

        m1_opacities.kappa_s[id_anux] = kBS_FourPi_hc3 /
                                        (c_light * J[id_anux]) *
                                        iso_integrals.integrand[id_anux];
    }

    return m1_opacities;
}

KOKKOS_INLINE_FUNCTION
M1Opacities
ComputeM1OpacitiesNotStimulated(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params)
{
    return ComputeM1OpacitiesGenericFormalism(quad_1d, quad_2d,
                                              my_grey_opacity_params, 0);
}

KOKKOS_INLINE_FUNCTION
M1Opacities ComputeM1Opacities(const MyQuadrature* quad_1d,
                               const MyQuadrature* quad_2d,
                               GreyOpacityParams* my_grey_opacity_params)
{
    return ComputeM1OpacitiesGenericFormalism(quad_1d, quad_2d,
                                              my_grey_opacity_params, 1);
}

/* Compute the integrands for the computation of the spectral emissivity and
 * inverse mean free path */
KOKKOS_INLINE_FUNCTION
MyQuadratureIntegrand SpectralIntegrand(BS_REAL* var, void* p)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    // energies and parameters
    BS_REAL nu_bar = var[0]; // [MeV]

    GreyOpacityParams* my_grey_opacity_params = (GreyOpacityParams*)p;
    MyEOSParams my_eos_params  = my_grey_opacity_params->eos_pars;
    OpacityFlags opacity_flags = my_grey_opacity_params->opacity_flags;
    OpacityParams opacity_pars = my_grey_opacity_params->opacity_pars;

    BS_REAL nu = my_grey_opacity_params->kernel_pars.pair_kernel_params.omega;

    BS_REAL block_factor[total_num_species]; // blocking factor

    // compute the neutrino & anti-neutrino distribution function
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];

    for (int idx = 0; idx < total_num_species; ++idx)
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
        my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = one;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.filter    = zero;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax      = zero;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu        = one;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime  = one;
        // pair_kernels_m1 = PairKernels(&my_eos_params,
        // &my_grey_opacity_params->kernel_pars.pair_kernel_params);
        pair_kernels_m1 = PairKernels(
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

    BS_REAL pro_term[total_num_species] = {0};

    if (opacity_pars.neglect_blocking == false)
    {
        for (int idx = 0; idx < total_num_species; ++idx)
        {
            block_factor[idx] = one - g_nu_bar[idx];
        }
    }
    else
    {
        for (int idx = 0; idx < total_num_species; ++idx)
        {
            block_factor[idx] = one;
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

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        pro_term[idx] += inelastic_kernels_m1.em[idx] * g_nu_bar[idx];
    }

    BS_REAL ann_term[total_num_species] = {0};

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

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        ann_term[idx] += inelastic_kernels_m1.abs[idx] * block_factor[idx];
    }

    BS_REAL integrand_1[total_num_species], integrand_2[total_num_species];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        integrand_1[idx] = POW2(nu_bar) * pro_term[idx];
        integrand_2[idx] = POW2(nu_bar) * ann_term[idx];
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
    const BS_REAL nu, MyQuadrature* quad_1d,
    GreyOpacityParams* my_grey_opacity_params)
{
    constexpr BS_REAL zero    = 0;
    constexpr BS_REAL one     = 1;
    constexpr BS_REAL four    = 4;
    constexpr BS_REAL four_pi = 4 * kBS_Pi;
    constexpr BS_REAL c_light = kBS_Clight;

    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;

    GreyOpacityParams local_grey_params = *my_grey_opacity_params;
    local_grey_params.opacity_flags.use_inelastic_scatt = 0;

    // set up 1d integration
    MyFunctionMultiD integrand_m1_1d;
    MyQuadratureIntegrand integrand_m1_1d_info = {.n = 8};
    integrand_m1_1d.function                   = &SpectralIntegrand;
    integrand_m1_1d.dim                        = 1;
    integrand_m1_1d.params                     = &local_grey_params;
    integrand_m1_1d.my_quadrature_integrand    = integrand_m1_1d_info;

    // compute the neutrino & anti-neutrino distribution function
    BS_REAL g_nu[total_num_species];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        g_nu[idx] = TotalNuF(nu, &my_grey_opacity_params->distr_pars, idx);
    }

    const BS_REAL eta_e = my_grey_opacity_params->eos_pars.mu_e /
                          my_grey_opacity_params->eos_pars.temp;

    constexpr BS_REAL temp_multiple = 0.5 * 4.364;

    BS_REAL s_pair[8], s_neps[8];

    for (int i = 0; i < 8; ++i)
    {
        // s[i] = 1.5 * my_grey_opacity_params->eos_pars.temp;
        s_pair[i] = temp_multiple * my_grey_opacity_params->eos_pars.temp;
        s_neps[i] = nu;
        // s[i] = my_grey_opacity_params->eos_pars.temp;
        // s[i] = 2.425E-03 * my_grey_opacity_params->eos_pars.temp;
        // s[i] =
        // 0.5 * my_grey_opacity_params->eos_pars.temp *
        // (FDI_p4(eta_e) / FDI_p3(eta_e) + FDI_p4(-eta_e) / FDI_p3(-eta_e));
    }

    MyQuadratureIntegrand integrals_pair_1d =
        GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_pair);

    MyQuadratureIntegrand integrals_neps_1d = {0};
    if (my_grey_opacity_params->opacity_flags.use_inelastic_scatt == 1)
    {
        local_grey_params.opacity_flags                     = {0};
        local_grey_params.opacity_flags.use_inelastic_scatt = 1;
        integrand_m1_1d.params = &local_grey_params;
        integrals_neps_1d =
            GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_neps);
    }

    MyOpacity abs_em_beta = {0};
    if (my_grey_opacity_params->opacity_flags.use_abs_em)
    {
        abs_em_beta = AbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                                 &my_grey_opacity_params->eos_pars); // [s^-1]
    }

    BS_REAL iso_scatt = zero;
    if (my_grey_opacity_params->opacity_flags.use_iso)
    {
        iso_scatt = IsoScattLegCoeff(nu, &my_grey_opacity_params->opacity_pars,
                                     &my_grey_opacity_params->eos_pars, 0);
    }

    
    SpectralOpacities sp_opacities;

    sp_opacities.j[id_nue] = abs_em_beta.em[id_nue] +
                             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[0] +
                                               integrals_neps_1d.integrand[0]);
    sp_opacities.j[id_anue] = abs_em_beta.em[id_anue] +
                              kBS_FourPi_hc3 * (integrals_pair_1d.integrand[1] +
                                                integrals_neps_1d.integrand[1]);
    sp_opacities.j[id_nux] = abs_em_beta.em[id_nux] +
                             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[2] +
                                               integrals_neps_1d.integrand[2]);
    sp_opacities.j[id_anux] = abs_em_beta.em[id_anux] +
                              kBS_FourPi_hc3 * (integrals_pair_1d.integrand[3] +
                                                integrals_neps_1d.integrand[3]);

    sp_opacities.kappa[id_nue] =
        (abs_em_beta.abs[id_nue] +
         kBS_FourPi_hc3 * (integrals_pair_1d.integrand[4] +
                           integrals_neps_1d.integrand[4])) /
        c_light;
    sp_opacities.kappa[id_anue] =
        (abs_em_beta.abs[id_anue] +
         kBS_FourPi_hc3 * (integrals_pair_1d.integrand[5] +
                           integrals_neps_1d.integrand[5])) /
        c_light;
    sp_opacities.kappa[id_nux] =
        (abs_em_beta.abs[id_nux] +
         kBS_FourPi_hc3 * (integrals_pair_1d.integrand[6] +
                           integrals_neps_1d.integrand[6])) /
        c_light;
    sp_opacities.kappa[id_anux] =
        (abs_em_beta.abs[id_anux] +
         kBS_FourPi_hc3 * (integrals_pair_1d.integrand[7] +
                           integrals_neps_1d.integrand[7])) /
        c_light;

    sp_opacities.j_s[id_nue]  = four_pi * POW2(nu) * g_nu[id_nue] * iso_scatt;
    sp_opacities.j_s[id_anue] = four_pi * POW2(nu) * g_nu[id_anue] * iso_scatt;
    sp_opacities.j_s[id_nux]  = four_pi * POW2(nu) * g_nu[id_nux] * iso_scatt;
    sp_opacities.j_s[id_anux] = four_pi * POW2(nu) * g_nu[id_anux] * iso_scatt;

    sp_opacities.kappa_s[id_nue] =
        four_pi * POW2(nu) * (one - g_nu[id_nue]) * iso_scatt / c_light;
    sp_opacities.kappa_s[id_anue] =
        four_pi * POW2(nu) * (one - g_nu[id_anue]) * iso_scatt / c_light;
    sp_opacities.kappa_s[id_nux] =
        four_pi * POW2(nu) * (one - g_nu[id_nux]) * iso_scatt / c_light;
    sp_opacities.kappa_s[id_anux] =
        four_pi * POW2(nu) * (one - g_nu[id_anux]) * iso_scatt / c_light;

    return sp_opacities;
}


// Version with stimulated absorption
KOKKOS_INLINE_FUNCTION
SpectralOpacities
ComputeSpectralOpacitiesStimulatedAbs(const BS_REAL nu, MyQuadrature* quad_1d,
                                      GreyOpacityParams* my_grey_opacity_params)
{
    constexpr BS_REAL c_light = kBS_Clight;

    SpectralOpacities spec_opacs = ComputeSpectralOpacitiesNotStimulatedAbs(
        nu, quad_1d, my_grey_opacity_params);

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        spec_opacs.kappa[idx] =
            spec_opacs.j[idx] / c_light + spec_opacs.kappa[idx];
        spec_opacs.kappa_s[idx] =
            spec_opacs.j_s[idx] / c_light + spec_opacs.kappa_s[idx];
    }

    return spec_opacs;
}

#endif // BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_
