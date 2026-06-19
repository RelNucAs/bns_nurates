#ifndef BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_HPP_
#define BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_HPP_

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

    if (grey_pars->opacity_pars.brem_implementation == BREM_BRT06)
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
    else if (grey_pars->opacity_pars.brem_implementation == BREM_HR98)
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
    else if (grey_pars->opacity_pars.brem_implementation == BREM_GP19)
    {
        for (int i = 0; i < 2 * n; ++i)
        {

            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[i];

            // compute the brem kernels
            brem_ker =
                BremKernelAbsGP19(&grey_pars->kernel_pars.brem_kernel_params,
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
                brem_ker = BremKernelAbsGP19(
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
    else
    {
        BS_ASSERT(false, "Unknown bremsstrahlung implementation: %d",
                  (int)grey_pars->opacity_pars.brem_implementation);
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
    BS_REAL nu_array[BS_N_MAX];
    M1MatrixKokkos2D out = {0};

    for (int i = 0; i < n; ++i)
    {
        nu_array[i]     = t * quad->points[i];
        nu_array[n + i] = t / quad->points[i];
    }

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
    // compute the neutrino & anti-neutrino distribution function
    // BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    // BS_REAL block_factor_nu[total_num_species],
    //    block_factor_nu_bar[total_num_species];
    //
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


// ===========================================================================
// FUSED 2D integrals (no materialized M1MatrixKokkos2D)
//
// The original code path computes a full [species][2n][2n] matrix
// (M1MatrixKokkos2D, 25.6 KB) on the stack via ComputeDoubleIntegrand /
// ComputeNEPSIntegrand, then reduces it with GaussLegendreIntegrate2DMatrix*.
// With one GPU thread per cell that 25.6 KB (x2: pair + NEPS) lives in local
// memory -> register spill -> near-zero occupancy.
//
// The functions below FUSE the fill and the reduction: every matrix entry is
// weighted and immediately accumulated into the ~16 running integrand values,
// so the big matrix is never stored.  Results are identical (to round-off) to
// ComputeDoubleIntegrand + GaussLegendreIntegrate2DMatrixForM1Coeffs and
// ComputeNEPSIntegrand + GaussLegendreIntegrate2DMatrixForNEPS, but each cell's
// working set is a handful of registers, restoring full occupancy while keeping
// the simple one-thread-per-cell parallelism.
// ===========================================================================

// Fused pair+brem double integral.  Equivalent to:
//   M1MatrixKokkos2D m = ComputeDoubleIntegrand(quad, t, gop, stim_abs);
//   GaussLegendreIntegrate2DMatrixForM1Coeffs(quad, &m, t, n2d, e2d);
// (accumulates into *n2d, *e2d which the caller initialised to zero).
KOKKOS_INLINE_FUNCTION
void AddPairBremDoubleIntegralFused(const MyQuadrature* quad, BS_REAL t,
                                    GreyOpacityParams* gop, const int stim_abs,
                                    MyQuadratureIntegrand* n2d,
                                    MyQuadratureIntegrand* e2d)
{
    constexpr BS_REAL one  = 1;
    const int n  = quad->nx;
    const int N2 = 2 * n;

    const bool use_pair         = (gop->opacity_flags.use_pair == 1);
    const bool use_brem         = (gop->opacity_flags.use_brem == 1);
    const bool neglect_blocking = gop->opacity_pars.neglect_blocking;
    const BremImpl brem_impl    = gop->opacity_pars.brem_implementation;
    // blocking partner: nue<->anue, nux<->anux
    const int pair_partner[total_num_species] = {id_anue, id_nue, id_anux, id_nux};

    BS_REAL acc_n_em[total_num_species] = {0}, acc_n_ab[total_num_species] = {0};
    BS_REAL acc_e_em[total_num_species] = {0}, acc_e_ab[total_num_species] = {0};

    for (int a = 0; a < N2; ++a)
    {
        const int    ia    = (a < n) ? a : a - n;
        const BS_REAL x_a  = quad->points[ia];
        const BS_REAL nu_a = (a < n) ? t * x_a : t / x_a;
        const BS_REAL Jac_a = (a < n) ? one : one / (x_a * x_a);
        const BS_REAL w_a  = quad->w[ia];

        BS_REAL g_a[total_num_species], bf_a[total_num_species];
        for (int s = 0; s < total_num_species; ++s)
        {
            g_a[s]  = TotalNuF(nu_a, &gop->distr_pars, s);
            bf_a[s] = neglect_blocking ? one : (one - g_a[s]);
        }

        for (int b = a; b < N2; ++b)
        {
            const int    ib    = (b < n) ? b : b - n;
            const BS_REAL x_b  = quad->points[ib];
            const BS_REAL nu_b = (b < n) ? t * x_b : t / x_b;
            const BS_REAL Jac_b = (b < n) ? one : one / (x_b * x_b);
            const BS_REAL w_b  = quad->w[ib];

            const BS_REAL nu_fourth = (nu_a * nu_a) * (nu_b * nu_b);
            const BS_REAL wJ        = w_a * w_b * Jac_a * Jac_b;

            BS_REAL g_b[total_num_species], bf_b[total_num_species];
            for (int s = 0; s < total_num_species; ++s)
            {
                g_b[s]  = TotalNuF(nu_b, &gop->distr_pars, s);
                bf_b[s] = neglect_blocking ? one : (one - g_b[s]);
            }

            // Raw (unweighted) kernels for entry [a][b] and its transpose [b][a]
            BS_REAL em_ab[total_num_species] = {0}, ab_ab[total_num_species] = {0};
            BS_REAL em_ba[total_num_species] = {0}, ab_ba[total_num_species] = {0};

            if (use_pair)
            {
                gop->kernel_pars.pair_kernel_params.cos_theta = one;
                gop->kernel_pars.pair_kernel_params.filter    = 0;
                gop->kernel_pars.pair_kernel_params.lmax      = 0;
                gop->kernel_pars.pair_kernel_params.mu        = one;
                gop->kernel_pars.pair_kernel_params.mu_prime  = one;
                gop->kernel_pars.pair_kernel_params.omega       = nu_a;
                gop->kernel_pars.pair_kernel_params.omega_prime = nu_b;
                MyKernelOutput pair_1, pair_2;
                PairKernels(&gop->eos_pars, &gop->kernel_pars.pair_kernel_params,
                            &pair_1, &pair_2);
                for (int s = 0; s < total_num_species; ++s)
                {
                    em_ab[s] += pair_1.em[s];  ab_ab[s] += pair_1.abs[s];
                    if (a != b) { em_ba[s] += pair_2.em[s]; ab_ba[s] += pair_2.abs[s]; }
                }
            }
            if (use_brem)
            {
                gop->kernel_pars.brem_kernel_params.omega       = nu_a;
                gop->kernel_pars.brem_kernel_params.omega_prime = nu_b;
                MyKernelOutput brem_ker;
                if (brem_impl == BREM_BRT06)
                {
                    brem_ker = BremKernelsBRT06(&gop->kernel_pars.brem_kernel_params,
                                                &gop->eos_pars);
                }
                else
                {
                    gop->kernel_pars.brem_kernel_params.l = 0;
                    gop->kernel_pars.brem_kernel_params.use_NN_medium_corr =
                        gop->opacity_pars.use_NN_medium_corr;
                    brem_ker = BremKernelsLegCoeff(&gop->kernel_pars.brem_kernel_params,
                                                   &gop->eos_pars);
                }
                for (int s = 0; s < total_num_species; ++s)
                {
                    em_ab[s] += brem_ker.em[0];  ab_ab[s] += brem_ker.abs[0];
                    if (a != b) { em_ba[s] += brem_ker.em[0]; ab_ba[s] += brem_ker.abs[0]; }
                }
            }

            // Weight (block factors + nu^4) and accumulate into the integrals.
            // n-contribution = wJ * value ; e-contribution = wJ * value * nu_row
            for (int s = 0; s < total_num_species; ++s)
            {
                const int pp = pair_partner[s];
                BS_REAL em_w, ab_w;
                // entry [a][b]: row = a (energy nu_a), col = b (energy nu_b)
                if (stim_abs == 1)
                {
                    em_w = nu_fourth * em_ab[s] * bf_b[pp];
                    ab_w = nu_fourth * g_a[s] * (em_ab[s] * bf_b[pp] + ab_ab[s] * g_b[pp]);
                }
                else
                {
                    em_w = nu_fourth * (one - g_a[s]) * em_ab[s] * bf_b[pp];
                    ab_w = nu_fourth * g_a[s] * ab_ab[s] * g_b[pp];
                }
                acc_n_em[s] += wJ * em_w;          acc_n_ab[s] += wJ * ab_w;
                acc_e_em[s] += wJ * em_w * nu_a;    acc_e_ab[s] += wJ * ab_w * nu_a;

                if (a != b)
                {
                    // entry [b][a]: row = b (energy nu_b), col = a (energy nu_a)
                    if (stim_abs == 1)
                    {
                        em_w = nu_fourth * em_ba[s] * bf_a[pp];
                        ab_w = nu_fourth * g_b[s] * (em_ba[s] * bf_a[pp] + ab_ba[s] * g_a[pp]);
                    }
                    else
                    {
                        em_w = nu_fourth * (one - g_b[s]) * em_ba[s] * bf_a[pp];
                        ab_w = nu_fourth * g_b[s] * ab_ba[s] * g_a[pp];
                    }
                    acc_n_em[s] += wJ * em_w;          acc_n_ab[s] += wJ * ab_w;
                    acc_e_em[s] += wJ * em_w * nu_b;    acc_e_ab[s] += wJ * ab_w * nu_b;
                }
            }
        }
    }

    const BS_REAL t_sqr = t * t;
    for (int s = 0; s < total_num_species; ++s)
    {
        n2d->integrand[0 + s]                 += acc_n_em[s] * t_sqr;
        n2d->integrand[total_num_species + s] += acc_n_ab[s] * t_sqr;
        e2d->integrand[0 + s]                 += acc_e_em[s] * t_sqr;
        e2d->integrand[total_num_species + s] += acc_e_ab[s] * t_sqr;
    }
}

// Helper: weighted NEPS kernel values for one (nu, nu_bar) energy pair.
// Fills the "[.][j]" cell (em_j/ab_j) and the "[.][n+j]" cell (em_nj/ab_nj),
// matching the inner body of ComputeNEPSIntegrand.
KOKKOS_INLINE_FUNCTION
void NEPSCellFused(GreyOpacityParams* gop, BS_REAL nu, BS_REAL nu_bar,
                   const int stim_abs, const bool neglect_blocking,
                   BS_REAL em_j[total_num_species],  BS_REAL ab_j[total_num_species],
                   BS_REAL em_nj[total_num_species], BS_REAL ab_nj[total_num_species])
{
    constexpr BS_REAL one = 1;
    const BS_REAL nu_fourth = (nu * nu) * (nu_bar * nu_bar);

    BS_REAL g_nu[total_num_species], g_nb[total_num_species];
    for (int s = 0; s < total_num_species; ++s)
    {
        g_nu[s] = TotalNuF(nu,     &gop->distr_pars, s);
        g_nb[s] = TotalNuF(nu_bar, &gop->distr_pars, s);
    }

    gop->kernel_pars.inelastic_kernel_params.omega       = nu;
    gop->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
    MyKernelOutput inel_1 =
        InelasticScattKernels(&gop->kernel_pars.inelastic_kernel_params, &gop->eos_pars);
    gop->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
    gop->kernel_pars.inelastic_kernel_params.omega_prime = nu;
    MyKernelOutput inel_2 =
        InelasticScattKernels(&gop->kernel_pars.inelastic_kernel_params, &gop->eos_pars);

    for (int s = 0; s < total_num_species; ++s)
    {
        const BS_REAL bf_nu = neglect_blocking ? one : (one - g_nu[s]);
        const BS_REAL bf_nb = neglect_blocking ? one : (one - g_nb[s]);
        const BS_REAL te1 = inel_1.em[s]  * g_nb[s];
        const BS_REAL ta1 = inel_1.abs[s] * bf_nb;
        const BS_REAL te2 = inel_2.em[s]  * g_nu[s];
        const BS_REAL ta2 = inel_2.abs[s] * bf_nu;
        if (stim_abs == 1)
        {
            ab_j[s]  = nu_fourth * g_nu[s] * (te1 + ta1);
            em_j[s]  = nu_fourth * te1;
            ab_nj[s] = nu_fourth * g_nb[s] * (te2 + ta2);
            em_nj[s] = nu_fourth * te2;
        }
        else
        {
            ab_j[s]  = nu_fourth * g_nu[s] * ta1;
            em_j[s]  = nu_fourth * (one - g_nu[s]) * te1;
            ab_nj[s] = nu_fourth * g_nb[s] * ta2;
            em_nj[s] = nu_fourth * (one - g_nb[s]) * te2;
        }
    }
}

// Fused inelastic-scattering (NEPS) double integral.  Equivalent to:
//   M1MatrixKokkos2D m = ComputeNEPSIntegrand(quad, t, gop, stim_abs);
//   GaussLegendreIntegrate2DMatrixForNEPS(quad, &m, t, n2d, e2d);
KOKKOS_INLINE_FUNCTION
void AddNEPSDoubleIntegralFused(const MyQuadrature* quad, BS_REAL t,
                                GreyOpacityParams* gop, const int stim_abs,
                                MyQuadratureIntegrand* n2d,
                                MyQuadratureIntegrand* e2d)
{
    constexpr BS_REAL half = 0.5, one = 1;
    const int n = quad->nx;
    const bool neglect_blocking = gop->opacity_pars.neglect_blocking;

    BS_REAL acc_n_em[total_num_species] = {0}, acc_n_ab[total_num_species] = {0};
    BS_REAL acc_e_em[total_num_species] = {0}, acc_e_ab[total_num_species] = {0};

    for (int i = 0; i < n; ++i)
    {
        const BS_REAL x_i  = quad->points[i];
        const BS_REAL w_i  = quad->w[i];
        const BS_REAL x3_i = x_i * x_i * x_i;

        for (int j = 0; j < n; ++j)
        {
            const BS_REAL x_j  = quad->points[j];
            const BS_REAL w_ij = w_i * quad->w[j];

            // variant 1 -> matrix cells [i][j] and [i][n+j]
            const BS_REAL nu1 = half * t * x_i * (one - x_j);
            const BS_REAL nb1 = half * t * x_i * (one + x_j);
            // variant 2 -> matrix cells [n+i][j] and [n+i][n+j]
            const BS_REAL nu2 = half * t * (one - x_j) / x_i;
            const BS_REAL nb2 = half * t * (one + x_j) / x_i;

            BS_REAL em_ij[total_num_species],  ab_ij[total_num_species];
            BS_REAL em_inj[total_num_species], ab_inj[total_num_species];
            BS_REAL em_nij[total_num_species], ab_nij[total_num_species];
            BS_REAL em_ninj[total_num_species],ab_ninj[total_num_species];
            NEPSCellFused(gop, nu1, nb1, stim_abs, neglect_blocking,
                          em_ij, ab_ij, em_inj, ab_inj);
            NEPSCellFused(gop, nu2, nb2, stim_abs, neglect_blocking,
                          em_nij, ab_nij, em_ninj, ab_ninj);

            for (int s = 0; s < total_num_species; ++s)
            {
                acc_n_em[s] += w_ij * (x_i * (em_ij[s] + em_inj[s]) +
                                       (em_nij[s] + em_ninj[s]) / x3_i);
                acc_n_ab[s] += w_ij * (x_i * (ab_ij[s] + ab_inj[s]) +
                                       (ab_nij[s] + ab_ninj[s]) / x3_i);
                acc_e_em[s] += w_ij * (x_i * (nu1 * em_ij[s] + nb1 * em_inj[s]) +
                                       (nu2 * em_nij[s] + nb2 * em_ninj[s]) / x3_i);
                acc_e_ab[s] += w_ij * (x_i * (nu1 * ab_ij[s] + nb1 * ab_inj[s]) +
                                       (nu2 * ab_nij[s] + nb2 * ab_ninj[s]) / x3_i);
            }
        }
    }

    const BS_REAL half_t_sqr = half * t * t;
    for (int s = 0; s < total_num_species; ++s)
    {
        n2d->integrand[0 + s]                 += acc_n_em[s] * half_t_sqr;
        n2d->integrand[total_num_species + s] += acc_n_ab[s] * half_t_sqr;
        e2d->integrand[0 + s]                 += acc_e_em[s] * half_t_sqr;
        e2d->integrand[total_num_species + s] += acc_e_ab[s] * half_t_sqr;
    }
}


/* Computes the opacities for the M1 code, with thermal and
 * non-thermal processes treated together.
 *
 * NEPS is treated together with other reactions.
 * NEPS is considered for computation of number emissivity (eta_0) 
 * and absorsivity (kappa_0_a).
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
        // and [MeV nm^-3]
        n[idx] = my_grey_opacity_params->m1_pars.n[idx];
        J[idx] = my_grey_opacity_params->m1_pars.J[idx];
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
        // Fused fill+reduce: no 25.6 KB M1MatrixKokkos2D on the stack.
        AddPairBremDoubleIntegralFused(quad_2d, s_pair, my_grey_opacity_params,
                                       stim_abs, &n_integrals_2d, &e_integrals_2d);
    }

    MyQuadratureIntegrand n_neps_2d = {0};
    MyQuadratureIntegrand e_neps_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_scatt == 1)
    {
        // Fused fill+reduce: no 25.6 KB M1MatrixKokkos2D on the stack.
        AddNEPSDoubleIntegralFused(quad_2d, four * s_neps, my_grey_opacity_params,
                                   stim_abs, &n_neps_2d, &e_neps_2d);
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


/* Computes the opacities for the M1 code, with thermal and
 * non-thermal processes treated separately.
 *
 * NEPS is treated separately from other reactions.
 * NEPS is NOT considered for computation of number emissivity (eta_0) 
 * and absorsivity (kappa_0_a).
 */
KOKKOS_INLINE_FUNCTION
M1OpacitiesNonThermalSeparated ComputeM1OpacitiesGenericFormalismNonThermalSeparated(
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
        // and [MeV nm^-3]
        n[idx] = my_grey_opacity_params->m1_pars.n[idx];
        J[idx] = my_grey_opacity_params->m1_pars.J[idx];
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
        // Fused fill+reduce: no 25.6 KB M1MatrixKokkos2D on the stack.
        AddPairBremDoubleIntegralFused(quad_2d, s_pair, my_grey_opacity_params,
                                       stim_abs, &n_integrals_2d, &e_integrals_2d);
    }

    MyQuadratureIntegrand n_neps_2d = {0};
    MyQuadratureIntegrand e_neps_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_scatt == 1)
    {
        // Fused fill+reduce: no 25.6 KB M1MatrixKokkos2D on the stack.
        AddNEPSDoubleIntegralFused(quad_2d, four * s_neps, my_grey_opacity_params,
                                   stim_abs, &n_neps_2d, &e_neps_2d);
    }

    M1OpacitiesNonThermalSeparated m1_opacities_non_th_separated = {0};

    /* Set all opacities to zero. They'll be left as 0 if the neutrino
    number/energy density is too low (to avoid inf/nan, since the number/energy
    density appears in the denominator). Note that emissivities do no need this
    precaution (and it also make sense physically: you can produce neutrinos
    even if there are none to start with). */
    constexpr BS_REAL zero = 0;
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        m1_opacities_non_th_separated.kappa_0_a[idx]      = zero;
        m1_opacities_non_th_separated.kappa_a_th[idx]     = zero;
        m1_opacities_non_th_separated.kappa_a_non_th[idx] = zero;
        m1_opacities_non_th_separated.kappa_s[idx]        = zero;
    }

    /* Electron neutrinos */
    m1_opacities_non_th_separated.eta_0[id_nue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[0] +
                          beta_n_em_integrals.integrand[id_nue]);

    m1_opacities_non_th_separated.eta_th[id_nue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[0] +
                          beta_j_em_integrals.integrand[id_nue]);

    m1_opacities_non_th_separated.eta_non_th[id_nue] =
        kBS_FourPi_hc3 * kBS_FourPi_hc3 * e_neps_2d.integrand[0];

    if (n[id_nue] > THRESHOLD_N)
    {
        m1_opacities_non_th_separated.kappa_0_a[id_nue] =
            kBS_FourPi_hc3 / (c_light * n[id_nue]) *
            (kBS_FourPi_hc3 * n_integrals_2d.integrand[4] +
             beta_n_abs_integrals.integrand[id_nue]);
    }
    if (J[id_nue] > THRESHOLD_J)
    {
        m1_opacities_non_th_separated.kappa_a_th[id_nue] =
            n[id_nue] == zero ?
                zero :
                kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                    (kBS_FourPi_hc3 * e_integrals_2d.integrand[4] +
                     beta_j_abs_integrals.integrand[id_nue]);

        m1_opacities_non_th_separated.kappa_a_non_th[id_nue] =
            n[id_nue] == zero ?
                zero :
                kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                    (kBS_FourPi_hc3 * e_neps_2d.integrand[4]);

        m1_opacities_non_th_separated.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                       iso_integrals.integrand[id_nue];
    }

    /* Electron anti-neutrinos */
    m1_opacities_non_th_separated.eta_0[id_anue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[1] +
                          beta_n_em_integrals.integrand[id_anue]);

    m1_opacities_non_th_separated.eta_th[id_anue] =
        kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[1] +
                          beta_j_em_integrals.integrand[id_anue]);

    m1_opacities_non_th_separated.eta_non_th[id_anue] =
        kBS_FourPi_hc3 * kBS_FourPi_hc3 * e_neps_2d.integrand[1];

    if (n[id_anue] > THRESHOLD_N)
    {
        m1_opacities_non_th_separated.kappa_0_a[id_anue] =
            kBS_FourPi_hc3 / (c_light * n[id_anue]) *
            (kBS_FourPi_hc3 * n_integrals_2d.integrand[5] +
             beta_n_abs_integrals.integrand[id_anue]);
    }
    if (J[id_anue] > THRESHOLD_J)
    {
        m1_opacities_non_th_separated.kappa_a_th[id_anue] =
            kBS_FourPi_hc3 / (c_light * J[id_anue]) *
            (kBS_FourPi_hc3 * e_integrals_2d.integrand[5] +
             beta_j_abs_integrals.integrand[id_anue]);

        m1_opacities_non_th_separated.kappa_a_non_th[id_anue] =
            kBS_FourPi_hc3 / (c_light * J[id_anue]) *
            (kBS_FourPi_hc3 * e_neps_2d.integrand[5]);

        m1_opacities_non_th_separated.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                        (c_light * J[id_anue]) *
                                        iso_integrals.integrand[id_anue];
    }

    /* Heavy neutrinos */
    m1_opacities_non_th_separated.eta_0[id_nux] =
        kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[2];

    m1_opacities_non_th_separated.eta_th[id_nux] =
        kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[2];

    m1_opacities_non_th_separated.eta_non_th[id_nux] =
        kBS_FourPi_hc3_sqr * e_neps_2d.integrand[2];

    if (n[id_nux] > THRESHOLD_N)
    {
        m1_opacities_non_th_separated.kappa_0_a[id_nux] =
            kBS_FourPi_hc3_sqr / (c_light * n[id_nux]) *
            n_integrals_2d.integrand[6];
    }
    if (J[id_nux] > THRESHOLD_J)
    {
        m1_opacities_non_th_separated.kappa_a_th[id_nux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
            e_integrals_2d.integrand[6];

        m1_opacities_non_th_separated.kappa_a_non_th[id_nux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
            e_neps_2d.integrand[6];

        m1_opacities_non_th_separated.kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J[id_nux]) *
                                       iso_integrals.integrand[id_nux];
    }

    /* Heavy anti-neutrinos */
    m1_opacities_non_th_separated.eta_0[id_anux] =
        kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[3];

    m1_opacities_non_th_separated.eta_th[id_anux] =
        kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[3];

    m1_opacities_non_th_separated.eta_non_th[id_anux] =
        kBS_FourPi_hc3_sqr * e_neps_2d.integrand[3];

    if (n[id_anux] > THRESHOLD_N)
    {
        m1_opacities_non_th_separated.kappa_0_a[id_anux] =
            n[id_anux] == zero ?
                zero :
                kBS_FourPi_hc3_sqr / (c_light * n[id_anux]) *
                    n_integrals_2d.integrand[7];
    }
    if (J[id_anux] > THRESHOLD_J)
    {
        m1_opacities_non_th_separated.kappa_a_th[id_anux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
            e_integrals_2d.integrand[7];

        m1_opacities_non_th_separated.kappa_a_non_th[id_anux] =
            kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
            e_neps_2d.integrand[7];

        m1_opacities_non_th_separated.kappa_s[id_anux] = kBS_FourPi_hc3 /
                                        (c_light * J[id_anux]) *
                                        iso_integrals.integrand[id_anux];
    }

    return m1_opacities_non_th_separated;
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


KOKKOS_INLINE_FUNCTION
M1OpacitiesNonThermalSeparated
ComputeM1OpacitiesNotStimulatedNonThermalSeparated(
                    MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                    GreyOpacityParams* my_grey_opacity_params)
{
    return ComputeM1OpacitiesGenericFormalismNonThermalSeparated(
                        quad_1d, quad_2d, my_grey_opacity_params, 0);
}


KOKKOS_INLINE_FUNCTION
M1OpacitiesNonThermalSeparated 
ComputeM1OpacitiesNonThermalSeparated(
            const MyQuadrature* quad_1d, const MyQuadrature* quad_2d, 
            GreyOpacityParams* my_grey_opacity_params)
{
    return ComputeM1OpacitiesGenericFormalismNonThermalSeparated(
                                quad_1d, quad_2d, my_grey_opacity_params, 1);
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
        if (opacity_pars.brem_implementation == BREM_BRT06)
        {
            brem_kernels_m1 = BremKernelsBRT06(
                &my_grey_opacity_params->kernel_pars.brem_kernel_params,
                &my_eos_params);
        }
        else if (opacity_pars.brem_implementation == BREM_HR98)
        {
            my_grey_opacity_params->kernel_pars.brem_kernel_params.l = 0;
            my_grey_opacity_params->kernel_pars.brem_kernel_params
                .use_NN_medium_corr =
                my_grey_opacity_params->opacity_pars.use_NN_medium_corr;
            brem_kernels_m1 = BremKernelsLegCoeff(
                &my_grey_opacity_params->kernel_pars.brem_kernel_params,
                &my_eos_params);
        }
        else if (opacity_pars.brem_implementation == BREM_GP19)
        {
            brem_kernels_m1 = BremKernelAbsGP19(
                &my_grey_opacity_params->kernel_pars.brem_kernel_params,
                &my_eos_params);
        }
        else
        {
            BS_ASSERT(false, "Unknown bremsstrahlung implementation: %d",
                      (int)opacity_pars.brem_implementation);
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

    // const BS_REAL eta_e = my_grey_opacity_params->eos_pars.mu_e /
    //                       my_grey_opacity_params->eos_pars.temp;

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

    sp_opacities.j[id_nue]  = abs_em_beta.em[id_nue] +
                              kBS_FourPi_hc3 * (integrals_pair_1d.integrand[0] +
                                                integrals_neps_1d.integrand[0]);
    sp_opacities.j[id_anue] = abs_em_beta.em[id_anue] +
                              kBS_FourPi_hc3 * (integrals_pair_1d.integrand[1] +
                                                integrals_neps_1d.integrand[1]);
    sp_opacities.j[id_nux]  = abs_em_beta.em[id_nux] +
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


#ifdef KOKKOS_CORE_HPP
// ===========================================================================
// Team-parallel opacity computation (GPU-optimized)
//
// One Kokkos team per grid cell. The M1MatrixKokkos2D scratch matrix lives in
// team shared memory instead of per-thread local memory, eliminating register
// spill. The 2D quadrature fill loops (pair+brem: 210 iters; NEPS: 100 iters)
// are distributed over TeamThreadRange. Each thread makes a register-local copy
// of MyKernelParams to avoid the race on omega/omega_prime. Pair+brem weighting
// (WeightNuNuBar + CommonWeights) is fused into the fill loop. Serial work
// (1D integrals, integration sums, assembly) runs via Kokkos::single(PerTeam).
// ===========================================================================

// Decode flat upper-triangular index (row<=col) for an N-wide matrix.
KOKKOS_INLINE_FUNCTION
void DecodeUpperTriIdx(int flat, int N, int &row, int &col) {
    row = 0;
    while (flat >= N - row) { flat -= N - row; ++row; }
    col = row + flat;
}

// ---------------------------------------------------------------------------
// ComputeDoubleIntegrandFillTeam
// Team-parallel version of ComputeDoubleIntegrand().
// Writes the fused (pair+brem+weight+nu^4) integrand into *mat (scratch).
// Emits a team_barrier before returning so *mat is visible to all threads.
// ---------------------------------------------------------------------------
template <class MemberType>
KOKKOS_INLINE_FUNCTION
void ComputeDoubleIntegrandFillTeam(const MemberType &team,
                                    const MyQuadrature *quad, BS_REAL t,
                                    GreyOpacityParams *gop,
                                    M1MatrixKokkos2D *mat,
                                    const int stim_abs) {
    const int n    = quad->nx;
    const int N2   = 2 * n;
    const int N_tri = N2 * (N2 + 1) / 2; // 210 for n=10

    // Zero scratch matrix (3200 doubles = 25600 bytes)
    const int mat_n = static_cast<int>(sizeof(M1MatrixKokkos2D) / sizeof(BS_REAL));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, mat_n), [&](const int idx) {
        reinterpret_cast<BS_REAL *>(mat)[idx] = BS_REAL(0);
    });
    team.team_barrier();

    const bool use_pair         = (gop->opacity_flags.use_pair == 1);
    const bool use_brem         = (gop->opacity_flags.use_brem == 1);
    const bool neglect_blocking = gop->opacity_pars.neglect_blocking;
    const BremImpl brem_impl    = gop->opacity_pars.brem_implementation;
    // pair_partner[s]: partner species for blocking (nue<->anue, nux<->anux)
    constexpr int pair_partner[4] = {id_anue, id_nue, id_anux, id_nux};

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_tri), [&](const int flat) {
        int row, col;
        DecodeUpperTriIdx(flat, N2, row, col);

        const BS_REAL nu_row = (row < n) ? t * quad->points[row]
                                         : t / quad->points[row - n];
        const BS_REAL nu_col = (col < n) ? t * quad->points[col]
                                         : t / quad->points[col - n];
        const BS_REAL nu_fourth = (nu_row * nu_row) * (nu_col * nu_col);

        // Thread-local copy avoids race on gop->kernel_pars.{pair,brem}_kernel_params.omega
        MyKernelParams kp = gop->kernel_pars;

        BS_REAL em_rc[total_num_species] = {BS_REAL(0)};
        BS_REAL ab_rc[total_num_species] = {BS_REAL(0)};
        BS_REAL em_cr[total_num_species] = {BS_REAL(0)};
        BS_REAL ab_cr[total_num_species] = {BS_REAL(0)};

        if (use_pair) {
            kp.pair_kernel_params.cos_theta  = BS_REAL(1);
            kp.pair_kernel_params.filter     = BS_REAL(0);
            kp.pair_kernel_params.lmax       = BS_REAL(0);
            kp.pair_kernel_params.mu         = BS_REAL(1);
            kp.pair_kernel_params.mu_prime   = BS_REAL(1);
            kp.pair_kernel_params.omega       = nu_row;
            kp.pair_kernel_params.omega_prime = nu_col;
            MyKernelOutput pair_1, pair_2;
            PairKernels(&gop->eos_pars, &kp.pair_kernel_params, &pair_1, &pair_2);
            for (int s = 0; s < total_num_species; ++s) {
                em_rc[s] += pair_1.em[s];
                ab_rc[s] += pair_1.abs[s];
                if (row != col) { em_cr[s] += pair_2.em[s]; ab_cr[s] += pair_2.abs[s]; }
            }
        }

        if (use_brem) {
            kp.brem_kernel_params.omega       = nu_row;
            kp.brem_kernel_params.omega_prime = nu_col;
            MyKernelOutput brem_ker;
            if (brem_impl == BREM_BRT06) {
                brem_ker = BremKernelsBRT06(&kp.brem_kernel_params, &gop->eos_pars);
            } else {
                kp.brem_kernel_params.l = 0;
                kp.brem_kernel_params.use_NN_medium_corr = gop->opacity_pars.use_NN_medium_corr;
                brem_ker = BremKernelsLegCoeff(&kp.brem_kernel_params, &gop->eos_pars);
            }
            for (int s = 0; s < total_num_species; ++s) {
                em_rc[s] += brem_ker.em[0]; ab_rc[s] += brem_ker.abs[0];
                if (row != col) { em_cr[s] += brem_ker.em[0]; ab_cr[s] += brem_ker.abs[0]; }
            }
        }

        BS_REAL g_row[total_num_species], g_col[total_num_species];
        BS_REAL bf_row[total_num_species], bf_col[total_num_species];
        for (int s = 0; s < total_num_species; ++s) {
            g_row[s]  = TotalNuF(nu_row, &gop->distr_pars, s);
            g_col[s]  = TotalNuF(nu_col, &gop->distr_pars, s);
            bf_row[s] = neglect_blocking ? BS_REAL(1) : BS_REAL(1) - g_row[s];
            bf_col[s] = neglect_blocking ? BS_REAL(1) : BS_REAL(1) - g_col[s];
        }

        if (row == col) {
            for (int s = 0; s < total_num_species; ++s) {
                const int pp = pair_partner[s];
                const BS_REAL em_w = em_rc[s] * bf_row[pp];
                const BS_REAL ab_w = ab_rc[s] * g_row[pp];
                if (stim_abs == 1) {
                    mat->m1_mat_em[s][row][row] = nu_fourth * em_w;
                    mat->m1_mat_ab[s][row][row] = nu_fourth * g_row[s] * (em_w + ab_w);
                } else {
                    mat->m1_mat_em[s][row][row] = nu_fourth * (BS_REAL(1) - g_row[s]) * em_w;
                    mat->m1_mat_ab[s][row][row] = nu_fourth * g_row[s] * ab_w;
                }
            }
        } else {
            for (int s = 0; s < total_num_species; ++s) {
                const int pp = pair_partner[s];
                const BS_REAL em_rc_w = em_rc[s] * bf_col[pp];
                const BS_REAL ab_rc_w = ab_rc[s] * g_col[pp];
                const BS_REAL em_cr_w = em_cr[s] * bf_row[pp];
                const BS_REAL ab_cr_w = ab_cr[s] * g_row[pp];
                if (stim_abs == 1) {
                    mat->m1_mat_em[s][row][col] = nu_fourth * em_rc_w;
                    mat->m1_mat_ab[s][row][col] = nu_fourth * g_row[s] * (em_rc_w + ab_rc_w);
                    mat->m1_mat_em[s][col][row] = nu_fourth * em_cr_w;
                    mat->m1_mat_ab[s][col][row] = nu_fourth * g_col[s] * (em_cr_w + ab_cr_w);
                } else {
                    mat->m1_mat_em[s][row][col] = nu_fourth * (BS_REAL(1)-g_row[s]) * em_rc_w;
                    mat->m1_mat_ab[s][row][col] = nu_fourth * g_row[s] * ab_rc_w;
                    mat->m1_mat_em[s][col][row] = nu_fourth * (BS_REAL(1)-g_col[s]) * em_cr_w;
                    mat->m1_mat_ab[s][col][row] = nu_fourth * g_col[s] * ab_cr_w;
                }
            }
        }
    });
    team.team_barrier();
}

// ---------------------------------------------------------------------------
// ComputeNEPSIntegrandFillTeam
// Team-parallel version of ComputeNEPSIntegrand().
// For each (i,j) in [0,n)x[0,n) two (nu,nu_bar) variants are computed,
// each writing to two unique matrix entries — all 4n^2 writes are race-free.
// ---------------------------------------------------------------------------
template <class MemberType>
KOKKOS_INLINE_FUNCTION
void ComputeNEPSIntegrandFillTeam(const MemberType &team,
                                   const MyQuadrature *quad, BS_REAL t,
                                   GreyOpacityParams *gop,
                                   M1MatrixKokkos2D *mat,
                                   const int stim_abs) {
    const int n    = quad->nx;
    const int N_sq = n * n;
    constexpr BS_REAL half = BS_REAL(0.5);
    constexpr BS_REAL one  = BS_REAL(1);

    const int mat_n = static_cast<int>(sizeof(M1MatrixKokkos2D) / sizeof(BS_REAL));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, mat_n), [&](const int idx) {
        reinterpret_cast<BS_REAL *>(mat)[idx] = BS_REAL(0);
    });
    team.team_barrier();

    const bool neglect_blocking = gop->opacity_pars.neglect_blocking;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N_sq), [&](const int flat) {
        const int i    = flat / n;
        const int j    = flat % n;
        const BS_REAL x_i = quad->points[i];
        const BS_REAL x_j = quad->points[j];
        MyKernelParams kp = gop->kernel_pars;

        // Variant 1: nu ~ x_i  → writes [i][j] and [i][n+j]
        {
            const BS_REAL nu     = half * t * x_i * (one - x_j);
            const BS_REAL nu_bar = half * t * x_i * (one + x_j);
            const BS_REAL nf     = (nu * nu) * (nu_bar * nu_bar);
            BS_REAL g_nu[total_num_species], g_nb[total_num_species];
            BS_REAL bf_nu[total_num_species], bf_nb[total_num_species];
            for (int s = 0; s < total_num_species; ++s) {
                g_nu[s]  = TotalNuF(nu,     &gop->distr_pars, s);
                g_nb[s]  = TotalNuF(nu_bar, &gop->distr_pars, s);
                bf_nu[s] = neglect_blocking ? one : one - g_nu[s];
                bf_nb[s] = neglect_blocking ? one : one - g_nb[s];
            }
            kp.inelastic_kernel_params.omega       = nu;
            kp.inelastic_kernel_params.omega_prime = nu_bar;
            MyKernelOutput inel_1 = InelasticScattKernels(&kp.inelastic_kernel_params, &gop->eos_pars);
            kp.inelastic_kernel_params.omega       = nu_bar;
            kp.inelastic_kernel_params.omega_prime = nu;
            MyKernelOutput inel_2 = InelasticScattKernels(&kp.inelastic_kernel_params, &gop->eos_pars);
            for (int s = 0; s < total_num_species; ++s) {
                const BS_REAL te1 = inel_1.em[s]  * g_nb[s];
                const BS_REAL ta1 = inel_1.abs[s] * bf_nb[s];
                const BS_REAL te2 = inel_2.em[s]  * g_nu[s];
                const BS_REAL ta2 = inel_2.abs[s] * bf_nu[s];
                if (stim_abs == 1) {
                    mat->m1_mat_em[s][i][j]     = nf * te1;
                    mat->m1_mat_ab[s][i][j]     = nf * g_nu[s] * (te1 + ta1);
                    mat->m1_mat_em[s][i][n + j] = nf * te2;
                    mat->m1_mat_ab[s][i][n + j] = nf * g_nb[s] * (te2 + ta2);
                } else {
                    mat->m1_mat_em[s][i][j]     = nf * (one - g_nu[s]) * te1;
                    mat->m1_mat_ab[s][i][j]     = nf * g_nu[s] * ta1;
                    mat->m1_mat_em[s][i][n + j] = nf * (one - g_nb[s]) * te2;
                    mat->m1_mat_ab[s][i][n + j] = nf * g_nb[s] * ta2;
                }
            }
        }

        // Variant 2: nu ~ 1/x_i  → writes [n+i][j] and [n+i][n+j]
        {
            const BS_REAL nu     = half * t * (one - x_j) / x_i;
            const BS_REAL nu_bar = half * t * (one + x_j) / x_i;
            const BS_REAL nf     = (nu * nu) * (nu_bar * nu_bar);
            BS_REAL g_nu[total_num_species], g_nb[total_num_species];
            BS_REAL bf_nu[total_num_species], bf_nb[total_num_species];
            for (int s = 0; s < total_num_species; ++s) {
                g_nu[s]  = TotalNuF(nu,     &gop->distr_pars, s);
                g_nb[s]  = TotalNuF(nu_bar, &gop->distr_pars, s);
                bf_nu[s] = neglect_blocking ? one : one - g_nu[s];
                bf_nb[s] = neglect_blocking ? one : one - g_nb[s];
            }
            kp.inelastic_kernel_params.omega       = nu;
            kp.inelastic_kernel_params.omega_prime = nu_bar;
            MyKernelOutput inel_1 = InelasticScattKernels(&kp.inelastic_kernel_params, &gop->eos_pars);
            kp.inelastic_kernel_params.omega       = nu_bar;
            kp.inelastic_kernel_params.omega_prime = nu;
            MyKernelOutput inel_2 = InelasticScattKernels(&kp.inelastic_kernel_params, &gop->eos_pars);
            for (int s = 0; s < total_num_species; ++s) {
                const BS_REAL te1 = inel_1.em[s]  * g_nb[s];
                const BS_REAL ta1 = inel_1.abs[s] * bf_nb[s];
                const BS_REAL te2 = inel_2.em[s]  * g_nu[s];
                const BS_REAL ta2 = inel_2.abs[s] * bf_nu[s];
                if (stim_abs == 1) {
                    mat->m1_mat_em[s][n + i][j]     = nf * te1;
                    mat->m1_mat_ab[s][n + i][j]     = nf * g_nu[s] * (te1 + ta1);
                    mat->m1_mat_em[s][n + i][n + j] = nf * te2;
                    mat->m1_mat_ab[s][n + i][n + j] = nf * g_nb[s] * (te2 + ta2);
                } else {
                    mat->m1_mat_em[s][n + i][j]     = nf * (one - g_nu[s]) * te1;
                    mat->m1_mat_ab[s][n + i][j]     = nf * g_nu[s] * ta1;
                    mat->m1_mat_em[s][n + i][n + j] = nf * (one - g_nb[s]) * te2;
                    mat->m1_mat_ab[s][n + i][n + j] = nf * g_nb[s] * ta2;
                }
            }
        }
    });
    team.team_barrier();
}

// ---------------------------------------------------------------------------
// ComputeM1OpacitiesTeam
// Team-parallel equivalent of ComputeM1Opacities() (always stim_abs=1).
// gop must already be filled in scratch (thread 0 + barrier before calling).
// *result is written by team thread 0; caller should team_barrier() after.
// ---------------------------------------------------------------------------
template <class MemberType>
KOKKOS_INLINE_FUNCTION
void ComputeM1OpacitiesTeam(const MemberType &team,
                             const MyQuadrature *quad_1d,
                             const MyQuadrature *quad_2d,
                             GreyOpacityParams *gop,
                             M1MatrixKokkos2D *mat,
                             M1Opacities *result) {
    constexpr int stim_abs = 1;
    constexpr BS_REAL c_light       = kBS_Clight;
    constexpr BS_REAL three_halves  = BS_REAL(1.5);
    constexpr BS_REAL five_sixths   = BS_REAL(5) / BS_REAL(6);
    constexpr BS_REAL five          = BS_REAL(5);
    constexpr BS_REAL temp_multiple = BS_REAL(0.5) * BS_REAL(4.364);

    const BS_REAL temp   = gop->eos_pars.temp;
    const BS_REAL eta_e  = gop->eos_pars.mu_e / temp;
    const BS_REAL s_pair = temp_multiple * temp;
    const BS_REAL s_nux  = three_halves * temp;
    const BS_REAL s_neps4 = BS_REAL(4) * temp_multiple * temp;

    // --- Phase 1: pair+brem fill (all threads) ---
    MyQuadratureIntegrand n2d = {0}, e2d = {0};
    if (gop->opacity_flags.use_pair || gop->opacity_flags.use_brem) {
        ComputeDoubleIntegrandFillTeam(team, quad_2d, s_pair, gop, mat, stim_abs);
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            GaussLegendreIntegrate2DMatrixForM1Coeffs(quad_2d, mat, s_pair, &n2d, &e2d);
        });
    }
    team.team_barrier();

    // --- Phase 2: NEPS fill (all threads, mat reused) ---
    MyQuadratureIntegrand nn = {0}, en = {0};
    if (gop->opacity_flags.use_inelastic_scatt) {
        ComputeNEPSIntegrandFillTeam(team, quad_2d, s_neps4, gop, mat, stim_abs);
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            GaussLegendreIntegrate2DMatrixForNEPS(quad_2d, mat, s_neps4, &nn, &en);
        });
    }
    team.team_barrier();

    // --- 1D integrals + assembly (thread 0) ---
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
        const BS_REAL *n_m1 = gop->m1_pars.n;
        const BS_REAL *J_m1 = gop->m1_pars.J;
        BS_REAL s_beta[total_num_species] = {0};
        BS_REAL s_iso[total_num_species]  = {0};

        if (eta_e > BS_REAL(-30.) && eta_e < BS_REAL(30.)) {
            s_beta[id_nue]  = temp * FDI_p5(eta_e)  / FDI_p4(eta_e);
            s_beta[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
        } else if (eta_e > BS_REAL(30.)) {
            s_beta[id_nue]  = temp * eta_e * five_sixths;
            s_beta[id_anue] = temp * five;
        } else {
            s_beta[id_nue]  = temp * five;
            s_beta[id_anue] = temp * eta_e * five_sixths;
        }
        for (int idx = 0; idx < total_num_species; ++idx)
            s_iso[idx] = (n_m1[idx] > THRESHOLD_N) ? (J_m1[idx] / n_m1[idx]) : s_nux;

        MyQuadratureIntegrand iso = {0};
        if (gop->opacity_flags.use_iso) {
            BS_REAL out_iso[total_num_species][BS_N_MAX];
            Scattering1DIntegrand(quad_1d, gop, s_iso, out_iso);
            iso = GaussLegendreIntegrate1DMatrix(quad_1d, total_num_species, out_iso, s_iso);
        }
        MyQuadratureIntegrand bne = {0}, bje = {0}, bna = {0}, bja = {0};
        if (gop->opacity_flags.use_abs_em) {
            BS_REAL obe[total_num_species][BS_N_MAX], oba[total_num_species][BS_N_MAX];
            Beta1DIntegrand(quad_1d, gop, s_beta, obe, oba, stim_abs);
            GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, obe, s_beta, &bne, &bje);
            GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, oba, s_beta, &bna, &bja);
        }

        *result = {0};
        constexpr BS_REAL z = BS_REAL(0);
        for (int s = 0; s < total_num_species; ++s) {
            result->kappa_0_a[s] = z; result->kappa_a[s] = z; result->kappa_s[s] = z;
        }

        // nue
        result->eta_0[id_nue] = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * (n2d.integrand[0] + nn.integrand[0]) + bne.integrand[id_nue]);
        result->eta[id_nue]   = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * (e2d.integrand[0] + en.integrand[0]) + bje.integrand[id_nue]);
        if (n_m1[id_nue] > THRESHOLD_N)
            result->kappa_0_a[id_nue] = kBS_FourPi_hc3 / (c_light * n_m1[id_nue]) *
                (kBS_FourPi_hc3 * (n2d.integrand[4] + nn.integrand[4]) + bna.integrand[id_nue]);
        if (J_m1[id_nue] > THRESHOLD_J) {
            result->kappa_a[id_nue] = (n_m1[id_nue] == z) ? z :
                kBS_FourPi_hc3 / (c_light * J_m1[id_nue]) *
                (kBS_FourPi_hc3 * (e2d.integrand[4] + en.integrand[4]) + bja.integrand[id_nue]);
            result->kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J_m1[id_nue]) * iso.integrand[id_nue];
        }

        // anue
        result->eta_0[id_anue] = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * (n2d.integrand[1] + nn.integrand[1]) + bne.integrand[id_anue]);
        result->eta[id_anue]   = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * (e2d.integrand[1] + en.integrand[1]) + bje.integrand[id_anue]);
        if (n_m1[id_anue] > THRESHOLD_N)
            result->kappa_0_a[id_anue] = kBS_FourPi_hc3 / (c_light * n_m1[id_anue]) *
                (kBS_FourPi_hc3 * (n2d.integrand[5] + nn.integrand[5]) + bna.integrand[id_anue]);
        if (J_m1[id_anue] > THRESHOLD_J) {
            result->kappa_a[id_anue] = kBS_FourPi_hc3 / (c_light * J_m1[id_anue]) *
                (kBS_FourPi_hc3 * (e2d.integrand[5] + en.integrand[5]) + bja.integrand[id_anue]);
            result->kappa_s[id_anue] = kBS_FourPi_hc3 / (c_light * J_m1[id_anue]) * iso.integrand[id_anue];
        }

        // nux
        result->eta_0[id_nux] = kBS_FourPi_hc3_sqr * (n2d.integrand[2] + nn.integrand[2]);
        result->eta[id_nux]   = kBS_FourPi_hc3_sqr * (e2d.integrand[2] + en.integrand[2]);
        if (n_m1[id_nux] > THRESHOLD_N)
            result->kappa_0_a[id_nux] = kBS_FourPi_hc3_sqr / (c_light * n_m1[id_nux]) *
                (n2d.integrand[6] + nn.integrand[6]);
        if (J_m1[id_nux] > THRESHOLD_J) {
            result->kappa_a[id_nux] = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_nux]) *
                (e2d.integrand[6] + en.integrand[6]);
            result->kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J_m1[id_nux]) * iso.integrand[id_nux];
        }

        // anux
        result->eta_0[id_anux] = kBS_FourPi_hc3_sqr * (n2d.integrand[3] + nn.integrand[3]);
        result->eta[id_anux]   = kBS_FourPi_hc3_sqr * (e2d.integrand[3] + en.integrand[3]);
        if (n_m1[id_anux] > THRESHOLD_N)
            result->kappa_0_a[id_anux] = (n_m1[id_anux] == z) ? z :
                kBS_FourPi_hc3_sqr / (c_light * n_m1[id_anux]) *
                (n2d.integrand[7] + nn.integrand[7]);
        if (J_m1[id_anux] > THRESHOLD_J) {
            result->kappa_a[id_anux] = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_anux]) *
                (e2d.integrand[7] + en.integrand[7]);
            result->kappa_s[id_anux] = kBS_FourPi_hc3 / (c_light * J_m1[id_anux]) * iso.integrand[id_anux];
        }
    });
}

// ---------------------------------------------------------------------------
// ComputeM1OpacitiesNonThermalSeparatedTeam
// Team-parallel version with NEPS contributions kept separate (non-thermal
// separation formalism). Mirrors ComputeM1OpacitiesGenericFormalismNonThermalSeparated.
// ---------------------------------------------------------------------------
template <class MemberType>
KOKKOS_INLINE_FUNCTION
void ComputeM1OpacitiesNonThermalSeparatedTeam(
        const MemberType &team,
        const MyQuadrature *quad_1d,
        const MyQuadrature *quad_2d,
        GreyOpacityParams *gop,
        M1MatrixKokkos2D *mat,
        M1OpacitiesNonThermalSeparated *result) {
    constexpr int stim_abs = 1;
    constexpr BS_REAL c_light       = kBS_Clight;
    constexpr BS_REAL three_halves  = BS_REAL(1.5);
    constexpr BS_REAL five_sixths   = BS_REAL(5) / BS_REAL(6);
    constexpr BS_REAL five          = BS_REAL(5);
    constexpr BS_REAL temp_multiple = BS_REAL(0.5) * BS_REAL(4.364);

    const BS_REAL temp   = gop->eos_pars.temp;
    const BS_REAL eta_e  = gop->eos_pars.mu_e / temp;
    const BS_REAL s_pair = temp_multiple * temp;
    const BS_REAL s_nux  = three_halves * temp;
    const BS_REAL s_neps4 = BS_REAL(4) * temp_multiple * temp;

    MyQuadratureIntegrand n2d = {0}, e2d = {0};
    if (gop->opacity_flags.use_pair || gop->opacity_flags.use_brem) {
        ComputeDoubleIntegrandFillTeam(team, quad_2d, s_pair, gop, mat, stim_abs);
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            GaussLegendreIntegrate2DMatrixForM1Coeffs(quad_2d, mat, s_pair, &n2d, &e2d);
        });
    }
    team.team_barrier();

    MyQuadratureIntegrand nn = {0}, en = {0};
    if (gop->opacity_flags.use_inelastic_scatt) {
        ComputeNEPSIntegrandFillTeam(team, quad_2d, s_neps4, gop, mat, stim_abs);
        Kokkos::single(Kokkos::PerTeam(team), [&]() {
            GaussLegendreIntegrate2DMatrixForNEPS(quad_2d, mat, s_neps4, &nn, &en);
        });
    }
    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
        const BS_REAL *n_m1 = gop->m1_pars.n;
        const BS_REAL *J_m1 = gop->m1_pars.J;
        BS_REAL s_beta[total_num_species] = {0};
        BS_REAL s_iso[total_num_species]  = {0};

        if (eta_e > BS_REAL(-30.) && eta_e < BS_REAL(30.)) {
            s_beta[id_nue]  = temp * FDI_p5(eta_e)  / FDI_p4(eta_e);
            s_beta[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
        } else if (eta_e > BS_REAL(30.)) {
            s_beta[id_nue]  = temp * eta_e * five_sixths;
            s_beta[id_anue] = temp * five;
        } else {
            s_beta[id_nue]  = temp * five;
            s_beta[id_anue] = temp * eta_e * five_sixths;
        }
        for (int idx = 0; idx < total_num_species; ++idx)
            s_iso[idx] = (n_m1[idx] > THRESHOLD_N) ? (J_m1[idx] / n_m1[idx]) : s_nux;

        MyQuadratureIntegrand iso = {0};
        if (gop->opacity_flags.use_iso) {
            BS_REAL out_iso[total_num_species][BS_N_MAX];
            Scattering1DIntegrand(quad_1d, gop, s_iso, out_iso);
            iso = GaussLegendreIntegrate1DMatrix(quad_1d, total_num_species, out_iso, s_iso);
        }
        MyQuadratureIntegrand bne = {0}, bje = {0}, bna = {0}, bja = {0};
        if (gop->opacity_flags.use_abs_em) {
            BS_REAL obe[total_num_species][BS_N_MAX], oba[total_num_species][BS_N_MAX];
            Beta1DIntegrand(quad_1d, gop, s_beta, obe, oba, stim_abs);
            GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, obe, s_beta, &bne, &bje);
            GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, oba, s_beta, &bna, &bja);
        }

        *result = {0};
        constexpr BS_REAL z = BS_REAL(0);
        for (int s = 0; s < total_num_species; ++s) {
            result->kappa_0_a[s] = z; result->kappa_a_th[s] = z;
            result->kappa_a_non_th[s] = z; result->kappa_s[s] = z;
        }

        // nue
        result->eta_0[id_nue]      = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * n2d.integrand[0] + bne.integrand[id_nue]);
        result->eta_th[id_nue]     = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * e2d.integrand[0] + bje.integrand[id_nue]);
        result->eta_non_th[id_nue] = kBS_FourPi_hc3_sqr * en.integrand[0];
        if (n_m1[id_nue] > THRESHOLD_N)
            result->kappa_0_a[id_nue] = kBS_FourPi_hc3 / (c_light * n_m1[id_nue]) *
                (kBS_FourPi_hc3 * n2d.integrand[4] + bna.integrand[id_nue]);
        if (J_m1[id_nue] > THRESHOLD_J) {
            result->kappa_a_th[id_nue]     = (n_m1[id_nue] == z) ? z :
                kBS_FourPi_hc3 / (c_light * J_m1[id_nue]) *
                (kBS_FourPi_hc3 * e2d.integrand[4] + bja.integrand[id_nue]);
            result->kappa_a_non_th[id_nue] = (n_m1[id_nue] == z) ? z :
                kBS_FourPi_hc3_sqr / (c_light * J_m1[id_nue]) * en.integrand[4];
            result->kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J_m1[id_nue]) * iso.integrand[id_nue];
        }

        // anue
        result->eta_0[id_anue]      = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * n2d.integrand[1] + bne.integrand[id_anue]);
        result->eta_th[id_anue]     = kBS_FourPi_hc3 *
            (kBS_FourPi_hc3 * e2d.integrand[1] + bje.integrand[id_anue]);
        result->eta_non_th[id_anue] = kBS_FourPi_hc3_sqr * en.integrand[1];
        if (n_m1[id_anue] > THRESHOLD_N)
            result->kappa_0_a[id_anue] = kBS_FourPi_hc3 / (c_light * n_m1[id_anue]) *
                (kBS_FourPi_hc3 * n2d.integrand[5] + bna.integrand[id_anue]);
        if (J_m1[id_anue] > THRESHOLD_J) {
            result->kappa_a_th[id_anue]     = kBS_FourPi_hc3 / (c_light * J_m1[id_anue]) *
                (kBS_FourPi_hc3 * e2d.integrand[5] + bja.integrand[id_anue]);
            result->kappa_a_non_th[id_anue] = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_anue]) * en.integrand[5];
            result->kappa_s[id_anue] = kBS_FourPi_hc3 / (c_light * J_m1[id_anue]) * iso.integrand[id_anue];
        }

        // nux
        result->eta_0[id_nux]      = kBS_FourPi_hc3_sqr * n2d.integrand[2];
        result->eta_th[id_nux]     = kBS_FourPi_hc3_sqr * e2d.integrand[2];
        result->eta_non_th[id_nux] = kBS_FourPi_hc3_sqr * en.integrand[2];
        if (n_m1[id_nux] > THRESHOLD_N)
            result->kappa_0_a[id_nux] = kBS_FourPi_hc3_sqr / (c_light * n_m1[id_nux]) * n2d.integrand[6];
        if (J_m1[id_nux] > THRESHOLD_J) {
            result->kappa_a_th[id_nux]     = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_nux]) * e2d.integrand[6];
            result->kappa_a_non_th[id_nux] = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_nux]) * en.integrand[6];
            result->kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J_m1[id_nux]) * iso.integrand[id_nux];
        }

        // anux
        result->eta_0[id_anux]      = kBS_FourPi_hc3_sqr * n2d.integrand[3];
        result->eta_th[id_anux]     = kBS_FourPi_hc3_sqr * e2d.integrand[3];
        result->eta_non_th[id_anux] = kBS_FourPi_hc3_sqr * en.integrand[3];
        if (n_m1[id_anux] > THRESHOLD_N)
            result->kappa_0_a[id_anux] = (n_m1[id_anux] == z) ? z :
                kBS_FourPi_hc3_sqr / (c_light * n_m1[id_anux]) * n2d.integrand[7];
        if (J_m1[id_anux] > THRESHOLD_J) {
            result->kappa_a_th[id_anux]     = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_anux]) * e2d.integrand[7];
            result->kappa_a_non_th[id_anux] = kBS_FourPi_hc3_sqr / (c_light * J_m1[id_anux]) * en.integrand[7];
            result->kappa_s[id_anux] = kBS_FourPi_hc3 / (c_light * J_m1[id_anux]) * iso.integrand[id_anux];
        }
    });
}

#endif // KOKKOS_CORE_HPP

#endif // BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_HPP_

