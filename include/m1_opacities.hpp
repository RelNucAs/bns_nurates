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
void ElBeta1DIntegrand(const MyQuadrature* quad, GreyOpacityParams* grey_pars,
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

            abs_em_beta = ElStimAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElStimAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElStimAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElStimAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_anue][n + i] = nu_sqr * abs_em_beta.em[id_anue];
            out_ab[id_anue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue];
        }
    }

    return;
}


KOKKOS_INLINE_FUNCTION
void MuonicBeta1DIntegrand(const MyQuadrature* quad, GreyOpacityParams* grey_pars,
                     const BS_REAL* t, const BS_REAL Q, BS_REAL out_em[][BS_N_MAX],
                     BS_REAL out_ab[][BS_N_MAX], const int stim_abs)
{
    const int n = quad->nx;
    BS_REAL nu, nu_sqr, g_nu;
    MyOpacity abs_em_beta;

    if (stim_abs == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            // nu_mu_1
            nu = (1 + quad->points[i]) * kBS_Mmu - Q;
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_num]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_num], i, quad->points[i]);
            nu_sqr = POW2(nu);

            g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_num);

            abs_em_beta = MuonStimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_num][i] = nu_sqr * abs_em_beta.em[id_num];
            // ab = em + ab (stimulated absorption)
            out_ab[id_num][i] = nu_sqr * g_nu * abs_em_beta.abs[id_num];

            // anu_mu_1
            nu = (1 + quad->points[i]) * kBS_Mmu + Q;
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anum]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anum], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anum);

            abs_em_beta = MuonStimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_anum][i] = nu_sqr * abs_em_beta.em[id_anum];
            // ab = em + ab (stimulated absorption)
            out_ab[id_anum][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anum];

            // nu_mu_2
            nu = t[id_num] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_num]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_num], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_num);

            abs_em_beta = MuonStimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_num][n + i] = nu_sqr * abs_em_beta.em[id_num];
            // ab = em + ab (stimulated absorption)
            out_ab[id_num][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_num];

            // anu_mu_2
            nu = t[id_anum] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anum]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anum], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anum);

            abs_em_beta = MuonStimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            out_em[id_anum][n + i] = nu_sqr * abs_em_beta.em[id_anum];

            // ab = em + ab (stimulated absorption)
            out_ab[id_anum][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anum];
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            // nu_mu_1
            nu = (1 + quad->points[i]) * kBS_Mmu - Q;
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_num]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_num], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_num);

            abs_em_beta = MuonAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_num][i] = nu_sqr * abs_em_beta.em[id_num];
            out_ab[id_num][i] = nu_sqr * g_nu * abs_em_beta.abs[id_num];

            // anu_mu_1
            nu = (1 + quad->points[i]) * kBS_Mmu + Q;
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anum]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anum], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anum);

            abs_em_beta = MuonAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_anum][i] = nu_sqr * abs_em_beta.em[id_anum];
            out_ab[id_anum][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anum];

            // nu_mu_2
            nu = t[id_num] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_num]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_num], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_num);

            abs_em_beta = MuonAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_num][n + i] = nu_sqr * abs_em_beta.em[id_num];
            out_ab[id_num][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_num];

            // anu_mu_2
            nu = t[id_anum] / quad->points[i];
            BS_ASSERT(nu >= 0,
                      "Neutrino energy is negative (nu=%e, t[id_anum]=%e, "
                      "quad->points[%d]=%e)",
                      nu, t[id_anum], i, quad->points[i]);
            nu_sqr = POW2(nu);
            g_nu   = TotalNuF(nu, &grey_pars->distr_pars, id_anum);

            abs_em_beta = MuonAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            out_em[id_anum][n + i] = nu_sqr * abs_em_beta.em[id_anum];
            out_ab[id_anum][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anum];
        }
    }

    return;
}


/*
KOKKOS_INLINE_FUNCTION
void AddElBetaReactionToIntegrand(int n, BS_REAL* nu_array,
                                GreyOpacityParams* grey_pars,
                                M1MatrixKokkos2D* out, const int stim_abs)
{
    BS_REAL nu;
    MyOpacity abs_em_beta;

    if (stim_abs == 1)
    {
        // i < 2 * n: we consider all the cases
        for (int i = 0; i < 2 * n; ++i) 
        {
            nu = nu_array[i];

            abs_em_beta = ElStimAbsOpacity(nu, &grey_pars->opacity_pars,
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

            abs_em_beta = ElAbsOpacity(nu, &grey_pars->opacity_pars,
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
*/


/*
KOKKOS_INLINE_FUNCTION
void AddMuonicBetaReactionToIntegrand(int n, BS_REAL* nu_array,
                                GreyOpacityParams* grey_pars,
                                M1MatrixKokkos2D* out, const int stim_abs)
{
    BS_REAL nu;
    MyOpacity abs_em_beta;

    if (stim_abs == 1)
    {
        // We consider only i<n, since we don't split the integral
        // and we have less points to be computed
        for (int i = 0; i < n; ++i)
        {
            nu = nu_array[i];

            abs_em_beta = MuonStimAbsOpacity(nu, &grey_pars->opacity_pars,
                                         &grey_pars->eos_pars); // [s^-1]

            for (int j = 0; i < n; ++j)
            {
                out->m1_mat_em[id_nux][i][j] += abs_em_beta.em[id_nux];
                out->m1_mat_em[id_anux][i][j] += abs_em_beta.em[id_anux];

                // ab = em + ab (stimulated absorption)
                out->m1_mat_ab[id_nux][i][j] += abs_em_beta.abs[id_nux];
                out->m1_mat_ab[id_anux][i][j] += abs_em_beta.abs[id_anux];
            }
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            nu = nu_array[i];

            abs_em_beta = MuonAbsOpacity(nu, &grey_pars->opacity_pars,
                                     &grey_pars->eos_pars); // [s^-1]

            for (int j = 0; i < n; ++j)
            {
                out->m1_mat_em[id_nux][i][j] += abs_em_beta.em[id_nux];
                out->m1_mat_em[id_anux][i][j] += abs_em_beta.em[id_anux];

                out->m1_mat_ab[id_nux][i][j] += abs_em_beta.abs[id_nux];
                out->m1_mat_ab[id_anux][i][j] += abs_em_beta.abs[id_anux];
            }
        }
    }
    return;
}
*/


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


/*
KOKKOS_INLINE_FUNCTION
void AddInelNEPSKernelsToIntegrand(int n, BS_REAL* nu_array,
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

        inel_1 = InelasticNEPSKernels(
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

            inel_1 = InelasticNEPSKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticNEPSKernels(
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
*/


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

        if constexpr (total_num_species == 4)
        {
            out->m1_mat_em[id_nue][i][i] *= block_factor_nu[id_anue];
            out->m1_mat_em[id_anue][i][i] *= block_factor_nu[id_nue];
            out->m1_mat_em[id_nux][i][i] *= block_factor_nu[id_anux];
            out->m1_mat_em[id_anux][i][i] *= block_factor_nu[id_nux];

            out->m1_mat_ab[id_nue][i][i] *= g_nu[id_anue];
            out->m1_mat_ab[id_anue][i][i] *= g_nu[id_nue];
            out->m1_mat_ab[id_nux][i][i] *= g_nu[id_anux];
            out->m1_mat_ab[id_anux][i][i] *= g_nu[id_nux];
        }
        else if constexpr (total_num_species == 6)
        {
            out->m1_mat_em[id_nue][i][i] *= block_factor_nu[id_anue];
            out->m1_mat_em[id_anue][i][i] *= block_factor_nu[id_nue];
            out->m1_mat_em[id_num][i][i] *= block_factor_nu[id_anum];
            out->m1_mat_em[id_anum][i][i] *= block_factor_nu[id_num];
            out->m1_mat_em[id_nut][i][i] *= block_factor_nu[id_anut];
            out->m1_mat_em[id_anut][i][i] *= block_factor_nu[id_nut];

            out->m1_mat_ab[id_nue][i][i] *= g_nu[id_anue];
            out->m1_mat_ab[id_anue][i][i] *= g_nu[id_nue];
            out->m1_mat_ab[id_num][i][i] *= g_nu[id_anum];
            out->m1_mat_ab[id_anum][i][i] *= g_nu[id_num];
            out->m1_mat_ab[id_nut][i][i] *= g_nu[id_anut];
            out->m1_mat_ab[id_anut][i][i] *= g_nu[id_nut];
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

            if constexpr (total_num_species == 4)
            {
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
            else if constexpr (total_num_species == 6)
            {
                out->m1_mat_em[id_nue][i][j] *= block_factor_nu_bar[id_anue];
                out->m1_mat_em[id_anue][i][j] *= block_factor_nu_bar[id_nue];
                out->m1_mat_em[id_num][i][j] *= block_factor_nu_bar[id_anum];
                out->m1_mat_em[id_anum][i][j] *= block_factor_nu_bar[id_num];
                out->m1_mat_em[id_nut][i][j] *= block_factor_nu_bar[id_anut];
                out->m1_mat_em[id_anut][i][j] *= block_factor_nu_bar[id_nut];

                out->m1_mat_em[id_nue][j][i] *= block_factor_nu[id_anue];
                out->m1_mat_em[id_anue][j][i] *= block_factor_nu[id_nue];
                out->m1_mat_em[id_num][j][i] *= block_factor_nu[id_anum];
                out->m1_mat_em[id_anum][j][i] *= block_factor_nu[id_num];
                out->m1_mat_em[id_nut][j][i] *= block_factor_nu[id_anut];
                out->m1_mat_em[id_anut][j][i] *= block_factor_nu[id_nut];

                out->m1_mat_ab[id_nue][i][j] *= g_nu_bar[id_anue];
                out->m1_mat_ab[id_anue][i][j] *= g_nu_bar[id_nue];
                out->m1_mat_ab[id_num][i][j] *= g_nu_bar[id_anum];
                out->m1_mat_ab[id_anum][i][j] *= g_nu_bar[id_num];
                out->m1_mat_ab[id_nut][i][j] *= g_nu_bar[id_anut];
                out->m1_mat_ab[id_anut][i][j] *= g_nu_bar[id_nut];

                out->m1_mat_ab[id_nue][j][i] *= g_nu[id_anue];
                out->m1_mat_ab[id_anue][j][i] *= g_nu[id_nue];
                out->m1_mat_ab[id_num][j][i] *= g_nu[id_anum];
                out->m1_mat_ab[id_anum][j][i] *= g_nu[id_num];
                out->m1_mat_ab[id_nut][j][i] *= g_nu[id_anut];
                out->m1_mat_ab[id_anut][j][i] *= g_nu[id_nut];
            }
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
    BS_REAL block_factor_nu[total_num_species], 
            block_factor_nu_bar[total_num_species];

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
                if (grey_pars->opacity_pars.neglect_blocking == false){
                    block_factor_nu[idx] = one - g_nu[idx];
                } else{
                    block_factor_nu[idx] = one;
                }

                out->m1_mat_ab[idx][i][i] *= nu_fourth * g_nu[idx];
                out->m1_mat_em[idx][i][i] *= nu_fourth * block_factor_nu[idx];
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
                    if (grey_pars->opacity_pars.neglect_blocking == false){
                        block_factor_nu_bar[idx] = one - g_nu_bar[idx];
                    } else{
                        block_factor_nu_bar[idx] = one;
                    }

                    out->m1_mat_ab[idx][i][j] *= nu_fourth * g_nu[idx];
                    out->m1_mat_ab[idx][j][i] *= nu_fourth * g_nu_bar[idx];

                    out->m1_mat_em[idx][i][j] *= nu_fourth * block_factor_nu[idx];
                    out->m1_mat_em[idx][j][i] *=
                        nu_fourth * block_factor_nu_bar[idx];
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
    //     AddInelNEPSKernelsToIntegrand(n, nu_array, grey_pars, &out);
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


// Compute NEPS Integrand for double GL integration
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

            // compute the NEPS kernels
            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

            inel_1 = InelasticNEPSKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticNEPSKernels(
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
                        nu_fourth * block_factor_nu[idx] * tmp_em_1;

                    out.m1_mat_ab[idx][i][n + j] =
                        nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                    out.m1_mat_em[idx][i][n + j] =
                        nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
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

            // compute the NEPS kernels
            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

            inel_1 = InelasticNEPSKernels(
                &grey_pars->kernel_pars.inelastic_kernel_params,
                &grey_pars->eos_pars);

            grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
            grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

            inel_2 = InelasticNEPSKernels(
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
                        nu_fourth * block_factor_nu[idx] * tmp_em_1;

                    out.m1_mat_ab[idx][n + i][n + j] =
                        nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                    out.m1_mat_em[idx][n + i][n + j] =
                        nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
                }
            }
        }
    }

    return out;
}


// Compute NMS Integrand for double GL integration
KOKKOS_INLINE_FUNCTION
M1MatrixKokkos2D ComputeNMSIntegrand(const MyQuadrature* quad, BS_REAL t,
                                      GreyOpacityParams* grey_pars,
                                      const int stim_abs)
{
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL one  = 1;

    BS_ASSERT((stim_abs == 0) || (stim_abs == 1));

    const int n = quad->nx;
    BS_REAL x_i, x_j;
    BS_REAL nu, nu_bar, nu_fourth;
    BS_REAL u1, u2;
    BS_REAL min1, min2;

    BS_REAL tmp_em_1, tmp_em_2;
    BS_REAL tmp_abs_1, tmp_abs_2;

    // compute the neutrino & anti-neutrino distribution function
    BS_REAL g_nu[total_num_species], g_nu_bar[total_num_species];
    BS_REAL block_factor_nu[total_num_species],
            block_factor_nu_bar[total_num_species];

    MyKernelOutput inel_1, inel_2;

    M1MatrixKokkos2D out = {0};

    if (grey_pars->opacity_pars.NMS_implementation == NMS_KernelInterp)
    {

        constexpr BS_REAL umin = NMS_w_min + NMS_wp_min;
        constexpr BS_REAL umax = NMS_w_max + NMS_wp_max;
        constexpr BS_REAL vmax = NMS_w_max - NMS_wp_min;

        for (int i = 0; i < n; ++i)
        {

            x_i = quad->points[i];

            u1 = umin + (t - umin) * x_i;
            u2 = t + (umax - t) * x_i;
            min1 = std::min(u1, vmax);
            min2 = std::min(u2, vmax);

            for (int j = 0; j < n; ++j)
            {

                x_j = quad->points[j];

                nu = half * (u1 - min1 * x_j);
                BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
                nu_bar = half * (u1 + min1 * x_j);
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

                // compute the NMS kernels
                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

                inel_1 = InelasticNMSKernels_DirectInterp(
                            &grey_pars->kernel_pars.inelastic_kernel_params,
                            &grey_pars->eos_pars);

                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

                inel_2 = InelasticNMSKernels_DirectInterp(
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
                            nu_fourth * block_factor_nu[idx] * tmp_em_1;

                        out.m1_mat_ab[idx][i][n + j] =
                            nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                        out.m1_mat_em[idx][i][n + j] =
                            nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
                    }
                }

                nu = half * (u2 - min2 * x_j);
                BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
                nu_bar = half * (u2 + min2 * x_j);
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

                // compute the NMS kernels
                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

                inel_1 = InelasticNMSKernels_DirectInterp(
                            &grey_pars->kernel_pars.inelastic_kernel_params,
                            &grey_pars->eos_pars);

                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

                inel_2 = InelasticNMSKernels_DirectInterp(
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
                            nu_fourth * block_factor_nu[idx] * tmp_em_1;

                        out.m1_mat_ab[idx][n + i][n + j] =
                            nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                        out.m1_mat_em[idx][n + i][n + j] =
                            nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
                    }
                }
            }
        }
    }

    else if (grey_pars->opacity_pars.NMS_implementation == NMS_SemiAnalytical)
    {

        constexpr BS_REAL umin = NMSParams_w_min + NMSParams_wp_min;
        constexpr BS_REAL umax = NMSParams_w_max + NMSParams_wp_max;
        constexpr BS_REAL vmax = NMSParams_w_max - NMSParams_wp_min;
    
        for (int i = 0; i < n; ++i)
        {

            x_i = quad->points[i];

            u1 = umin + (t - umin) * x_i;
            u2 = t + (umax - t) * x_i;
            min1 = std::min(u1, vmax);
            min2 = std::min(u2, vmax);

            for (int j = 0; j < n; ++j)
            {

                x_j = quad->points[j];

                nu = half * (u1 - min1 * x_j);
                BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
                nu_bar = half * (u1 + min1 * x_j);
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

                // compute the NMS kernels
                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

                inel_1 = InelasticNMSKernels_SemiAnalytical(
                            &grey_pars->kernel_pars.inelastic_kernel_params,
                            &grey_pars->eos_pars);

                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

                inel_2 = InelasticNMSKernels_SemiAnalytical(
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
                            nu_fourth * block_factor_nu[idx] * tmp_em_1;

                        out.m1_mat_ab[idx][i][n + j] =
                            nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                        out.m1_mat_em[idx][i][n + j] =
                            nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
                    }
                }

                nu = half * (u2 - min2 * x_j);
                BS_ASSERT(nu >= 0, "Neutrino energy is negative (nu=%e)", nu);
                nu_bar = half * (u2 + min2 * x_j);
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

                // compute the NMS kernels
                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;

                inel_1 = InelasticNMSKernels_SemiAnalytical(
                            &grey_pars->kernel_pars.inelastic_kernel_params,
                            &grey_pars->eos_pars);

                grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
                grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;

                inel_2 = InelasticNMSKernels_SemiAnalytical(
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
                            nu_fourth * block_factor_nu[idx] * tmp_em_1;

                        out.m1_mat_ab[idx][n + i][n + j] =
                            nu_fourth * g_nu_bar[idx] * tmp_abs_2;
                        out.m1_mat_em[idx][n + i][n + j] =
                            nu_fourth * block_factor_nu_bar[idx] * tmp_em_2;
                    }
                }
            }
        }
    }

    else  //fallback
    {
        BS_ASSERT(false, "Unknown NMS implementation: %d",
                  (int)grey_pars->opacity_pars.NMS_implementation);
    }

    return out;
}


// Compute Muon Decay Integrand for double GL integration
// (integrands for energy quantities are built directly during integration)
template <int stim_abs>  //template to speed-up the "stim_abs" if-statement
KOKKOS_INLINE_FUNCTION
M1MatrixKokkos2D ComputeMuonDecayIntegrand(const MyQuadrature* quad, BS_REAL t,
                                           GreyOpacityParams* grey_pars)
{
    const int n = quad->nx;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL max = MuonDecay_wnumu_max;
    constexpr BS_REAL min = MuonDecay_wnumu_min;
    constexpr int combinations[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    const int idx_shifter[4][2] = {{0, 0}, {0, n}, {n, 0}, {n, n}};

    BS_ASSERT((stim_abs == 0) || (stim_abs == 1));

    BS_REAL x_i, x_j;
    BS_REAL w_numu[2], w_anue[2];
    BS_REAL nu_fourth;
    BS_REAL tmp_em_numu, tmp_em_anue;
    BS_REAL tmp_abs_numu, tmp_abs_anue;

    // define neutrino & anti-neutrino distribution function
    BS_REAL g_numu[2], g_anue[2];
    BS_REAL block_factor_numu[2], block_factor_anue[2];

    MyKernelOutput decay_kernels;

    M1MatrixKokkos2D out = {0};

    for (int i = 0; i < n; ++i)
    {

        x_i = quad->points[i];

        w_numu[0] = min + (t - min) * x_i;
        w_numu[1] = t + (max - t) * x_i;
        BS_ASSERT(w_numu[0] >= 0, "Neutrino energy is negative (nu=%e)", w_numu[0]);
        BS_ASSERT(w_numu[1] >= 0, "Neutrino energy is negative (nu=%e)", w_numu[1]);

        for (int j = 0; j < n; ++j)
        {

            x_j = quad->points[j];

            w_anue[0] = min + (t - min) * x_j;
            w_anue[1] = t + (max - t) * x_j;
            BS_ASSERT(w_anue[0] >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      w_anue[0]);
            BS_ASSERT(w_anue[1] >= 0, "Neutrino energy is negative (nu_bar=%e)",
                      w_anue[1]);

            // Calculate distributions and block factors
            g_numu[0] = TotalNuF(w_numu[0], &grey_pars->distr_pars, id_num);
            g_numu[1] = TotalNuF(w_numu[1], &grey_pars->distr_pars, id_num);
            g_anue[0] = TotalNuF(w_anue[0], &grey_pars->distr_pars, id_anue);
            g_anue[1] = TotalNuF(w_anue[1], &grey_pars->distr_pars, id_anue);

            if (grey_pars->opacity_pars.neglect_blocking == false)
            {
                block_factor_numu[0] = one - g_numu[0];
                block_factor_numu[1] = one - g_numu[1];
                block_factor_anue[0] = one - g_anue[0];
                block_factor_anue[1] = one - g_anue[1];
            }
            else
            {
                block_factor_numu[0] = one;
                block_factor_numu[1] = one;
                block_factor_anue[0] = one;
                block_factor_anue[1] = one;
            }
            
            // Compute the Muon Decay Integrands:

            // Cycle over all combinations:
            for (int k = 0; k < 4; ++k)
            {
                int idx_numu = combinations[k][0];
                int idx_anue = combinations[k][1];
                int N_i = idx_shifter[k][0];
                int N_j = idx_shifter[k][1];

                grey_pars->kernel_pars.muon_decay_kernel_params.omega_numu = w_numu[idx_numu];
                grey_pars->kernel_pars.muon_decay_kernel_params.omega_anue = w_anue[idx_anue];

                decay_kernels = MuonDecayKernels(&grey_pars->kernel_pars.muon_decay_kernel_params,
                                                 &grey_pars->eos_pars);
                nu_fourth = POW2(w_numu[idx_numu]) * POW2(w_anue[idx_anue]);

                // Take always the id_num case of the kernel
                tmp_em_numu  = decay_kernels.em[id_num] * block_factor_anue[idx_anue];
                tmp_abs_numu = decay_kernels.abs[id_num] * g_anue[idx_anue];
                tmp_em_anue  = decay_kernels.em[id_num] * block_factor_numu[idx_numu];
                tmp_abs_anue = decay_kernels.abs[id_num] * g_numu[idx_numu];

                if constexpr (stim_abs == 1)
                {
                    out.m1_mat_ab[id_num][i + N_i][j + N_j] =
                        nu_fourth * g_numu[idx_numu] * (tmp_em_numu + tmp_abs_numu);
                    out.m1_mat_em[id_num][i + N_i][j + N_j] = nu_fourth * tmp_em_numu;
                    out.m1_mat_ab[id_anue][i + N_i][j + N_j] =
                        nu_fourth * g_anue[idx_anue] * (tmp_em_anue + tmp_abs_anue);
                    out.m1_mat_em[id_anue][i + N_i][j + N_j] = nu_fourth * tmp_em_anue;
                }
                else
                {
                    out.m1_mat_ab[id_num][i + N_i][j + N_j] =
                        nu_fourth * g_numu[idx_numu] * tmp_abs_numu;
                    out.m1_mat_em[id_num][i + N_i][j + N_j] =
                        nu_fourth * block_factor_numu[idx_numu] * tmp_em_numu;
                    out.m1_mat_ab[id_anue][i + N_i][j + N_j] = 
                                                    out.m1_mat_ab[id_num][i + N_i][j + N_j];
                    out.m1_mat_em[id_anue][i + N_i][j + N_j] = 
                                                    out.m1_mat_em[id_num][i + N_i][j + N_j];
                }
            }
        }
    }
    return out;
}


/* Computes the opacities for the M1 code, with thermal and
 * non-thermal processes treated together.
 *
 * NEPS and NMS are treated together with other reactions.
 * NEPS and NMS are considered for computation of number emissivity (eta_0) 
 * and absorsivity (kappa_0_a).
 */
KOKKOS_INLINE_FUNCTION
M1Opacities ComputeM1OpacitiesGenericFormalism(
    const MyQuadrature* quad_1d, const MyQuadrature* quad_2d,
    GreyOpacityParams* my_grey_opacity_params, const int stim_abs)
{

    constexpr int id_nue_kappa = total_num_species + id_nue;
    constexpr int id_anue_kappa = total_num_species + id_anue;
    constexpr int id_nux_kappa = total_num_species + id_nux;
    constexpr int id_anux_kappa = total_num_species + id_anux;
    constexpr int id_num_kappa = total_num_species + id_num;
    constexpr int id_anum_kappa = total_num_species + id_anum;
    constexpr int id_nut_kappa = total_num_species + id_nut;
    constexpr int id_anut_kappa = total_num_species + id_anut;

    constexpr BS_REAL four    = 4;
    constexpr BS_REAL two     = 2;
    constexpr BS_REAL c_light = kBS_Clight;
    BS_REAL umin, umax, vmax;
    BS_REAL dQ = kBS_Q;
    BS_REAL dU = 0.;
    // Mean field corrections
    if (my_grey_opacity_params->opacity_pars.use_dU)
        dU = my_grey_opacity_params->eos_pars.dU; // [MeV]
    if (my_grey_opacity_params->opacity_pars.use_dm_eff)
        dQ = my_grey_opacity_params->eos_pars.dm_eff; // [MeV]
    const BS_REAL Q = dU + dQ;

    // Extremals for NMS integration
    if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_KernelInterp)
    {
        umin = NMS_w_min + NMS_wp_min;
        umax = NMS_w_max + NMS_wp_max;
        vmax = NMS_w_max - NMS_wp_min;
    }
    else if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_SemiAnalytical)
    {
        umin = NMSParams_w_min + NMSParams_wp_min;
        umax = NMSParams_w_max + NMSParams_wp_max;
        vmax = NMSParams_w_max - NMSParams_wp_min;
    }


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
    const BS_REAL s_nux  = three_halves * temp;
    const BS_REAL s_neps = temp_multiple * temp;
    constexpr BS_REAL s_mudec = kBS_Mmu / 3.; 
    // s_nms depends on the grid boundaries:
    const BS_REAL s_nms = std::max(two * umin, std::min(four * s_neps, five_sixths * umax));

    BS_REAL s_beta_el[total_num_species] = {0}, s_beta_muon[total_num_species] = {0},
            s_iso[total_num_species] = {0};

    if (eta_e > -30. and eta_e < 30.)
    {
        s_beta_el[id_nue]  = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
        s_beta_el[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    }
    else if (eta_e > 30.)
    {
        s_beta_el[id_nue]  = temp * eta_e * five_sixths;
        s_beta_el[id_anue] = temp * five;
    }
    else
    {
        s_beta_el[id_nue]  = temp * five;
        s_beta_el[id_anue] = temp * eta_e * five_sixths;
    }

    s_beta_muon[id_num] = two * kBS_Mmu - Q;
    s_beta_muon[id_anum] = two * kBS_Mmu + Q;

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


    MyQuadratureIntegrand el_beta_n_em_integrals  = {0};
    MyQuadratureIntegrand el_beta_j_em_integrals  = {0};
    MyQuadratureIntegrand el_beta_n_abs_integrals = {0};
    MyQuadratureIntegrand el_beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_abs_em == 1)
    {
        BS_REAL out_el_beta_em[total_num_species][BS_N_MAX];
        BS_REAL out_el_beta_ab[total_num_species][BS_N_MAX];

        ElBeta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta_el, out_el_beta_em,
                        out_el_beta_ab, stim_abs);

        ElBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_el_beta_em,
                                                 s_beta_el, &el_beta_n_em_integrals,
                                                 &el_beta_j_em_integrals);
        ElBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_el_beta_ab,
                                                 s_beta_el, &el_beta_n_abs_integrals,
                                                 &el_beta_j_abs_integrals);
    }


    MyQuadratureIntegrand muon_beta_n_em_integrals  = {0};
    MyQuadratureIntegrand muon_beta_j_em_integrals  = {0};
    MyQuadratureIntegrand muon_beta_n_abs_integrals = {0};
    MyQuadratureIntegrand muon_beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_muonic_beta == 1)
    {
        BS_REAL out_muon_beta_em[total_num_species][BS_N_MAX];
        BS_REAL out_muon_beta_ab[total_num_species][BS_N_MAX];

        MuonicBeta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta_muon, Q, out_muon_beta_em,
                            out_muon_beta_ab, stim_abs);

        MuonicBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, out_muon_beta_em,
                                                 s_beta_muon, Q, &muon_beta_n_em_integrals,
                                                 &muon_beta_j_em_integrals);
        MuonicBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, out_muon_beta_ab,
                                                 s_beta_muon, Q, &muon_beta_n_abs_integrals,
                                                 &muon_beta_j_abs_integrals);
    }


    MyQuadratureIntegrand n_integrals_2d = {0};
    MyQuadratureIntegrand e_integrals_2d = {0};

    if ((my_grey_opacity_params->opacity_flags.use_pair == 1) ||
        (my_grey_opacity_params->opacity_flags.use_brem == 1))
    {
        M1MatrixKokkos2D out_pair = ComputeDoubleIntegrand(
            quad_2d, s_pair, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForM1Coeffs(
            quad_2d, &out_pair, s_pair, &n_integrals_2d, &e_integrals_2d);
    }


    MyQuadratureIntegrand n_neps_2d = {0};
    MyQuadratureIntegrand e_neps_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_NEPS == 1)
    {
        M1MatrixKokkos2D out_inel = ComputeNEPSIntegrand(
            quad_2d, four * s_neps, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForNEPS(quad_2d, &out_inel, four * s_neps,
                                              &n_neps_2d, &e_neps_2d);
    }


    MyQuadratureIntegrand n_nms_2d = {0};
    MyQuadratureIntegrand e_nms_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_NMS == 1)
    {
        M1MatrixKokkos2D out_inel = ComputeNMSIntegrand(
            quad_2d, s_nms, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForNMS(quad_2d, &out_inel, s_nms,
                                              umin, umax, vmax,
                                              &n_nms_2d, &e_nms_2d);
    }


    MyQuadratureIntegrand n_mudec_2d = {0};
    MyQuadratureIntegrand e_mudec_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_muon_decay == 1)
    {
        M1MatrixKokkos2D out_mudec;

        // use templates to speed-up the "stim_abs" if-statement
        if (stim_abs == 1){
            out_mudec = ComputeMuonDecayIntegrand<1>(
                            quad_2d, s_mudec, my_grey_opacity_params);
        } else {
            out_mudec = ComputeMuonDecayIntegrand<0>(
                            quad_2d, s_mudec, my_grey_opacity_params);
        }

        GaussLegendreIntegrate2DMatrixForMuonDecay(quad_2d, &out_mudec,
                                                   &n_mudec_2d, &e_mudec_2d);
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

    if constexpr (total_num_species == 4)
    {

        /* Electron neutrinos */
        m1_opacities.eta_0[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_nue] +
                                                n_neps_2d.integrand[id_nue]) +
                            el_beta_n_em_integrals.integrand[id_nue]);
        m1_opacities.eta[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_nue] +
                                                e_neps_2d.integrand[id_nue]) +
                            el_beta_j_em_integrals.integrand[id_nue]);
        if (n[id_nue] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_nue] =
                kBS_FourPi_hc3 / (c_light * n[id_nue]) *
                (kBS_FourPi_hc3 *
                    (n_integrals_2d.integrand[id_nue_kappa] + 
                     n_neps_2d.integrand[id_nue_kappa]) +
                el_beta_n_abs_integrals.integrand[id_nue]);
        }
        if (J[id_nue] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                        (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_nue_kappa] +
                                           e_neps_2d.integrand[id_nue_kappa]) +
                        el_beta_j_abs_integrals.integrand[id_nue]);
            m1_opacities.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                        iso_integrals.integrand[id_nue];
        }

        /* Electron anti-neutrinos */
        m1_opacities.eta_0[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anue] +
                                                n_neps_2d.integrand[id_anue]) +
                            el_beta_n_em_integrals.integrand[id_anue]);
        m1_opacities.eta[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anue] +
                                                e_neps_2d.integrand[id_anue]) +
                            el_beta_j_em_integrals.integrand[id_anue]);
        if (n[id_anue] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * n[id_anue]) *
                (kBS_FourPi_hc3 *
                    (n_integrals_2d.integrand[id_anue_kappa] + 
                     n_neps_2d.integrand[id_anue_kappa]) +
                el_beta_n_abs_integrals.integrand[id_anue]);
        }
        if (J[id_anue] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * J[id_anue]) *
                (kBS_FourPi_hc3 *
                    (e_integrals_2d.integrand[id_anue_kappa] + 
                     e_neps_2d.integrand[id_anue_kappa]) +
                el_beta_j_abs_integrals.integrand[id_anue]);
            m1_opacities.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anue]) *
                                            iso_integrals.integrand[id_anue];
        }

        /* Heavy neutrinos */
        m1_opacities.eta_0[id_nux] =
            kBS_FourPi_hc3_sqr * (n_integrals_2d.integrand[id_nux] + 
                                  n_neps_2d.integrand[id_nux]);
        m1_opacities.eta[id_nux] =
            kBS_FourPi_hc3_sqr * (e_integrals_2d.integrand[id_nux] + 
                                  e_neps_2d.integrand[id_nux]);
        if (n[id_nux] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_nux] =
                kBS_FourPi_hc3_sqr / (c_light * n[id_nux]) *
                (n_integrals_2d.integrand[id_nux_kappa] + 
                 n_neps_2d.integrand[id_nux_kappa]);
        }
        if (J[id_nux] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_nux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
                (e_integrals_2d.integrand[id_nux_kappa] + 
                 e_neps_2d.integrand[id_nux_kappa]);
            m1_opacities.kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J[id_nux]) *
                                        iso_integrals.integrand[id_nux];
        }

        /* Heavy anti-neutrinos */
        m1_opacities.eta_0[id_anux] =
            kBS_FourPi_hc3_sqr * (n_integrals_2d.integrand[id_anux] + 
                                  n_neps_2d.integrand[id_anux]);
        m1_opacities.eta[id_anux] =
            kBS_FourPi_hc3_sqr * (e_integrals_2d.integrand[id_anux] + 
                                  e_neps_2d.integrand[id_anux]);
        if (n[id_anux] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_anux] =
                n[id_anux] == zero ?
                    zero :
                    kBS_FourPi_hc3_sqr / (c_light * n[id_anux]) *
                    (n_integrals_2d.integrand[id_anux_kappa] + 
                     n_neps_2d.integrand[id_anux_kappa]);
        }
        if (J[id_anux] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_anux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
                (e_integrals_2d.integrand[id_anux_kappa] + 
                 e_neps_2d.integrand[id_anux_kappa]);
            m1_opacities.kappa_s[id_anux] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anux]) *
                                            iso_integrals.integrand[id_anux];
        }
    }

    else if constexpr (total_num_species == 6)
    {

        /* Electron neutrinos */
        m1_opacities.eta_0[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_nue] +
                                                n_neps_2d.integrand[id_nue] +
                                                n_nms_2d.integrand[id_nue]) +
                              el_beta_n_em_integrals.integrand[id_nue]);
        m1_opacities.eta[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_nue] +
                                                e_neps_2d.integrand[id_nue] +
                                                e_nms_2d.integrand[id_nue]) +
                              el_beta_j_em_integrals.integrand[id_nue]);
        if (n[id_nue] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_nue] =
                kBS_FourPi_hc3 / (c_light * n[id_nue]) *
                (kBS_FourPi_hc3 *
                    (n_integrals_2d.integrand[id_nue_kappa] + 
                    n_neps_2d.integrand[id_nue_kappa] + 
                    n_nms_2d.integrand[id_nue_kappa]) +
                 el_beta_n_abs_integrals.integrand[id_nue]);
        }
        if (J[id_nue] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                        (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_nue_kappa] +
                                        e_neps_2d.integrand[id_nue_kappa] + 
                                        e_nms_2d.integrand[id_nue_kappa]) +
                         el_beta_j_abs_integrals.integrand[id_nue]);
            m1_opacities.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                        iso_integrals.integrand[id_nue];
        }

        /* Electron anti-neutrinos */
        m1_opacities.eta_0[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anue] +
                                                n_neps_2d.integrand[id_anue] +
                                                n_nms_2d.integrand[id_anue] +
                                                n_mudec_2d.integrand[id_anue]) +
                              el_beta_n_em_integrals.integrand[id_anue]);
        m1_opacities.eta[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anue] +
                                                e_neps_2d.integrand[id_anue] +
                                                e_nms_2d.integrand[id_anue] +
                                                e_mudec_2d.integrand[id_anue]) +
                              el_beta_j_em_integrals.integrand[id_anue]);
        if (n[id_anue] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * n[id_anue]) *
                (kBS_FourPi_hc3 *
                    (n_integrals_2d.integrand[id_anue_kappa] + 
                    n_neps_2d.integrand[id_anue_kappa] +
                    n_nms_2d.integrand[id_anue_kappa] +
                    n_mudec_2d.integrand[id_anue_kappa]) +
                 el_beta_n_abs_integrals.integrand[id_anue]);
        }
        if (J[id_anue] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * J[id_anue]) *
                (kBS_FourPi_hc3 *
                    (e_integrals_2d.integrand[id_anue_kappa] + 
                    e_neps_2d.integrand[id_anue_kappa] +
                    e_nms_2d.integrand[id_anue_kappa] +
                    e_mudec_2d.integrand[id_anue_kappa]) +
                 el_beta_j_abs_integrals.integrand[id_anue]);
            m1_opacities.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anue]) *
                                            iso_integrals.integrand[id_anue];
        }

        /* Muon neutrinos */
        m1_opacities.eta_0[id_num] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_num] + 
                                                n_neps_2d.integrand[id_num] +
                                                n_nms_2d.integrand[id_num] +
                                                n_mudec_2d.integrand[id_num]) +
                              muon_beta_n_em_integrals.integrand[id_num]);
        m1_opacities.eta[id_num] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_num] + 
                                                e_neps_2d.integrand[id_num] +
                                                e_nms_2d.integrand[id_num] +
                                                e_mudec_2d.integrand[id_num]) +
                              muon_beta_j_em_integrals.integrand[id_num]);
        if (n[id_num] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_num] =
                kBS_FourPi_hc3 / (c_light * n[id_num]) *
                (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_num_kappa] + 
                                n_neps_2d.integrand[id_num_kappa] +
                                n_nms_2d.integrand[id_num_kappa] +
                                n_mudec_2d.integrand[id_num_kappa]) +
                 muon_beta_n_abs_integrals.integrand[id_num]);
        }
        if (J[id_num] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_num] =
                kBS_FourPi_hc3 / (c_light * J[id_num]) *
                (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_num_kappa] + 
                                e_neps_2d.integrand[id_num_kappa] +
                                e_nms_2d.integrand[id_num_kappa] +
                                e_mudec_2d.integrand[id_num_kappa]) +
                 muon_beta_j_abs_integrals.integrand[id_num]);
            m1_opacities.kappa_s[id_num] = kBS_FourPi_hc3 / (c_light * J[id_num]) *
                                        iso_integrals.integrand[id_num];
        }

        /* Muon anti-neutrinos */
        m1_opacities.eta_0[id_anum] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anum] + 
                                                n_neps_2d.integrand[id_anum] +
                                                n_nms_2d.integrand[id_anum]) + 
                              muon_beta_n_em_integrals.integrand[id_anum]);
        m1_opacities.eta[id_anum] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anum] + 
                                                e_neps_2d.integrand[id_anum] +
                                                e_nms_2d.integrand[id_anum]) + 
                              muon_beta_j_em_integrals.integrand[id_anum]);
        if (n[id_anum] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_anum] =
                n[id_anum] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * n[id_anum]) *
                        (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anum_kappa] + 
                                        n_neps_2d.integrand[id_anum_kappa] +
                                        n_nms_2d.integrand[id_anum_kappa]) +
                         muon_beta_n_abs_integrals.integrand[id_anum]);
        }
        if (J[id_anum] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_anum] =
                kBS_FourPi_hc3 / (c_light * J[id_anum]) *
                (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anum_kappa] + 
                                e_neps_2d.integrand[id_anum_kappa] +
                                e_nms_2d.integrand[id_anum_kappa]) +
                 muon_beta_j_abs_integrals.integrand[id_anum]);

            m1_opacities.kappa_s[id_anum] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anum]) *
                                            iso_integrals.integrand[id_anum];
        }

        /* Tau neutrinos */
        m1_opacities.eta_0[id_nut] =
            kBS_FourPi_hc3_sqr * (n_integrals_2d.integrand[id_nut] + 
                                  n_neps_2d.integrand[id_nut] +
                                  n_nms_2d.integrand[id_nut]);
        m1_opacities.eta[id_nut] =
            kBS_FourPi_hc3_sqr * (e_integrals_2d.integrand[id_nut] + 
                                  e_neps_2d.integrand[id_nut] +
                                  e_nms_2d.integrand[id_nut]);
        if (n[id_nut] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_nut] =
                kBS_FourPi_hc3_sqr / (c_light * n[id_nut]) *
                (n_integrals_2d.integrand[id_nut_kappa] + 
                 n_neps_2d.integrand[id_nut_kappa] +
                 n_nms_2d.integrand[id_nut_kappa]);
        }
        if (J[id_nut] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_nut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nut]) *
                (e_integrals_2d.integrand[id_nut_kappa] + 
                 e_neps_2d.integrand[id_nut_kappa] +
                 e_nms_2d.integrand[id_nut_kappa]);
            m1_opacities.kappa_s[id_nut] = kBS_FourPi_hc3 / (c_light * J[id_nut]) *
                                        iso_integrals.integrand[id_nut];
        }

        /* Tau anti-neutrinos */
        m1_opacities.eta_0[id_anut] =
            kBS_FourPi_hc3_sqr * (n_integrals_2d.integrand[id_anut] + 
                                  n_neps_2d.integrand[id_anut] +
                                  n_nms_2d.integrand[id_anut]);
        m1_opacities.eta[id_anut] =
            kBS_FourPi_hc3_sqr * (e_integrals_2d.integrand[id_anut] + 
                                  e_neps_2d.integrand[id_anut] +
                                  e_nms_2d.integrand[id_anut]);
        if (n[id_anut] > THRESHOLD_N)
        {
            m1_opacities.kappa_0_a[id_anut] =
                n[id_anut] == zero ?
                    zero :
                    kBS_FourPi_hc3_sqr / (c_light * n[id_anut]) *
                        (n_integrals_2d.integrand[id_anut_kappa] + 
                         n_neps_2d.integrand[id_anut_kappa] +
                         n_nms_2d.integrand[id_anut_kappa]);
        }
        if (J[id_anut] > THRESHOLD_J)
        {
            m1_opacities.kappa_a[id_anut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anut]) *
                (e_integrals_2d.integrand[id_anut_kappa] + 
                 e_neps_2d.integrand[id_anut_kappa] +
                 e_nms_2d.integrand[id_anut_kappa]);

            m1_opacities.kappa_s[id_anut] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anut]) *
                                            iso_integrals.integrand[id_anut];
        }
    }

    else  //fallback
    {
        BS_ASSERT(false, "Wrong number of nu species: %d", total_num_species);
    }

    return m1_opacities;
}


/* Computes the opacities for the M1 code, with thermal and
 * non-thermal processes treated separately.
 *
 * NEPS and NMS are treated separately from other reactions.
 * NEPS and NMS are NOT considered for computation of number emissivity (eta_0) 
 * and absorsivity (kappa_0_a).
 */
KOKKOS_INLINE_FUNCTION
M1OpacitiesNonThermalSeparated ComputeM1OpacitiesGenericFormalismNonThermalSeparated(
                const MyQuadrature* quad_1d, const MyQuadrature* quad_2d,
                GreyOpacityParams* my_grey_opacity_params, const int stim_abs)
{

    constexpr int id_nue_kappa = total_num_species + id_nue;
    constexpr int id_anue_kappa = total_num_species + id_anue;
    constexpr int id_nux_kappa = total_num_species + id_nux;
    constexpr int id_anux_kappa = total_num_species + id_anux;
    constexpr int id_num_kappa = total_num_species + id_num;
    constexpr int id_anum_kappa = total_num_species + id_anum;
    constexpr int id_nut_kappa = total_num_species + id_nut;
    constexpr int id_anut_kappa = total_num_species + id_anut;

    constexpr BS_REAL four = 4;
    constexpr BS_REAL two = 2;
    constexpr BS_REAL c_light = kBS_Clight;
    BS_REAL umin, umax, vmax;
    BS_REAL dQ = kBS_Q;
    BS_REAL dU = 0.;
    // Mean field corrections
    if (my_grey_opacity_params->opacity_pars.use_dU)
        dU = my_grey_opacity_params->eos_pars.dU; // [MeV]
    if (my_grey_opacity_params->opacity_pars.use_dm_eff)
        dQ = my_grey_opacity_params->eos_pars.dm_eff; // [MeV]
    const BS_REAL Q = dU + dQ;

    // Extremals for NMS integration
    if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_KernelInterp)
    {
        umin = NMS_w_min + NMS_wp_min;
        umax = NMS_w_max + NMS_wp_max;
        vmax = NMS_w_max - NMS_wp_min;
    }
    else if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_SemiAnalytical)
    {
        umin = NMSParams_w_min + NMSParams_wp_min;
        umax = NMSParams_w_max + NMSParams_wp_max;
        vmax = NMSParams_w_max - NMSParams_wp_min;
    }

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
    const BS_REAL s_nux  = three_halves * temp;
    const BS_REAL s_neps = temp_multiple * temp;
    constexpr BS_REAL s_mudec = kBS_Mmu / 3.;
    // s_nms depends on the grid boundaries:
    const BS_REAL s_nms = std::max(two * umin, std::min(four * s_neps, five_sixths * umax));

    BS_REAL s_beta_el[total_num_species] = {0}, s_beta_muon[total_num_species];
    BS_REAL s_iso[total_num_species] = {0};

    if (eta_e > -30. and eta_e < 30.)
    {
        s_beta_el[id_nue]  = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
        s_beta_el[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    }
    else if (eta_e > 30.)
    {
        s_beta_el[id_nue]  = temp * eta_e * five_sixths;
        s_beta_el[id_anue] = temp * five;
    }
    else
    {
        s_beta_el[id_nue]  = temp * five;
        s_beta_el[id_anue] = temp * eta_e * five_sixths;
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        s_iso[idx] = (n[idx] > THRESHOLD_N) ? (J[idx] / n[idx]) : s_nux;
    }

    s_beta_muon[id_num] = two * kBS_Mmu - Q;
    s_beta_muon[id_anum] = two * kBS_Mmu + Q;

    MyQuadratureIntegrand iso_integrals = {0};
    if (my_grey_opacity_params->opacity_flags.use_iso == 1)
    {
        BS_REAL out_iso[total_num_species][BS_N_MAX];
        Scattering1DIntegrand(quad_1d, my_grey_opacity_params, s_iso, out_iso);
        iso_integrals = GaussLegendreIntegrate1DMatrix(
            quad_1d, total_num_species, out_iso, s_iso);
    }

    MyQuadratureIntegrand el_beta_n_em_integrals  = {0};
    MyQuadratureIntegrand el_beta_j_em_integrals  = {0};
    MyQuadratureIntegrand el_beta_n_abs_integrals = {0};
    MyQuadratureIntegrand el_beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_abs_em == 1)
    {
        BS_REAL out_el_beta_em[total_num_species][BS_N_MAX];
        BS_REAL out_el_beta_ab[total_num_species][BS_N_MAX];

        ElBeta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta_el, out_el_beta_em,
                        out_el_beta_ab, stim_abs);

        ElBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_el_beta_em,
                                                 s_beta_el, &el_beta_n_em_integrals,
                                                 &el_beta_j_em_integrals);
        ElBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_el_beta_ab,
                                                 s_beta_el, &el_beta_n_abs_integrals,
                                                 &el_beta_j_abs_integrals);
    }


    MyQuadratureIntegrand muon_beta_n_em_integrals  = {0};
    MyQuadratureIntegrand muon_beta_j_em_integrals  = {0};
    MyQuadratureIntegrand muon_beta_n_abs_integrals = {0};
    MyQuadratureIntegrand muon_beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_muonic_beta == 1)
    {
        BS_REAL out_muon_beta_em[total_num_species][BS_N_MAX];
        BS_REAL out_muon_beta_ab[total_num_species][BS_N_MAX];

        MuonicBeta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta_muon, Q, out_muon_beta_em,
                            out_muon_beta_ab, stim_abs);

        MuonicBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, out_muon_beta_em,
                                                 s_beta_muon, Q, &muon_beta_n_em_integrals,
                                                 &muon_beta_j_em_integrals);
        MuonicBetaGaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, out_muon_beta_ab,
                                                 s_beta_muon, Q, &muon_beta_n_abs_integrals,
                                                 &muon_beta_j_abs_integrals);
    }


    MyQuadratureIntegrand n_integrals_2d = {0};
    MyQuadratureIntegrand e_integrals_2d = {0};

    if ((my_grey_opacity_params->opacity_flags.use_pair == 1) ||
        (my_grey_opacity_params->opacity_flags.use_brem == 1))
    {
        M1MatrixKokkos2D out_pair = ComputeDoubleIntegrand(
            quad_2d, s_pair, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForM1Coeffs(
            quad_2d, &out_pair, s_pair, &n_integrals_2d, &e_integrals_2d);
    }

    MyQuadratureIntegrand n_neps_2d = {0};
    MyQuadratureIntegrand e_neps_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_NEPS == 1)
    {
        M1MatrixKokkos2D out_inel = ComputeNEPSIntegrand(
            quad_2d, four * s_neps, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForNEPS(quad_2d, &out_inel, four * s_neps,
                                              &n_neps_2d, &e_neps_2d);
    }

    MyQuadratureIntegrand n_nms_2d = {0};
    MyQuadratureIntegrand e_nms_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_inelastic_NMS == 1)
    {
        M1MatrixKokkos2D out_inel = ComputeNMSIntegrand(
            quad_2d, s_nms, my_grey_opacity_params, stim_abs);
        GaussLegendreIntegrate2DMatrixForNMS(quad_2d, &out_inel, s_nms,
                                              umin, umax, vmax,
                                              &n_nms_2d, &e_nms_2d);
    }


    MyQuadratureIntegrand n_mudec_2d = {0};
    MyQuadratureIntegrand e_mudec_2d = {0};

    if (my_grey_opacity_params->opacity_flags.use_muon_decay == 1)
    {
        M1MatrixKokkos2D out_mudec;

        // use templates to speed-up the "stim_abs" if-statement
        if (stim_abs == 1){
            out_mudec = ComputeMuonDecayIntegrand<1>(
                            quad_2d, s_mudec, my_grey_opacity_params);
        } else {
            out_mudec = ComputeMuonDecayIntegrand<0>(
                            quad_2d, s_mudec, my_grey_opacity_params);
        }

        GaussLegendreIntegrate2DMatrixForMuonDecay(quad_2d, &out_mudec,
                                                   &n_mudec_2d, &e_mudec_2d);
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


    if constexpr (total_num_species == 4)
    {
        /* Electron neutrinos */
        m1_opacities_non_th_separated.eta_0[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_nue] +
                              el_beta_n_em_integrals.integrand[id_nue]);

        m1_opacities_non_th_separated.eta_th[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_nue] +
                              el_beta_j_em_integrals.integrand[id_nue]);

        m1_opacities_non_th_separated.eta_non_th[id_nue] =
            kBS_FourPi_hc3_sqr * e_neps_2d.integrand[id_nue];

        if (n[id_nue] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_nue] =
                kBS_FourPi_hc3 / (c_light * n[id_nue]) *
                (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_nue_kappa] +
                 el_beta_n_abs_integrals.integrand[id_nue]);
        }
        if (J[id_nue] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                        (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_nue_kappa] +
                         el_beta_j_abs_integrals.integrand[id_nue]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3_sqr / (c_light * J[id_nue]) * e_neps_2d.integrand[id_nue_kappa];

            m1_opacities_non_th_separated.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                        iso_integrals.integrand[id_nue];
        }

        /* Electron anti-neutrinos */
        m1_opacities_non_th_separated.eta_0[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_anue] +
                              el_beta_n_em_integrals.integrand[id_anue]);

        m1_opacities_non_th_separated.eta_th[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_anue] +
                              el_beta_j_em_integrals.integrand[id_anue]);

        m1_opacities_non_th_separated.eta_non_th[id_anue] =
            kBS_FourPi_hc3_sqr * e_neps_2d.integrand[id_anue];

        if (n[id_anue] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * n[id_anue]) *
                (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_anue_kappa] +
                 el_beta_n_abs_integrals.integrand[id_anue]);
        }
        if (J[id_anue] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_anue] =
                kBS_FourPi_hc3 / (c_light * J[id_anue]) *
                (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_anue_kappa] +
                 el_beta_j_abs_integrals.integrand[id_anue]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_anue] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anue]) * e_neps_2d.integrand[id_anue_kappa];

            m1_opacities_non_th_separated.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anue]) *
                                            iso_integrals.integrand[id_anue];
        }

        /* Heavy neutrinos */
        m1_opacities_non_th_separated.eta_0[id_nux] =
            kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[id_nux];

        m1_opacities_non_th_separated.eta_th[id_nux] =
            kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[id_nux];

        m1_opacities_non_th_separated.eta_non_th[id_nux] =
            kBS_FourPi_hc3_sqr * e_neps_2d.integrand[id_nux];

        if (n[id_nux] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_nux] =
                kBS_FourPi_hc3_sqr / (c_light * n[id_nux]) *
                n_integrals_2d.integrand[id_nux_kappa];
        }
        if (J[id_nux] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_nux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
                e_integrals_2d.integrand[id_nux_kappa];

            m1_opacities_non_th_separated.kappa_a_non_th[id_nux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nux]) *
                e_neps_2d.integrand[6];

            m1_opacities_non_th_separated.kappa_s[id_nux] = kBS_FourPi_hc3 / (c_light * J[id_nux]) *
                                        iso_integrals.integrand[id_nux];
        }

        /* Heavy anti-neutrinos */
        m1_opacities_non_th_separated.eta_0[id_anux] =
            kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[id_anux];

        m1_opacities_non_th_separated.eta_th[id_anux] =
            kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[id_anux];

        m1_opacities_non_th_separated.eta_non_th[id_anux] =
            kBS_FourPi_hc3_sqr * e_neps_2d.integrand[id_anux];

        if (n[id_anux] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_anux] =
                n[id_anux] == zero ?
                    zero :
                    kBS_FourPi_hc3_sqr / (c_light * n[id_anux]) *
                    n_integrals_2d.integrand[id_anux_kappa];
        }
        if (J[id_anux] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_anux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
                e_integrals_2d.integrand[id_anux_kappa];

            m1_opacities_non_th_separated.kappa_a_non_th[id_anux] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anux]) *
                e_neps_2d.integrand[id_anux_kappa];

            m1_opacities_non_th_separated.kappa_s[id_anux] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anux]) *
                                            iso_integrals.integrand[id_anux];
        }

    }

    else if constexpr (total_num_species == 6)
    {
        /* Electron neutrinos */
        m1_opacities_non_th_separated.eta_0[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_nue] +
                              el_beta_n_em_integrals.integrand[id_nue]);

        m1_opacities_non_th_separated.eta_th[id_nue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_nue] +
                              el_beta_j_em_integrals.integrand[id_nue]);

        m1_opacities_non_th_separated.eta_non_th[id_nue] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_nue] +
                                  e_nms_2d.integrand[id_nue]);

        if (n[id_nue] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_nue] =
                kBS_FourPi_hc3 / (c_light * n[id_nue]) *
                (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_nue_kappa] +
                 el_beta_n_abs_integrals.integrand[id_nue]);
        }
        if (J[id_nue] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                        (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_nue_kappa] +
                         el_beta_j_abs_integrals.integrand[id_nue]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_nue] =
                n[id_nue] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                        (kBS_FourPi_hc3 * (e_neps_2d.integrand[id_nue_kappa] +
                                           e_nms_2d.integrand[id_nue_kappa]));

            m1_opacities_non_th_separated.kappa_s[id_nue] = kBS_FourPi_hc3 / (c_light * J[id_nue]) *
                                        iso_integrals.integrand[id_nue];
        }

        /* Electron anti-neutrinos */
        m1_opacities_non_th_separated.eta_0[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anue] +
                                                n_mudec_2d.integrand[id_anue]) +
                              el_beta_n_em_integrals.integrand[id_anue]);

        m1_opacities_non_th_separated.eta_th[id_anue] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anue] +
                                                e_mudec_2d.integrand[id_anue]) +
                              el_beta_j_em_integrals.integrand[id_anue]);

        m1_opacities_non_th_separated.eta_non_th[id_anue] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_anue] +
                                  e_nms_2d.integrand[id_anue]);

        if (n[id_anue] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_anue] =
                kBS_FourPi_hc3 / (c_light * n[id_anue]) *
                (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_anue_kappa] +
                                   n_mudec_2d.integrand[id_anue_kappa]) +
                 el_beta_n_abs_integrals.integrand[id_anue]);
        }
        if (J[id_anue] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_anue] =
                kBS_FourPi_hc3 / (c_light * J[id_anue]) *
                (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_anue_kappa] +
                                   e_mudec_2d.integrand[id_anue_kappa]) +
                 el_beta_j_abs_integrals.integrand[id_anue]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_anue] =
                kBS_FourPi_hc3 / (c_light * J[id_anue]) *
                (kBS_FourPi_hc3 * (e_neps_2d.integrand[id_anue_kappa] +
                                   e_nms_2d.integrand[id_anue_kappa]));

            m1_opacities_non_th_separated.kappa_s[id_anue] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anue]) *
                                            iso_integrals.integrand[id_anue];
        }

        /* Muon neutrinos */
        m1_opacities_non_th_separated.eta_0[id_num] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_num] +
                                                n_mudec_2d.integrand[id_num]) +
                              muon_beta_n_em_integrals.integrand[id_num]);

        m1_opacities_non_th_separated.eta_th[id_num] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_num] +
                                                e_mudec_2d.integrand[id_num]) +
                              muon_beta_j_em_integrals.integrand[id_num]);

        m1_opacities_non_th_separated.eta_non_th[id_num] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_num] +
                                  e_nms_2d.integrand[id_num]);

        if (n[id_num] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_num] =
                kBS_FourPi_hc3 / (c_light * n[id_num]) *
                (kBS_FourPi_hc3 * (n_integrals_2d.integrand[id_num_kappa] +
                                   n_mudec_2d.integrand[id_num_kappa]) +
                 muon_beta_n_abs_integrals.integrand[id_num]);
        }
        if (J[id_num] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_num] =
                kBS_FourPi_hc3 / (c_light * J[id_num]) *
                (kBS_FourPi_hc3 * (e_integrals_2d.integrand[id_num_kappa] +
                                   e_mudec_2d.integrand[id_num_kappa]) +
                 muon_beta_j_abs_integrals.integrand[id_num]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_num] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_num]) *
                (e_neps_2d.integrand[id_num_kappa] + e_nms_2d.integrand[id_num_kappa]);

            m1_opacities_non_th_separated.kappa_s[id_num] = kBS_FourPi_hc3 / (c_light * J[id_num]) *
                                        iso_integrals.integrand[id_num];
        }

        /* Muon anti-neutrinos */
        m1_opacities_non_th_separated.eta_0[id_anum] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_anum] +
                              muon_beta_n_em_integrals.integrand[id_anum]);

        m1_opacities_non_th_separated.eta_th[id_anum] =
            kBS_FourPi_hc3 * (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_anum] +
                              muon_beta_j_em_integrals.integrand[id_anum]);

        m1_opacities_non_th_separated.eta_non_th[id_anum] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_anum] +
                                  e_nms_2d.integrand[id_anum]);

        if (n[id_anum] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_anum] =
                n[id_anum] == zero ?
                    zero :
                    kBS_FourPi_hc3 / (c_light * n[id_anum]) *
                        (kBS_FourPi_hc3 * n_integrals_2d.integrand[id_anum_kappa] +
                         muon_beta_n_abs_integrals.integrand[id_anum]);
        }
        if (J[id_anum] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_anum] =
                kBS_FourPi_hc3 / (c_light * J[id_anum]) *
                (kBS_FourPi_hc3 * e_integrals_2d.integrand[id_anum_kappa] +
                 muon_beta_j_abs_integrals.integrand[id_anum]);

            m1_opacities_non_th_separated.kappa_a_non_th[id_anum] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anum]) *
                (e_neps_2d.integrand[id_anum_kappa] + e_nms_2d.integrand[id_anum_kappa]);

            m1_opacities_non_th_separated.kappa_s[id_anum] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anum]) *
                                            iso_integrals.integrand[id_anum];
        }

        /* Tau neutrinos */
        m1_opacities_non_th_separated.eta_0[id_nut] =
            kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[id_nut];

        m1_opacities_non_th_separated.eta_th[id_nut] =
            kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[id_nut];

        m1_opacities_non_th_separated.eta_non_th[id_nut] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_nut] +
                                  e_nms_2d.integrand[id_nut]);

        if (n[id_nut] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_nut] =
                kBS_FourPi_hc3_sqr / (c_light * n[id_nut]) *
                    n_integrals_2d.integrand[id_nut_kappa];
        }
        if (J[id_nut] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_nut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nut]) *
                    e_integrals_2d.integrand[id_nut_kappa];

            m1_opacities_non_th_separated.kappa_a_non_th[id_nut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_nut]) *
                (e_neps_2d.integrand[id_nut_kappa] + e_nms_2d.integrand[id_nut_kappa]);

            m1_opacities_non_th_separated.kappa_s[id_nut] = kBS_FourPi_hc3 / (c_light * J[id_nut]) *
                                        iso_integrals.integrand[id_nut];
        }

        /* Tau anti-neutrinos */
        m1_opacities_non_th_separated.eta_0[id_anut] =
            kBS_FourPi_hc3_sqr * n_integrals_2d.integrand[id_anut];

        m1_opacities_non_th_separated.eta_th[id_anut] =
            kBS_FourPi_hc3_sqr * e_integrals_2d.integrand[id_anut];

        m1_opacities_non_th_separated.eta_non_th[id_anut] =
            kBS_FourPi_hc3_sqr * (e_neps_2d.integrand[id_anut] +
                                  e_nms_2d.integrand[id_anut]);

        if (n[id_anut] > THRESHOLD_N)
        {
            m1_opacities_non_th_separated.kappa_0_a[id_anut] =
                n[id_anut] == zero ?
                    zero :
                    kBS_FourPi_hc3_sqr / (c_light * n[id_anut]) *
                        n_integrals_2d.integrand[id_anut_kappa];
        }
        if (J[id_anut] > THRESHOLD_J)
        {
            m1_opacities_non_th_separated.kappa_a_th[id_anut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anut]) *
                    e_integrals_2d.integrand[id_anut_kappa];

            m1_opacities_non_th_separated.kappa_a_non_th[id_anut] =
                kBS_FourPi_hc3_sqr / (c_light * J[id_anut]) *
                (e_neps_2d.integrand[id_anut_kappa] + e_nms_2d.integrand[id_anut_kappa]);

            m1_opacities_non_th_separated.kappa_s[id_anut] = kBS_FourPi_hc3 /
                                            (c_light * J[id_anut]) *
                                            iso_integrals.integrand[id_anut];
        }

    }

    else  //fallback
    {
        BS_ASSERT(false, "Wrong number of nu species: %d", total_num_species);
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
        0}; //{.em_e = 0., .abs_e = 0., .em_x = 0., .abs_x = 0.}
            //{.em_e = 0., .abs_e = 0., .em_m = 0., .abs_m = 0., .em_t = 0., .abs_t = 0.}
    if (opacity_flags.use_pair)
    {
        my_grey_opacity_params->kernel_pars.pair_kernel_params.omega_prime =
            nu_bar;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.cos_theta = one;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.filter    = zero;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.lmax      = zero;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu        = one;
        my_grey_opacity_params->kernel_pars.pair_kernel_params.mu_prime  = one;

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
    MyKernelOutput inelastic_NEPS_kernels_m1 = {0};
    if (opacity_flags.use_inelastic_NEPS)
    {
        my_grey_opacity_params->kernel_pars.inelastic_kernel_params
            .omega_prime     = nu_bar;
        inelastic_NEPS_kernels_m1 = InelasticNEPSKernels(
            &my_grey_opacity_params->kernel_pars.inelastic_kernel_params,
            &my_grey_opacity_params->eos_pars);
    }

    // compute the inelastic NMS kernels
    MyKernelOutput inelastic_NMS_kernels_m1 = {0};
    if (opacity_flags.use_inelastic_NMS)
    {
        my_grey_opacity_params->kernel_pars.inelastic_kernel_params
            .omega_prime     = nu_bar;
        if (opacity_pars.NMS_implementation == NMS_KernelInterp)
        {
            inelastic_NMS_kernels_m1 = InelasticNMSKernels_DirectInterp(
                &my_grey_opacity_params->kernel_pars.inelastic_kernel_params,
                &my_grey_opacity_params->eos_pars);
        }
        else if (opacity_pars.NMS_implementation == NMS_SemiAnalytical)
        {
            inelastic_NMS_kernels_m1 = InelasticNMSKernels_SemiAnalytical(
                &my_grey_opacity_params->kernel_pars.inelastic_kernel_params,
                &my_grey_opacity_params->eos_pars);
        }
        else
        {
            BS_ASSERT(false, "Unknown NMS implementation: %d",
                      (int)opacity_pars.NMS_implementation);
        }
    }

    // Compute the muon decay kernels
    MyKernelOutput muon_decay_kernels_m1 = {0};
    if (opacity_flags.use_muon_decay)
    {
        my_grey_opacity_params->kernel_pars.muon_decay_kernel_params
            .omega_anue = nu_bar;
        muon_decay_kernels_m1 = MuonDecayKernels(
            &my_grey_opacity_params->kernel_pars.muon_decay_kernel_params,
            &my_grey_opacity_params->eos_pars);
    }

    // Block factor
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

    // Total production term
    BS_REAL pro_term[total_num_species] = {0};

    pro_term[id_nue] =
        (pair_kernels_m1.em[id_nue] + brem_kernels_m1.em[id_nue]) *
        block_factor[id_anue];
    pro_term[id_anue] =
        (pair_kernels_m1.em[id_anue] + brem_kernels_m1.em[id_anue]) *
        block_factor[id_nue];

    if constexpr (total_num_species == 4){

        pro_term[id_nux] =
            (pair_kernels_m1.em[id_nux] + brem_kernels_m1.em[id_nux]) *
            block_factor[id_anux];
        pro_term[id_anux] =
            (pair_kernels_m1.em[id_anux] + brem_kernels_m1.em[id_anux]) *
            block_factor[id_nux];

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            pro_term[idx] += inelastic_NEPS_kernels_m1.em[idx] * g_nu_bar[idx];
        }
    }
    else if constexpr (total_num_species == 6){

        pro_term[id_num] =
            (pair_kernels_m1.em[id_num] + brem_kernels_m1.em[id_num]) *
            block_factor[id_anum];
        pro_term[id_anum] =
            (pair_kernels_m1.em[id_anum] + brem_kernels_m1.em[id_anum]) *
            block_factor[id_num];
        pro_term[id_nut] =
            (pair_kernels_m1.em[id_nut] + brem_kernels_m1.em[id_nut]) *
            block_factor[id_anut];
        pro_term[id_anut] =
            (pair_kernels_m1.em[id_anut] + brem_kernels_m1.em[id_anut]) *
            block_factor[id_nut];

        //muon decay contribution
        pro_term[id_num] += muon_decay_kernels_m1.em[id_num] * 
                            block_factor[id_anue];
        pro_term[id_anue] += muon_decay_kernels_m1.em[id_anue] * 
                            block_factor[id_num];

        //inelastic scattering contribution
        for (int idx = 0; idx < total_num_species; ++idx)
        {
            pro_term[idx] += (inelastic_NEPS_kernels_m1.em[idx]
                            + inelastic_NMS_kernels_m1.em[idx]) * g_nu_bar[idx];
        }
    }

    // Total annihilation term
    BS_REAL ann_term[total_num_species] = {0};

    ann_term[id_nue] =
        (pair_kernels_m1.abs[id_nue] + brem_kernels_m1.abs[id_nue]) *
        g_nu_bar[id_anue];
    ann_term[id_anue] =
        (pair_kernels_m1.abs[id_anue] + brem_kernels_m1.abs[id_anue]) *
        g_nu_bar[id_nue];

    if constexpr (total_num_species == 4)
    {
        ann_term[id_nux] =
            (pair_kernels_m1.abs[id_nux] + brem_kernels_m1.abs[id_nux]) *
            g_nu_bar[id_anux];
        ann_term[id_anux] =
            (pair_kernels_m1.abs[id_anux] + brem_kernels_m1.abs[id_anux]) *
            g_nu_bar[id_nux];

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            ann_term[idx] += inelastic_NEPS_kernels_m1.abs[idx] * block_factor[idx];
        }
    }

    else if constexpr (total_num_species == 6)
    {
        ann_term[id_num] =
            (pair_kernels_m1.abs[id_num] + brem_kernels_m1.abs[id_num]) *
            g_nu_bar[id_anum];
        ann_term[id_anum] =
            (pair_kernels_m1.abs[id_anum] + brem_kernels_m1.abs[id_anum]) *
            g_nu_bar[id_num];
        ann_term[id_nut] =
            (pair_kernels_m1.abs[id_nut] + brem_kernels_m1.abs[id_nut]) *
            g_nu_bar[id_anut];
        ann_term[id_anut] =
            (pair_kernels_m1.abs[id_anut] + brem_kernels_m1.abs[id_anut]) *
            g_nu_bar[id_nut];

        //muon decay contribution
        ann_term[id_num] += muon_decay_kernels_m1.abs[id_num] * 
                            g_nu_bar[id_anue];
        ann_term[id_anue] += muon_decay_kernels_m1.abs[id_anue] * 
                            g_nu_bar[id_num];

        //inelastic scattering contribution
        for (int idx = 0; idx < total_num_species; ++idx)
        {
            ann_term[idx] += (inelastic_NEPS_kernels_m1.abs[idx]
                            + inelastic_NMS_kernels_m1.abs[idx]) * block_factor[idx];
        }
    }
    //////////////////////////////////////////////
    ////// ONLY FOR COMPARISON WITH NULIB ////////
    //////////////////////////////////////////////
    if (kirchoff_flag)
    {
        ann_term[id_nue] +=
            pair_kernels_m1.em[id_nue] * g_nu_bar[id_anue] / g_nu[id_nue];
        ann_term[id_anue] +=
            pair_kernels_m1.em[id_anue] * g_nu_bar[id_nue] / g_nu[id_anue];

        if constexpr (total_num_species == 4)
        {
            ann_term[id_nux] +=
                pair_kernels_m1.em[id_nux] * g_nu_bar[id_anux] / g_nu[id_nux];
            ann_term[id_anux] +=
                pair_kernels_m1.em[id_anux] * g_nu_bar[id_nux] / g_nu[id_anux];
        }
        else if constexpr (total_num_species == 6)
        {
            ann_term[id_num] +=
                pair_kernels_m1.em[id_num] * g_nu_bar[id_anum] / g_nu[id_num];
            ann_term[id_anum] +=
                pair_kernels_m1.em[id_anum] * g_nu_bar[id_num] / g_nu[id_anum];
                ann_term[id_nut] +=
                pair_kernels_m1.em[id_nut] * g_nu_bar[id_anut] / g_nu[id_nut];
            ann_term[id_anut] +=
                pair_kernels_m1.em[id_anut] * g_nu_bar[id_nut] / g_nu[id_anut];
        }
    }
    //////////////////////////////////////////////

    BS_REAL integrand_1[total_num_species], integrand_2[total_num_species];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        integrand_1[idx] = POW2(nu_bar) * pro_term[idx];
        integrand_2[idx] = POW2(nu_bar) * ann_term[idx];
    }

    MyQuadratureIntegrand result = {.n = double_total_num_species};

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        result.integrand[idx] = integrand_1[idx];
        result.integrand[total_num_species + idx] = integrand_2[idx];
    }

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
    constexpr BS_REAL two_over_three = 2. / 3.;
    constexpr BS_REAL four_pi = 4 * kBS_Pi;
    constexpr BS_REAL c_light = kBS_Clight;
    constexpr BS_REAL Mmu_over_three = kBS_Mmu / 3.;
    BS_REAL wmin, wmax;

    constexpr int id_nue_abs = total_num_species + id_nue;
    constexpr int id_anue_abs = total_num_species + id_anue;
    constexpr int id_nux_abs = total_num_species + id_nux;
    constexpr int id_anux_abs = total_num_species + id_anux;
    constexpr int id_num_abs = total_num_species + id_num;
    constexpr int id_anum_abs = total_num_species + id_anum;
    constexpr int id_nut_abs = total_num_species + id_nut;
    constexpr int id_anut_abs = total_num_species + id_anut;

    // Extremals for NMS integration
    if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_KernelInterp)
    {
        wmin = NMS_w_min;
        wmax = NMS_w_max;
    }
    else if (my_grey_opacity_params->opacity_pars.NMS_implementation == NMS_SemiAnalytical)
    {
        wmin = NMSParams_w_min;
        wmax = NMSParams_w_max;
    }

    my_grey_opacity_params->kernel_pars.pair_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.brem_kernel_params.omega      = nu;
    my_grey_opacity_params->kernel_pars.inelastic_kernel_params.omega = nu;
    my_grey_opacity_params->kernel_pars.muon_decay_kernel_params.omega_numu = nu;

    GreyOpacityParams local_grey_params = *my_grey_opacity_params;
    local_grey_params.opacity_flags.use_inelastic_NEPS = 0;
    local_grey_params.opacity_flags.use_inelastic_NMS = 0;
    local_grey_params.opacity_flags.use_muon_decay = 0;

    // set up 1d integration
    MyFunctionMultiD integrand_m1_1d;
    MyQuadratureIntegrand integrand_m1_1d_info = {.n = double_total_num_species};
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

    constexpr BS_REAL temp_multiple = 0.5 * 4.364;

    // s can vary in terms of nu flavors and type of integrand (em or abs)
    BS_REAL s_pair[double_total_num_species], s_neps[double_total_num_species];
    BS_REAL s_mudec[double_total_num_species], s_nms[double_total_num_species];

    for (int i = 0; i < double_total_num_species; ++i)
    {
        s_pair[i] = temp_multiple * my_grey_opacity_params->eos_pars.temp;
        s_neps[i] = nu;
        s_nms[i] = std::max(3. * wmin, std::min(nu, two_over_three * wmax));
        s_mudec[i] = Mmu_over_three;
    }

    MyQuadratureIntegrand integrals_pair_1d =
        GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_pair);

    MyQuadratureIntegrand integrals_neps_1d = {0};
    if (my_grey_opacity_params->opacity_flags.use_inelastic_NEPS == 1)
    {
        local_grey_params.opacity_flags                     = {0};
        local_grey_params.opacity_flags.use_inelastic_NEPS = 1;
        integrand_m1_1d.params = &local_grey_params;
        integrals_neps_1d =
            GaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_neps);
    }

    MyQuadratureIntegrand integrals_nms_1d = {0};
    if (my_grey_opacity_params->opacity_flags.use_inelastic_NMS == 1)
    {
        local_grey_params.opacity_flags                     = {0};
        local_grey_params.opacity_flags.use_inelastic_NMS = 1;
        integrand_m1_1d.params = &local_grey_params;
        integrals_nms_1d =
            MuonReactionsGaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_nms, wmin, wmax);
    }

    MyQuadratureIntegrand integrals_muon_decay_1d = {0};
    if (my_grey_opacity_params->opacity_flags.use_muon_decay == 1)
    {
        local_grey_params.opacity_flags                     = {0};
        local_grey_params.opacity_flags.use_muon_decay = 1;
        integrand_m1_1d.params = &local_grey_params;
        integrals_muon_decay_1d =
            MuonReactionsGaussLegendreIntegrate1D(quad_1d, &integrand_m1_1d, s_mudec, 
                                                    MuonDecay_wanue_min, MuonDecay_wanue_max);
    }

    MyOpacity abs_em_beta = {0};
    if (my_grey_opacity_params->opacity_flags.use_abs_em)
    {
        abs_em_beta = ElAbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                                 &my_grey_opacity_params->eos_pars); // [s^-1]
    }

    MyOpacity abs_em_muonic_beta = {0};
    if (my_grey_opacity_params->opacity_flags.use_muonic_beta)
    {
        abs_em_muonic_beta = MuonAbsOpacity(nu, &my_grey_opacity_params->opacity_pars,
                                 &my_grey_opacity_params->eos_pars); // [s^-1]
    }

    BS_REAL iso_scatt = zero;
    if (my_grey_opacity_params->opacity_flags.use_iso)
    {
        iso_scatt = IsoScattLegCoeff(nu, &my_grey_opacity_params->opacity_pars,
                                     &my_grey_opacity_params->eos_pars, 0);
    }

    SpectralOpacities sp_opacities;

    if constexpr (total_num_species == 4)
    {
        // Emissivities
        sp_opacities.j[id_nue]  = abs_em_beta.em[id_nue] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nue] +
                                                    integrals_neps_1d.integrand[id_nue]);
        sp_opacities.j[id_anue] = abs_em_beta.em[id_anue] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anue] +
                                                    integrals_neps_1d.integrand[id_anue]);
        sp_opacities.j[id_nux]  = kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nux] +
                                                    integrals_neps_1d.integrand[id_nux]);
        sp_opacities.j[id_anux] = kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anux] +
                                                    integrals_neps_1d.integrand[id_anux]);

        // Absorsivities
        sp_opacities.kappa[id_nue] =
            (abs_em_beta.abs[id_nue] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nue_abs] +
                               integrals_neps_1d.integrand[id_nue_abs])) /
            c_light;
        sp_opacities.kappa[id_anue] =
            (abs_em_beta.abs[id_anue] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anue_abs] +
                               integrals_neps_1d.integrand[id_anue_abs])) /
            c_light;
        sp_opacities.kappa[id_nux] =
            kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nux_abs] +
                              integrals_neps_1d.integrand[id_nux_abs]) /
            c_light;
        sp_opacities.kappa[id_anux] =
            kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anux_abs] +
                              integrals_neps_1d.integrand[id_anux_abs]) /
            c_light;

        
        // Scattering quantities
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
    }
    
    else if constexpr (total_num_species == 6)
    {
        // Emissivities
        sp_opacities.j[id_nue]  = abs_em_beta.em[id_nue] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nue] +
                                                    integrals_neps_1d.integrand[id_nue] +
                                                    integrals_nms_1d.integrand[id_nue]);
        sp_opacities.j[id_anue] = abs_em_beta.em[id_anue] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anue] +
                                                    integrals_neps_1d.integrand[id_anue] +
                                                    integrals_nms_1d.integrand[id_anue]  +
                                                    integrals_muon_decay_1d.integrand[id_anue]);
        sp_opacities.j[id_num]  = abs_em_muonic_beta.em[id_num] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_num] +
                                                    integrals_neps_1d.integrand[id_num] +
                                                    integrals_nms_1d.integrand[id_num]  +
                                                    integrals_muon_decay_1d.integrand[id_num]);
        sp_opacities.j[id_anum] = abs_em_muonic_beta.em[id_anum] +
                                  kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anum] +
                                                    integrals_neps_1d.integrand[id_anum] +
                                                    integrals_nms_1d.integrand[id_anum]);
        sp_opacities.j[id_nut]  = kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nut] +
                                                    integrals_neps_1d.integrand[id_nut] +
                                                    integrals_nms_1d.integrand[id_nut]);
        sp_opacities.j[id_anut] = kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anut] +
                                                    integrals_neps_1d.integrand[id_anut] +
                                                    integrals_nms_1d.integrand[id_anut]);

        // Absorsivities
        sp_opacities.kappa[id_nue] =
            (abs_em_beta.abs[id_nue] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nue_abs] +
                               integrals_neps_1d.integrand[id_nue_abs] +
                               integrals_nms_1d.integrand[id_nue_abs])) /
            c_light;
        sp_opacities.kappa[id_anue] =
            (abs_em_beta.abs[id_anue] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anue_abs] +
                               integrals_neps_1d.integrand[id_anue_abs] +
                               integrals_nms_1d.integrand[id_anue_abs]  +
                               integrals_muon_decay_1d.integrand[id_anue_abs])) /
            c_light;
        sp_opacities.kappa[id_num] =
            (abs_em_muonic_beta.abs[id_num] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_num_abs] +
                               integrals_neps_1d.integrand[id_num_abs] +
                               integrals_nms_1d.integrand[id_num_abs]  +
                               integrals_muon_decay_1d.integrand[id_num_abs])) /
            c_light;
        sp_opacities.kappa[id_anum] =
            (abs_em_muonic_beta.abs[id_anum] +
             kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anum_abs] +
                               integrals_neps_1d.integrand[id_anum_abs] +
                               integrals_nms_1d.integrand[id_anum_abs])) /
            c_light;
        sp_opacities.kappa[id_nut] =
            kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_nut_abs] +
                               integrals_neps_1d.integrand[id_nut_abs] +
                               integrals_nms_1d.integrand[id_nut_abs]) /
            c_light;
        sp_opacities.kappa[id_anut] =
            kBS_FourPi_hc3 * (integrals_pair_1d.integrand[id_anut_abs] +
                               integrals_neps_1d.integrand[id_anut_abs] +
                               integrals_nms_1d.integrand[id_anut_abs]) /
            c_light;

        // Scattering quantities
        sp_opacities.j_s[id_nue]  = four_pi * POW2(nu) * g_nu[id_nue] * iso_scatt;
        sp_opacities.j_s[id_anue] = four_pi * POW2(nu) * g_nu[id_anue] * iso_scatt;
        sp_opacities.j_s[id_num]  = four_pi * POW2(nu) * g_nu[id_num] * iso_scatt;
        sp_opacities.j_s[id_anum] = four_pi * POW2(nu) * g_nu[id_anum] * iso_scatt;
        sp_opacities.j_s[id_nut]  = four_pi * POW2(nu) * g_nu[id_nut] * iso_scatt;
        sp_opacities.j_s[id_anut] = four_pi * POW2(nu) * g_nu[id_anut] * iso_scatt;

        sp_opacities.kappa_s[id_nue] =
            four_pi * POW2(nu) * (one - g_nu[id_nue]) * iso_scatt / c_light;
        sp_opacities.kappa_s[id_anue] =
            four_pi * POW2(nu) * (one - g_nu[id_anue]) * iso_scatt / c_light;
        sp_opacities.kappa_s[id_num] =
            four_pi * POW2(nu) * (one - g_nu[id_num]) * iso_scatt / c_light;
        sp_opacities.kappa_s[id_anum] =
            four_pi * POW2(nu) * (one - g_nu[id_anum]) * iso_scatt / c_light;
        sp_opacities.kappa_s[id_nut] =
            four_pi * POW2(nu) * (one - g_nu[id_nut]) * iso_scatt / c_light;
        sp_opacities.kappa_s[id_anut] =
            four_pi * POW2(nu) * (one - g_nu[id_anut]) * iso_scatt / c_light;
    }
    
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

#endif // BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_HPP_
