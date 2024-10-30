#ifndef BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_
#define BNS_NURATES_SRC_OPACITIES_M1_OPACITIES_H_

//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  m1_opacities.hpp
//  \brief header file for all integration routines

#include "bns_nurates.hpp"
#include "opacities.hpp"
#include "kernels.hpp"
#include "integration.hpp"
#include "distribution.hpp"

/* Thresholds on the neutrino energy number and energy density. If values are
below the thresholds, absorption opacities or scattering opacities are set to 0.
*/
static const double THRESHOLD_N = 1e0;
static const double THRESHOLD_J = 1e-4;

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
                             GreyOpacityParams* grey_pars, const double* t,
                             double out[][n_max])
{
    static const double four_pi = 4. * kPi;

    const int n = quad->nx;

    double nu, iso_scatt, aux;

    double g_nu[total_num_species];

    for (int i = 0; i < n; i ++) {
        for (int idx = 0; idx < total_num_species; idx ++) {
            nu = t[idx] * quad->points[i];

            // compute the neutrino distribution function
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            iso_scatt = IsoScattTotal(nu, &grey_pars->opacity_pars,
                          &grey_pars->eos_pars);

            aux = four_pi * nu * nu * nu * nu * nu * iso_scatt;

            out[idx][i] = g_nu[idx] * aux;
       
            nu =  t[idx] / quad->points[i];

            // compute the neutrino distribution function
            g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            iso_scatt = IsoScattTotal(nu, &grey_pars->opacity_pars,
                          &grey_pars->eos_pars);

            aux = four_pi * nu * nu * nu * nu * nu * iso_scatt;

            out[idx][n + i] = g_nu[idx] * aux;
        }
    }   
    
    return;
}

KOKKOS_INLINE_FUNCTION
void Beta1DIntegrand(const MyQuadrature* quad,
                    GreyOpacityParams* grey_pars, const double *t,
                    double out_em[][n_max], double out_ab[][n_max], const int stim_abs) {
        const int n = quad->nx;
        double nu, nu_sqr, g_nu;
        MyOpacity abs_em_beta;
        
        if (stim_abs == 1) {
            for (int i = 0; i < n; i ++) {
                nu = t[id_nue] * quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

                abs_em_beta =
                    StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                
                out_em[id_nue][i] = nu_sqr * abs_em_beta.em[id_nue];
                out_ab[id_nue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)

                nu = t[id_anue] * quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

                abs_em_beta =
                    StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_anue][i] = nu_sqr * abs_em_beta.em[id_anue];
                out_ab[id_anue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue]; // ab = em + ab (stimulated absorption)

                nu = t[id_nue] / quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_nue);
                
                abs_em_beta =
                    StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_nue][n + i] = nu_sqr * abs_em_beta.em[id_nue];
                out_ab[id_nue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)

                nu = t[id_anue] / quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_anue);
               
                abs_em_beta =
                    StimAbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_anue][n + i] = nu_sqr * abs_em_beta.em[id_anue];
                out_ab[id_anue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue]; // ab = em + ab (stimulated absorption)
            }       
        } else {
            for (int i = 0; i < n; i ++) {
                nu = t[id_nue] * quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_nue);

                abs_em_beta =
                    AbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                
                out_em[id_nue][i] = nu_sqr * abs_em_beta.em[id_nue];
                out_ab[id_nue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)

                nu = t[id_anue] * quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_anue);

                abs_em_beta =
                    AbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_anue][i] = nu_sqr * abs_em_beta.em[id_anue];
                out_ab[id_anue][i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue]; // ab = em + ab (stimulated absorption)

                nu = t[id_nue] / quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_nue);
                
                abs_em_beta =
                    AbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_nue][n + i] = nu_sqr * abs_em_beta.em[id_nue];
                out_ab[id_nue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)

                nu = t[id_anue] / quad->points[i];
                nu_sqr = nu * nu;
                g_nu = TotalNuF(nu, &grey_pars->distr_pars, id_anue);
               
                abs_em_beta =
                    AbsOpacity(nu, &grey_pars->opacity_pars,
                                   &grey_pars->eos_pars); // [s^-1]

                out_em[id_anue][n + i] = nu_sqr * abs_em_beta.em[id_anue];
                out_ab[id_anue][n + i] = nu_sqr * g_nu * abs_em_beta.abs[id_anue]; // ab = em + ab (stimulated absorption)
            }       
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void AddBetaReactionToIntegrand(int n, double* nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D *out, const int stim_abs) {
        double nu;
        MyOpacity abs_em_beta;
        
        if (stim_abs == 1) {
        for (int i = 0; i < 2 * n; i ++) {
            nu = nu_array[i];
  
            abs_em_beta =
                StimAbsOpacity(nu, &grey_pars->opacity_pars,
                               &grey_pars->eos_pars); // [s^-1]
           
            for (int j = 0; i < 2 * n; j ++) {
                out->m1_mat_em[id_nue][i][j] += abs_em_beta.em[id_nue];
                out->m1_mat_em[id_anue][i][j] += abs_em_beta.em[id_anue];

                out->m1_mat_ab[id_nue][i][j] += abs_em_beta.abs[id_nue]; // ab = em + ab (stimulated absorption)
                out->m1_mat_ab[id_anue][i][j] += abs_em_beta.abs[id_anue];
            }       
        }
    } else {
        for (int i = 0; i < 2 * n; i ++) {
            nu = nu_array[i];
  
            abs_em_beta =
                AbsOpacity(nu, &grey_pars->opacity_pars,
                               &grey_pars->eos_pars); // [s^-1]
           
            for (int j = 0; i < 2 * n; j ++) {
                out->m1_mat_em[id_nue][i][j] += abs_em_beta.em[id_nue];
                out->m1_mat_em[id_anue][i][j] += abs_em_beta.em[id_anue];

                out->m1_mat_ab[id_nue][i][j] += abs_em_beta.abs[id_nue]; // ab = ab (not stimulated absorption)
                out->m1_mat_ab[id_anue][i][j] += abs_em_beta.abs[id_anue];
            }       
        }
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void AddPairKernelsToIntegrand(int n, double* nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D* out) {
      MyKernelOutput pair_1, pair_2;

      grey_pars->kernel_pars.pair_kernel_params.cos_theta = 1.;
      grey_pars->kernel_pars.pair_kernel_params.filter    = 0.;
      grey_pars->kernel_pars.pair_kernel_params.lmax      = 0;
      grey_pars->kernel_pars.pair_kernel_params.mu        = 1.;
      grey_pars->kernel_pars.pair_kernel_params.mu_prime  = 1.;

      for (int i = 0; i < 2 * n; i++)
      {
        grey_pars->kernel_pars.pair_kernel_params.omega = nu_array[i];
        grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[i];

        PairKernelsM1Test(&grey_pars->eos_pars,
                          &grey_pars->kernel_pars.pair_kernel_params,
                          &pair_1, &pair_2);
        
        for (int idx = 0; idx < total_num_species; idx++)
            {
                out->m1_mat_em[idx][i][i] += pair_1.em[idx];
                out->m1_mat_ab[idx][i][i] += pair_1.abs[idx];
            }

        for (int j = i + 1; j < 2 * n; j++)
        {
            grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[j];

            PairKernelsM1Test(&grey_pars->eos_pars,
                              &grey_pars->kernel_pars.pair_kernel_params,
                              &pair_1, &pair_2);

            for (int idx = 0; idx < total_num_species; idx++)
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
void AddBremKernelsToIntegrand(int n, double* nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D* out) {
    MyKernelOutput brem_ker;

        if (grey_pars->opacity_pars.use_BRT_brem == true)
        {
          for (int i = 0; i < 2 * n; i++)
          {

            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[i];
            
            // compute the brem kernels
            brem_ker = BremKernelsBRT06(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
            
            for (int idx = 0; idx < total_num_species; idx++)
            {
              out->m1_mat_em[idx][i][i] += brem_ker.em[0];
              out->m1_mat_ab[idx][i][i] += brem_ker.abs[0];
            }
            
            for (int j = i + 1; j < 2 * n; j++)
            {
              grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[j];
            
              // compute the brem kernels
              brem_ker = BremKernelsBRT06(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
            
              for (int idx = 0; idx < total_num_species; idx++)
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

          for (int i = 0; i < 2 * n; i++)
          {

            grey_pars->kernel_pars.brem_kernel_params.omega       = nu_array[i];
            grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[i];
            
            // compute the brem kernels
            brem_ker = BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
            
            for (int idx = 0; idx < total_num_species; idx++)
            {
              out->m1_mat_em[idx][i][i] += brem_ker.em[0];
              out->m1_mat_ab[idx][i][i] += brem_ker.abs[0];
            }

            for (int j = i + 1; j < 2 * n; j++)
            {
              // compute the brem kernels
              grey_pars->kernel_pars.brem_kernel_params.omega_prime = nu_array[j];
              brem_ker = BremKernelsLegCoeff(&grey_pars->kernel_pars.brem_kernel_params, &grey_pars->eos_pars);
            
              for (int idx = 0; idx < total_num_species; idx++)
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
void AddInelKernelsToIntegrand(int n, double* nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D* out) {
        double nu, nu_bar;
        double g_nu[total_num_species], g_nu_bar[total_num_species];
        double block_factor_nu[total_num_species], block_factor_nu_bar[total_num_species];
       MyKernelOutput inel_1, inel_2;

       for (int i = 0; i < 2 * n; i++)
       {
         nu = nu_array[i];

         for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                
                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx]     = 1. - g_nu[idx];
                }
                else
                {
                    block_factor_nu[idx]     = 1.;
                }
            }

           // compute the pair kernels
           grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
           grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;
           
           inel_1 = InelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars);
           
           for (int idx = 0; idx < total_num_species; idx++)
            {
                out->m1_mat_em[idx][i][i] += inel_1.em[idx] * g_nu[idx];
                out->m1_mat_ab[idx][i][i] += inel_1.abs[idx] * block_factor_nu[idx];
            }
          

         for (int j = i + 1; j < 2 * n; j++)
         {
           nu_bar = nu_array[j];
           
           for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);
                
                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu_bar[idx] = 1. - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu_bar[idx] = 1.;
                }
            }
           
           // compute the pair kernels
           grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu;
           grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu_bar;
           
           inel_1 = InelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars);
           
           grey_pars->kernel_pars.inelastic_kernel_params.omega       = nu_bar;
           grey_pars->kernel_pars.inelastic_kernel_params.omega_prime = nu;
   
           inel_2 = InelasticScattKernels(&grey_pars->kernel_pars.inelastic_kernel_params, &grey_pars->eos_pars);
           
           for (int idx = 0; idx < total_num_species; idx++)
            {
                out->m1_mat_em[idx][i][j] += inel_1.em[idx] * g_nu_bar[idx];
                out->m1_mat_em[idx][j][i] += inel_2.em[idx] * g_nu[idx];
                
                out->m1_mat_ab[idx][i][j] += inel_1.abs[idx] * block_factor_nu_bar[idx];
                out->m1_mat_ab[idx][j][i] += inel_2.abs[idx] * block_factor_nu[idx];
            }
          }
       }
       return;
    } 

KOKKOS_INLINE_FUNCTION
void WeightNuNuBarReactionsWithDistr(int n, double *nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D *out) {
    double nu, nu_bar;
    double g_nu[total_num_species], g_nu_bar[total_num_species];
    double block_factor_nu[total_num_species], block_factor_nu_bar[total_num_species];

   for (int i = 0; i < 2 * n; i++)
    {
        nu = nu_array[i];

        for (int idx = 0; idx < total_num_species; idx++)
        {
          g_nu[idx] = TotalNuF(nu, &grey_pars->distr_pars, idx);

            if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu[idx] = 1. - g_nu[idx];
                }
                else
                {
                    block_factor_nu[idx] = 1.;
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

        for (int j = i + 1; j < 2 * n; j++)
        {
            nu_bar = nu_array[j];

            for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);

                if (grey_pars->opacity_pars.neglect_blocking == false)
                {
                    block_factor_nu_bar[idx] = 1. - g_nu_bar[idx];
                }
                else
                {
                    block_factor_nu_bar[idx] = 1.;
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
void AddCommonWeightsToIntegrand(int n, double* nu_array, GreyOpacityParams *grey_pars, M1MatrixKokkos2D *out, int stim_abs) {
    double nu, nu_bar, nu_squared, nu_fourth;
    double g_nu[total_num_species], g_nu_bar[total_num_species];

    assert((stim_abs == 0) || (stim_abs == 1));

    if (stim_abs == 1) {
     for (int i = 0; i < 2 * n; i++)
       {
         nu = nu_array[i];
         nu_squared = nu * nu;
         nu_fourth = nu_squared * nu_squared;
         
         for (int idx = 0; idx < total_num_species; idx++)
         {
            g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                
            out->m1_mat_ab[idx][i][i] = nu_fourth * g_nu[idx] *
                    (out->m1_mat_em[idx][i][i] + out->m1_mat_ab[idx][i][i]);
            out->m1_mat_em[idx][i][i] *= nu_fourth;
          }

         for (int j = i + 1; j < 2 * n; j++)
         {
            nu_bar = nu_array[j];
            nu_fourth = nu_squared * nu_bar * nu_bar;

            for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);
                
                out->m1_mat_ab[idx][i][j] = nu_fourth * g_nu[idx] *
                    (out->m1_mat_em[idx][i][j] + out->m1_mat_ab[idx][i][j]);
                out->m1_mat_ab[idx][j][i] =
                    nu_fourth * g_nu_bar[idx] *
                    (out->m1_mat_em[idx][j][i] + out->m1_mat_ab[idx][j][i]);
                
                out->m1_mat_em[idx][i][j] *= nu_fourth;
                out->m1_mat_em[idx][j][i] *= nu_fourth;

            }
        }
    }
    } else {
    for (int i = 0; i < 2 * n; i++)
       {
         nu = nu_array[i];
         nu_squared = nu * nu;
         nu_fourth = nu_squared * nu_squared;
         
         for (int idx = 0; idx < total_num_species; idx++)
         {
            g_nu[idx]     = TotalNuF(nu, &grey_pars->distr_pars, idx);
                
            out->m1_mat_ab[idx][i][i] *= nu_fourth * g_nu[idx];
            out->m1_mat_em[idx][i][i] *= nu_fourth * (1. - g_nu[idx]);
          }

         for (int j = i + 1; j < 2 * n; j++)
         {
            nu_bar = nu_array[j];
            nu_fourth = nu_squared * nu_bar * nu_bar;

            for (int idx = 0; idx < total_num_species; idx++)
            {
                g_nu_bar[idx] = TotalNuF(nu_bar, &grey_pars->distr_pars, idx);
                
                out->m1_mat_ab[idx][i][j] *= nu_fourth * g_nu[idx];
                out->m1_mat_ab[idx][j][i] *= nu_fourth * g_nu_bar[idx];
                
                out->m1_mat_em[idx][i][j] *= nu_fourth * (1. - g_nu[idx]);
                out->m1_mat_em[idx][j][i] *= nu_fourth * (1. - g_nu_bar[idx]);

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
 * Note that there are no double integrals for the computation of the scattering
 * coefficient.
 */
KOKKOS_INLINE_FUNCTION
M1MatrixKokkos2D ComputeDoubleIntegrand(const MyQuadrature* quad, double t,
                                  GreyOpacityParams* grey_pars, const int stim_abs)
{
    const int n = quad->nx;
    double nu, nu_bar;
    double nu_array[n_max];
    M1MatrixKokkos2D out = {0};

    for (int i = 0; i < n; i++)
    {
        nu_array[i]     = t * quad->points[i];
        nu_array[n + i] = t / quad->points[i];
    }

    // compute the neutrino & anti-neutrino distribution function
    double g_nu[total_num_species], g_nu_bar[total_num_species];
    double block_factor_nu[total_num_species],
        block_factor_nu_bar[total_num_species];

    if (grey_pars->opacity_flags.use_pair == 1)
    {
      AddPairKernelsToIntegrand(n, nu_array, grey_pars, &out);
    }
     

    if (grey_pars->opacity_flags.use_brem == 1)
    {
        AddBremKernelsToIntegrand(n, nu_array, grey_pars, &out);
    }

    if ((grey_pars->opacity_flags.use_pair == 1) || (grey_pars->opacity_flags.use_brem == 1))
    {
      WeightNuNuBarReactionsWithDistr(n, nu_array, grey_pars, &out);
    }   


    if (grey_pars->opacity_flags.use_inelastic_scatt == 1)
    {
       AddInelKernelsToIntegrand(n, nu_array, grey_pars, &out);
    }

    /*
    if (grey_pars->opacity_flags.use_abs_em == 1)
    {
      AddBetaReactionToIntegrand(n, nu_array, grey_pars, &out, stim_abs);
    }
    */

            /*
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
            */

            //////////////////////////////////////////////

  AddCommonWeightsToIntegrand(n, nu_array, grey_pars, &out, stim_abs);

  return out;
}


/* Computes the opacities for the M1 code
 *
 */
KOKKOS_INLINE_FUNCTION
M1Opacities 
ComputeM1OpacitiesGenericFormalism(const MyQuadrature* quad_1d, const MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params, const int stim_abs)
{

    // compute some constants
    static const double four_pi_hc3 =
        (4. * kPi) / (kHClight * kHClight * kHClight); // [MeV^-3 cm^-3]
    static const double four_pi_hc3_sqr =
        four_pi_hc3 * four_pi_hc3; // [MeV^-6 cm^-6]

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

    // @TODO: choose this appropriately
    //const double s_pair = 0.5 * 4.364 * temp;
    const double s_pair = temp * (FDI_p4(eta_e) / FDI_p3(eta_e) + FDI_p4(-eta_e) / FDI_p3(-eta_e));
    const double s_nux  = 1.5 * temp;
    double s_beta[total_num_species] = {0}, s_iso[total_num_species] = {0};

    s_beta[id_nue] = temp * FDI_p5(eta_e) / FDI_p4(eta_e);
    s_beta[id_anue] = temp * FDI_p5(-eta_e) / FDI_p4(-eta_e);
    
    for (int idx = 0; idx < total_num_species; idx++)
    {
        s_iso[idx] = (n[idx] > THRESHOLD_N) ? (J[idx] / n[idx]) : s_nux;
    }


    MyQuadratureIntegrand iso_integrals = {0};
    if (my_grey_opacity_params->opacity_flags.use_iso == 1) {
        double out_iso[total_num_species][n_max];
        Scattering1DIntegrand(quad_1d, my_grey_opacity_params, s_iso, out_iso);
        iso_integrals = GaussLegendreIntegrate1DMatrix(quad_1d, total_num_species, out_iso, s_iso);
    }

    MyQuadratureIntegrand beta_n_em_integrals = {0}, beta_j_em_integrals = {0};
    MyQuadratureIntegrand beta_n_abs_integrals = {0}, beta_j_abs_integrals = {0};

    if (my_grey_opacity_params->opacity_flags.use_abs_em == 1) {
        double out_beta_em[total_num_species][n_max], out_beta_ab[total_num_species][n_max];
        Beta1DIntegrand(quad_1d, my_grey_opacity_params, s_beta, out_beta_em, out_beta_ab, stim_abs);
        GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_beta_em, s_beta, &beta_n_em_integrals, &beta_j_em_integrals);
        GaussLegendreIntegrate1DMatrixOnlyNumber(quad_1d, 2, out_beta_ab, s_beta, &beta_n_abs_integrals, &beta_j_abs_integrals);
    }


    MyQuadratureIntegrand n_integrals_2d = {0};
    MyQuadratureIntegrand e_integrals_2d = {0};
    M1MatrixKokkos2D out = ComputeDoubleIntegrand(quad_2d, s_pair, my_grey_opacity_params, stim_abs);
    GaussLegendreIntegrate2DMatrixForM1Coeffs(quad_2d, &out, s_pair, &n_integrals_2d, &e_integrals_2d);

    M1Opacities m1_opacities = {0};

    /* Set all opacities to zero. They'll be left as 0 if the neutrino
    number/energy density is too low (to avoid inf/nan, since the number/energy
    density appears in the denominator). Note that emissivities do no need this
    precaution (and it also make sense physically: you can produce neutrinos
    even if there are none to start with). */
    for (int idx = 0; idx < total_num_species; idx++)
    {
        m1_opacities.kappa_0_a[idx] = 0.;
        m1_opacities.kappa_a[idx]   = 0.;
        m1_opacities.kappa_s[idx]   = 0.;
    }

    /* Electron neutrinos */
    m1_opacities.eta_0[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[0] + beta_n_em_integrals.integrand[id_nue]);
    m1_opacities.eta[id_nue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[0] + beta_j_em_integrals.integrand[id_nue]);
    if (n[id_nue] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_nue] =
            four_pi_hc3 / (kClight * n[id_nue]) *
            (four_pi_hc3 * n_integrals_2d.integrand[4] +
             beta_n_abs_integrals.integrand[id_nue]);
    }
    if (J[id_nue] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_nue] =
            n[id_nue] == 0. ? 0. :
                              four_pi_hc3 / (kClight * J[id_nue]) *
                                  (four_pi_hc3 * e_integrals_2d.integrand[4] +
                                   beta_j_abs_integrals.integrand[id_nue]);
        m1_opacities.kappa_s[id_nue] =
            four_pi_hc3 / (kClight * J[id_nue]) * iso_integrals.integrand[id_nue];
    }

    /* Electron anti-neutrinos */
    m1_opacities.eta_0[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * n_integrals_2d.integrand[1] + beta_n_em_integrals.integrand[id_anue]);
    m1_opacities.eta[id_anue] =
        four_pi_hc3 *
        (four_pi_hc3 * e_integrals_2d.integrand[1] + beta_j_em_integrals.integrand[id_anue]);
    if (n[id_anue] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_anue] =
            four_pi_hc3 / (kClight * n[id_anue]) *
            (four_pi_hc3 * n_integrals_2d.integrand[5] +
             beta_n_abs_integrals.integrand[id_anue]);
    }
    if (J[id_anue] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_anue] =
            four_pi_hc3 / (kClight * J[id_anue]) *
            (four_pi_hc3 * e_integrals_2d.integrand[5] +
             beta_j_abs_integrals.integrand[id_anue]);
        m1_opacities.kappa_s[id_anue] =
            four_pi_hc3 / (kClight * J[id_anue]) * iso_integrals.integrand[id_anue];
    }

    /* Heavy neutrinos */
    m1_opacities.eta_0[id_nux] = four_pi_hc3_sqr * n_integrals_2d.integrand[2];
    m1_opacities.eta[id_nux]   = four_pi_hc3_sqr * e_integrals_2d.integrand[2];
    if (n[id_nux] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_nux] = four_pi_hc3_sqr /
                                         (kClight * n[id_nux]) *
                                         n_integrals_2d.integrand[6];
    }
    if (J[id_nux] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_nux] = four_pi_hc3_sqr / (kClight * J[id_nux]) *
                                       e_integrals_2d.integrand[6];
        m1_opacities.kappa_s[id_nux] =
            four_pi_hc3 / (kClight * J[id_nux]) * iso_integrals.integrand[id_nux];
    }

    /* Heavy anti-neutrinos */
    m1_opacities.eta_0[id_anux] = four_pi_hc3_sqr * n_integrals_2d.integrand[3];
    m1_opacities.eta[id_anux]   = four_pi_hc3_sqr * e_integrals_2d.integrand[3];
    if (n[id_anux] > THRESHOLD_N)
    {
        m1_opacities.kappa_0_a[id_anux] =
            n[id_anux] == 0. ? 0. :
                               four_pi_hc3_sqr / (kClight * n[id_anux]) *
                                   n_integrals_2d.integrand[7];
    }
    if (J[id_anux] > THRESHOLD_J)
    {
        m1_opacities.kappa_a[id_anux] = four_pi_hc3_sqr /
                                        (kClight * J[id_anux]) *
                                        e_integrals_2d.integrand[7];

        m1_opacities.kappa_s[id_anux] =
            four_pi_hc3 / (kClight * J[id_anux]) * iso_integrals.integrand[id_anux];
    }

    return m1_opacities;
}

KOKKOS_INLINE_FUNCTION
M1Opacities 
ComputeM1OpacitiesNotStimulated(MyQuadrature* quad_1d, MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params) {
    return ComputeM1OpacitiesGenericFormalism(quad_1d, quad_2d, my_grey_opacity_params, 0);
                                }

KOKKOS_INLINE_FUNCTION
M1Opacities 
ComputeM1Opacities(const MyQuadrature* quad_1d, const MyQuadrature* quad_2d,
                                GreyOpacityParams* my_grey_opacity_params) {
    return ComputeM1OpacitiesGenericFormalism(quad_1d, quad_2d, my_grey_opacity_params, 1);
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
inline
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
