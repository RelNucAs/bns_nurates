//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.hpp
//  \brief contains pair kernels and associated helper functions

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"

#ifdef KOKKOS_FLAG
#define func_tgamma(x) Kokkos::tgamma(x)
#else
#define func_tgamma(x) tgamma(x)
#endif

#ifdef GSL_INCLUDES_H_
#endif // GSL_INCLUDES_H_


KOKKOS_INLINE_FUNCTION
void PairPsi(const int l, const BS_REAL y, const BS_REAL z, const BS_REAL eta,
             BS_REAL* psi_out)
{
    BS_ASSERT(l == 0);

    BS_ASSERT(isfinite(y) and y >= 0.);
    BS_ASSERT(isfinite(y) and z >= 0.);
    BS_ASSERT(isfinite(eta));

    /* The check on eta is necessary because, if eta is very large, the electron
    phase space is full and the reaction is suppressed. The kernel should in
    principle go to 0 automatically, however in some cases numerical
    cancellation results in small but negative (and so unphysical) rates. This
    check fixes that. */
    /* TODO: The threshold of 200 is somewhat arbitrary, can we find some better
    motivated number? */
    if (eta - y - z > 200.)
    {
        psi_out[0] = 0.;
        psi_out[1] = 0.;
    }
    else
    {
        const BS_REAL FDI_p1_emy  = FDI_p1(eta - y);
        const BS_REAL FDI_p1_emz  = FDI_p1(eta - z);
        const BS_REAL FDI_p1_epy  = FDI_p1(eta + y);
        const BS_REAL FDI_p1_epz  = FDI_p1(eta + z);
        const BS_REAL FDI_p2_emy  = FDI_p2(eta - y);
        const BS_REAL FDI_p2_emz  = FDI_p2(eta - z);
        const BS_REAL FDI_p2_epy  = FDI_p2(eta + y);
        const BS_REAL FDI_p2_epz  = FDI_p2(eta + z);
        const BS_REAL FDI_p3_e    = FDI_p3(eta);
        const BS_REAL FDI_p3_emy  = FDI_p3(eta - y);
        const BS_REAL FDI_p3_emz  = FDI_p3(eta - z);
        const BS_REAL FDI_p3_epy  = FDI_p3(eta + y);
        const BS_REAL FDI_p3_epz  = FDI_p3(eta + z);
        const BS_REAL FDI_p3_emyz = FDI_p3(eta - y - z);
        const BS_REAL FDI_p3_epyz = FDI_p3(eta + y + z);
        const BS_REAL FDI_p4_e    = FDI_p4(eta);
        const BS_REAL FDI_p4_emy  = FDI_p4(eta - y);
        const BS_REAL FDI_p4_emz  = FDI_p4(eta - z);
        const BS_REAL FDI_p4_epy  = FDI_p4(eta + y);
        const BS_REAL FDI_p4_epz  = FDI_p4(eta + z);
        const BS_REAL FDI_p4_emyz = FDI_p4(eta - y - z);
        const BS_REAL FDI_p4_epyz = FDI_p4(eta + y + z);
        const BS_REAL FDI_p5_emy  = FDI_p5(eta - y);
        const BS_REAL FDI_p5_emz  = FDI_p5(eta - z);
        const BS_REAL FDI_p5_epy  = FDI_p5(eta + y);
        const BS_REAL FDI_p5_epz  = FDI_p5(eta + z);
        const BS_REAL FDI_p5_emyz = FDI_p5(eta - y - z);
        const BS_REAL FDI_p5_epyz = FDI_p5(eta + y + z);

        const BS_REAL x0  = 20 * FDI_p4_emy;
        const BS_REAL x1  = 20 * FDI_p4_epz;
        const BS_REAL x2  = 120 * FDI_p2_emy - 120 * FDI_p2_epz;
        const BS_REAL x3  = -40 * FDI_p3_emy + 40 * FDI_p3_epz;
        const BS_REAL x4  = 40 * FDI_p3_e;
        const BS_REAL x5  = 40 * FDI_p3_emyz - x4;
        const BS_REAL x6  = -20 * FDI_p4_e;
        const BS_REAL x7  = 20 * FDI_p4_emyz + x6;
        const BS_REAL x8  = -40 * FDI_p3_epyz + x4;
        const BS_REAL x9  = 20 * FDI_p4_epyz + x6;
        const BS_REAL x10 = -4 * FDI_p5_emy + 4 * FDI_p5_emyz - 4 * FDI_p5_emz +
                            4 * FDI_p5_epy - 4 * FDI_p5_epyz + 4 * FDI_p5_epz;
        const BS_REAL x11 = 20 * FDI_p4_epy;
        const BS_REAL x12 = 20 * FDI_p4_emz;
        const BS_REAL x13 = -40 * FDI_p3_emz + 40 * FDI_p3_epy;
        const BS_REAL x14 = 120 * FDI_p2_emz - 120 * FDI_p2_epy;

        const BS_REAL aux = 15 * POW2(y * z);

        psi_out[0] =
            x10 +
            y * (-x0 + x1 + x7 +
                 y * (x3 + x5 +
                      z * (x2 + z * (-120 * FDI_p1_emy + 120 * FDI_p1_epz))) +
                 z * (80 * FDI_p3_emy - 80 * FDI_p3_epz - x2 * z)) +
            z * (x0 - x1 + x9 + z * (x3 + x8));

        psi_out[1] =
            x10 + y * (-x11 + x12 + x9 + y * (x13 + x8)) +
            z * (x11 - x12 + x7 +
                 y * (80 * FDI_p3_emz - 80 * FDI_p3_epy - x14 * y) +
                 z * (x13 + x5 +
                      y * (x14 + y * (-120 * FDI_p1_emz + 120 * FDI_p1_epy))));

        psi_out[0] /= aux;
        psi_out[1] /= aux;
    }
}

/* Calculate Phi_l(y,z) from Eqn. (10) of Pons et. al. (1998)
 *
 * Inputs:
 *      l:            mode number
 *      omega:        neutrino energy [MeV]
 *      omega_prime:  anti-neutrino energy [MeV]
 *      temp:         temperature [MeV]
 *      e_x:          neutrino species type (0: elentron, 1: mu/tau)
 *
 * Output:
 *      Phi_l(y,z) = (G^2 temp^2)/(pi (1 - e^{y+z})) [alpha1 Psi_l(y,z) + alpha2
 * Psi_l(z,y)]
 */
KOKKOS_INLINE_FUNCTION
void PairPhi(const BS_REAL omega, const BS_REAL omega_prime, const int l,
             const BS_REAL eta, const BS_REAL T, BS_REAL* phi_out)
{
    const BS_REAL y = omega / T;
    const BS_REAL z = omega_prime / T;

    const BS_REAL aux = kBS_Pair_Phi * POW2(T) / (1. - SafeExp(y + z));

    BS_REAL pair_psi[2] = {0.};

    PairPsi(l, y, z, eta, pair_psi);

    phi_out[0] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[0] +
                  POW2(kBS_Pair_Alpha2_0) * pair_psi[1]) *
                 aux;
    phi_out[1] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[1] +
                  POW2(kBS_Pair_Alpha2_0) * pair_psi[0]) *
                 aux;
    phi_out[2] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[0] +
                  POW2(kBS_Pair_Alpha2_1) * pair_psi[1]) *
                 aux;
    phi_out[3] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[1] +
                  POW2(kBS_Pair_Alpha2_1) * pair_psi[0]) *
                 aux;
}

KOKKOS_INLINE_FUNCTION
MyKernelOutput PairKernels(const MyEOSParams* eos_pars,
                           const PairKernelParams* kernel_pars)
{
    // EOS specific parameters
    const BS_REAL T   = eos_pars->temp;
    const BS_REAL eta = eos_pars->mu_e / T;

    // kernel specific parameters
    const BS_REAL omega       = kernel_pars->omega;
    const BS_REAL omega_prime = kernel_pars->omega_prime;

    BS_REAL pair_phi[4] = {0.};

    PairPhi(omega, omega_prime, 0, eta, T, pair_phi);

    MyKernelOutput pair_kernel;

    pair_kernel.em[id_nue]  = 0.5 * pair_phi[0];
    pair_kernel.em[id_anue] = 0.5 * pair_phi[1];
    pair_kernel.em[id_nux]  = 0.5 * pair_phi[2];
    pair_kernel.em[id_anux] = 0.5 * pair_phi[3];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        pair_kernel.abs[idx] =
            SafeExp((omega + omega_prime) / T) * pair_kernel.em[idx];
    }

    return pair_kernel;
}

KOKKOS_INLINE_FUNCTION
void PairKernels(const MyEOSParams* eos_pars,
                 const PairKernelParams* kernel_pars, MyKernelOutput* out_for,
                 MyKernelOutput* out_inv)
{
    *out_for = PairKernels(eos_pars, kernel_pars);

    out_inv->em[id_nue]  = out_for->em[id_anue];
    out_inv->em[id_anue] = out_for->em[id_nue];
    out_inv->em[id_nux]  = out_for->em[id_anux];
    out_inv->em[id_anux] = out_for->em[id_nux];

    out_inv->abs[id_nue]  = out_for->abs[id_anue];
    out_inv->abs[id_anue] = out_for->abs[id_nue];
    out_inv->abs[id_nux]  = out_for->abs[id_anux];
    out_inv->abs[id_anux] = out_for->abs[id_nux];
}

KOKKOS_INLINE_FUNCTION
void PairKernelsTable(const int n, const BS_REAL* nu_array,
                      GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput pair_1, pair_2;

    grey_pars->kernel_pars.pair_kernel_params.cos_theta = 1.;
    grey_pars->kernel_pars.pair_kernel_params.filter    = 0.;
    grey_pars->kernel_pars.pair_kernel_params.lmax      = 0;
    grey_pars->kernel_pars.pair_kernel_params.mu        = 1.;
    grey_pars->kernel_pars.pair_kernel_params.mu_prime  = 1.;

    for (int i = 0; i < n; ++i)
    {
        grey_pars->kernel_pars.pair_kernel_params.omega = nu_array[i];

        for (int j = i; j < n; ++j)
        {
            grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[j];

            PairKernels(&grey_pars->eos_pars,
                        &grey_pars->kernel_pars.pair_kernel_params, &pair_1,
                        &pair_2);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] = pair_1.em[idx];
                out->m1_mat_em[idx][j][i] = pair_2.em[idx];

                out->m1_mat_ab[idx][i][j] = pair_1.abs[idx];
                out->m1_mat_ab[idx][j][i] = pair_2.abs[idx];
            }
        }
    }

    return;
}
