//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.hpp
//  \brief contains pair kernels and associated helper functions

#include <math.h>
#include <assert.h>
#include <stddef.h>
#include "functions.hpp"
#include "constants.hpp"

#ifdef KOKKOS_FLAG
#define func_tgamma(x) Kokkos::tgamma(x)
#else
#define func_tgamma(x) tgamma(x)
#endif

#ifdef GSL_INCLUDES_H_
#endif // GSL_INCLUDES_H_

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define POW0(X) ((1))
#define POW1(X) ((X))
#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))
#define POW5(X) ((X) * (X) * (X) * (X) * (X))

//===========================================
// Optimized functions to speed the code up
//===========================================

/* Compute T_l(alpha) as defined in Appendix B of Pons et. al. (1998)
 * upto a given accuracy
 *
 * T_l(alpha) = sum_{n=1..inf} (-1)^{n+1} e^{-n alpha} / n^l
 *
 * This is well-defined for alpha > 0 and l >= 1.
 *
 * Inputs:
 *    l:          mode number (>= 1)
 *    alpha:      function argument (> 0)
 *    tolerance:  a measure of accuracy of the computation (1e-6 is a good
 * number)
 *
 * Output:
 *    T_l(alpha) = sum_{n=1..inf} (-1)^{n+1} e^{-n alpha} / n^l
 */
KOKKOS_INLINE_FUNCTION
double PairTWithAlpha0(int l)
{
    switch (l)
    {
    case 1:
        return 0.693147180559945309417; // log(2.0)
        break;

    // else : pow(2., -l) * (pow(2., l) - 2.) * gsl_sf_zeta_int(l);  // Computed
    // in Mathematica
    case 2:
        return 0.822467033424113218236; // kPi * kPi / 12.
        break;

    case 3:
        return 0.901542677369695714050; // 3. * zeta(3.) / 4.
        break;

    case 4:
        return 0.947032829497245917577; // 7. * kPi * kPi * kPi * kPi / 720.
        break;

    case 5:
        return 0.972119770446909305936; // 15. * zeta(5.) / 16.
        break;

    case 6:
        return 0.985551091297435104098; // 31. * pow(kPi, 6) / 30240.
        break;

    default:
        printf(
            "PairTWithAlpha0 (kernel_pair.c): l = %d must be within 1 and 6\n",
            l);
        return -42.;
    }
}

KOKKOS_INLINE_FUNCTION
double PairTFitted(int l, double alpha)
{
    const double a[6][13] = {
        {0.0009591818519410957, -0.007025445963539029, 0.02413571687574133,
         -0.05283857730366366, 0.08565184150347487, -0.1148974833969069,
         0.14001269129489452, -0.16612181245182328, 0.19993311519859502,
         -0.2499952138069261, 0.333333164848656, -0.4999999980100303,
         0.9999999999994182},
        {
            9.806117191850437e-05,
            -0.0007303693955082037,
            0.0025729046032377245,
            -0.005853104482681513,
            0.010060162627142736,
            -0.014705257422400938,
            0.02014664666291523,
            -0.027727326966761626,
            0.039993773761235564,
            -0.062499552625431457,
            0.11111109531105254,
            -0.24999999981287394,
            0.999999999999945,
        },
        {9.83660257000452e-06, -7.455399425901051e-05, 0.00026968069130985794,
         -0.0006390877179467623, 0.0011691212967386872, -0.0018705673951822446,
         0.0028917745803742146, -0.0046250335703596395, 0.007999430192056203,
         -0.015624958913956509, 0.03703703558177793, -0.12499999998270364,
         0.999999999999995},
        {9.63655550546077e-07, -7.4501714680668755e-06, 2.775504190088471e-05,
         -6.878346378088752e-05, 0.00013456162975559815,
         -0.00023676741517985415, 0.00041435638829548566,
         -0.0007711865782589587, 0.0015999477373502997, -0.003906246205609882,
         0.012345678877131755, -0.06249999999837849, 0.9999999999999992},
        {9.440437566332406e-08, -7.417123619081298e-07, 2.838371714509535e-06,
         -7.350534216150028e-06, 1.5399946770677647e-05,
         -2.9877545939654987e-05, 5.931172050750023e-05,
         -0.00012856362411722033, 0.00031999523335656476,
         -0.0009765621364823982, 0.004115226323079645, -0.031249999999776252,
         1.0000000000000002},
        {1.1258637058816464e-08, -8.606136518301832e-08, 3.2215497280909555e-07,
         -8.319730818363639e-07, 1.8034586648497246e-06, -3.792399721761286e-06,
         8.496951883750312e-06, -2.1433700743647905e-05, 6.400013840911159e-05,
         -0.00024414064372257427, 0.001371742113499251, -0.015625000000008864,
         1.0000000000000002}};

    double fit_poly = 1.;

    double sum       = 0.;
    double exp_alpha = exp(-alpha);
    double exp_pow   = 1.;

    if (alpha < 25.)
    {
        if (alpha == 0)
            return PairTWithAlpha0(l);

        for (int i = 0; i < 13; i++)
        {
            sum += exp_pow * a[l - 1][12 - i];
            exp_pow = exp_pow * exp_alpha;
        }

        fit_poly = sum;
    }


    return exp_alpha * fit_poly;
}

/* Compute F_k(eta,x1) as defined in Appendix B of Pons et. al. (1998)
 *
 * Inputs:
 *    k:    an integer
 *    eta:  double
 *    x1:   double
 *
 * Output:
 *  F_k:    as defined below
 *
 * eta < 0:         F_k(eta,x1) = k! [T_{k+1}(-eta) - sum_{l=0..k}
 * T_{k+1-l}(x1-eta) x1^l / l!] 0 <= eta <= x1:  F_k(eta,x1) = eta^{k+1}/(k+1)
 *                              + k! [2 sum_{l=0..int((k-1)/2)} T_{2l+2}(0)
 * eta^{k-1-2l}/(k-1-2l)!
 *                                    + (-1)^k T_{k+1}(eta)
 *                                    - sum_{l=0..k} T_{k+1-l}(x1-eta) x1^l /
 * l!] x1 < eta:        F_k(eta,x1) = x1^{k+1}/(k+1)
 *                              + k! [(-1)^k T_{k+1}(eta)
 *                                     - sum_{l=0..k} (-1)^{k-l}
 * T_{k+1-l}(eta-x1) x1^l / l!]
 */
KOKKOS_INLINE_FUNCTION
double PairFOptimized(int k, double eta, double x1)
{

    double result = 0.;
    int fact      = 1;
    double tmp    = 1.;

    if (eta < 0.)
    {
        double sum = 0.;
        sum += PairTFitted(k + 1, x1 - eta); // l=0
        for (int l = 1; l <= k; l++)
        {
            tmp  = tmp * x1;
            fact = fact * l;
            sum += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
        }
        result = fact * (PairTFitted(k + 1, fabs(eta)) - sum);
    }
    else if (k != 0 && 0. <= eta && eta <= x1)
    {

        double eta_to_k_minus_one = pow(eta, k - 1.);
        double eta_to_k_plus_one  = eta_to_k_minus_one * eta * eta;

        double sum = 0.;
        tmp        = 1.;
        sum += PairTWithAlpha0(2) / func_tgamma(k); // l=0
        for (int l = 1; l <= (int)((k - 1.) / 2.); l++)
        {
            tmp = tmp / (eta * eta);
            sum += PairTWithAlpha0(2 * l + 2) * tmp /
                   func_tgamma(k - 1 - 2 * l + 1);
        }

        double sum_2 = 0.;
        sum_2 += PairTFitted(k + 1, x1 - eta); // l=0
        fact = 1;
        tmp  = 1.;
        for (int l = 1; l <= k; l++)
        {
            fact = fact * l;
            tmp  = tmp * x1;
            sum_2 += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
        }

        result = eta_to_k_plus_one / (k + 1.) +
                 fact * (2.0 * eta_to_k_minus_one * sum +
                         pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum_2);
    }
    else if (k == 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        sum += PairTFitted(k + 1, x1 - eta);
        for (int l = 1; l <= k; l++)
        {
            fact = fact * l;
            tmp  = tmp * x1;
            sum += PairTFitted(k + 1 - l, x1 - eta) * tmp / fact;
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 fact * (pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum);
    }
    else
    {

        double sum = 0.;
        sum += pow(-1., k) * PairTFitted(k + 1, eta - x1);
        for (int l = 1; l <= k; l++)
        {
            fact = fact * l;
            tmp  = tmp * x1;
            sum +=
                pow(-1., k - l) * PairTFitted(k + 1 - l, eta - x1) * tmp / fact;
        }

        tmp    = tmp * x1;
        result = tmp / (k + 1.) +
                 fact * (pow(-1., k) * PairTFitted(k + 1, fabs(eta)) - sum);
    }

    return result;
}

KOKKOS_INLINE_FUNCTION
void PairPsiOptimized(int l, double y, double z, double eta, double* psi_out)
{
    assert(l == 0);

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
        const double FDI_p1_emy  = FDI_p1(eta - y);
        const double FDI_p1_emz  = FDI_p1(eta - z);
        const double FDI_p1_epy  = FDI_p1(eta + y);
        const double FDI_p1_epz  = FDI_p1(eta + z);
        const double FDI_p2_emy  = FDI_p2(eta - y);
        const double FDI_p2_emz  = FDI_p2(eta - z);
        const double FDI_p2_epy  = FDI_p2(eta + y);
        const double FDI_p2_epz  = FDI_p2(eta + z);
        const double FDI_p3_e    = FDI_p3(eta);
        const double FDI_p3_emy  = FDI_p3(eta - y);
        const double FDI_p3_emz  = FDI_p3(eta - z);
        const double FDI_p3_epy  = FDI_p3(eta + y);
        const double FDI_p3_epz  = FDI_p3(eta + z);
        const double FDI_p3_emyz = FDI_p3(eta - y - z);
        const double FDI_p3_epyz = FDI_p3(eta + y + z);
        const double FDI_p4_e    = FDI_p4(eta);
        const double FDI_p4_emy  = FDI_p4(eta - y);
        const double FDI_p4_emz  = FDI_p4(eta - z);
        const double FDI_p4_epy  = FDI_p4(eta + y);
        const double FDI_p4_epz  = FDI_p4(eta + z);
        const double FDI_p4_emyz = FDI_p4(eta - y - z);
        const double FDI_p4_epyz = FDI_p4(eta + y + z);
        const double FDI_p5_emy  = FDI_p5(eta - y);
        const double FDI_p5_emz  = FDI_p5(eta - z);
        const double FDI_p5_epy  = FDI_p5(eta + y);
        const double FDI_p5_epz  = FDI_p5(eta + z);
        const double FDI_p5_emyz = FDI_p5(eta - y - z);
        const double FDI_p5_epyz = FDI_p5(eta + y + z);

        const double x0  = 20 * FDI_p4_emy;
        const double x1  = 20 * FDI_p4_epz;
        const double x2  = 120 * FDI_p2_emy - 120 * FDI_p2_epz;
        const double x3  = -40 * FDI_p3_emy + 40 * FDI_p3_epz;
        const double x4  = 40 * FDI_p3_e;
        const double x5  = 40 * FDI_p3_emyz - x4;
        const double x6  = -20 * FDI_p4_e;
        const double x7  = 20 * FDI_p4_emyz + x6;
        const double x8  = -40 * FDI_p3_epyz + x4;
        const double x9  = 20 * FDI_p4_epyz + x6;
        const double x10 = -4 * FDI_p5_emy + 4 * FDI_p5_emyz - 4 * FDI_p5_emz +
                           4 * FDI_p5_epy - 4 * FDI_p5_epyz + 4 * FDI_p5_epz;
        const double x11 = 20 * FDI_p4_epy;
        const double x12 = 20 * FDI_p4_emz;
        const double x13 = -40 * FDI_p3_emz + 40 * FDI_p3_epy;
        const double x14 = 120 * FDI_p2_emz - 120 * FDI_p2_epy;

        const double aux = (15 * POW2(y) * POW2(z));

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

    return;
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
void PairPhiOptimized(const double omega, const double omega_prime, int l,
                      double eta, double temp, double* phi_out)
{
    static const double kPairPhi = kGSqr / kPi;

    const double y = omega / temp;
    const double z = omega_prime / temp;

    const double phi_prefactor = kPairPhi * POW2(temp);
    const double phi_denom     = 1. - exp(y + z);

    double pair_psi[2] = {0.};

    PairPsiOptimized(l, y, z, eta, pair_psi);

    phi_out[0] =
        phi_prefactor *
        (POW2(kAlpha1_0) * pair_psi[0] + POW2(kAlpha2_0) * pair_psi[1]) /
        phi_denom;
    phi_out[1] =
        phi_prefactor *
        (POW2(kAlpha1_0) * pair_psi[1] + POW2(kAlpha2_0) * pair_psi[0]) /
        phi_denom;
    phi_out[2] =
        phi_prefactor *
        (POW2(kAlpha1_1) * pair_psi[0] + POW2(kAlpha2_1) * pair_psi[1]) /
        phi_denom;
    phi_out[3] =
        phi_prefactor *
        (POW2(kAlpha1_1) * pair_psi[1] + POW2(kAlpha2_1) * pair_psi[0]) /
        phi_denom;

    return;
}

KOKKOS_INLINE_FUNCTION
MyKernelOutput PairKernelsOptimized(MyEOSParams* eos_pars,
                                    PairKernelParams* kernel_pars)
{
    // EOS specific parameters
    const double eta  = eos_pars->mu_e / eos_pars->temp;
    const double temp = eos_pars->temp;

    // kernel specific parameters
    const double omega       = kernel_pars->omega;
    const double omega_prime = kernel_pars->omega_prime;

    double pair_phi[4] = {0.};

    PairPhiOptimized(omega, omega_prime, 0, eta, temp, pair_phi);

    MyKernelOutput pair_kernel;

    pair_kernel.em[id_nue]  = 0.5 * pair_phi[0];
    pair_kernel.em[id_anue] = 0.5 * pair_phi[1];
    pair_kernel.em[id_nux]  = 0.5 * pair_phi[2];
    pair_kernel.em[id_anux] = 0.5 * pair_phi[3];

    for (int idx = 0; idx < total_num_species; idx++)
    {
        pair_kernel.abs[idx] =
            SafeExp((omega + omega_prime) / temp) * pair_kernel.em[idx];
    }

    return pair_kernel;
}

KOKKOS_INLINE_FUNCTION
void PairKernelsM1Test(MyEOSParams* eos_pars, PairKernelParams* kernel_pars,
                       MyKernelOutput* out_for, MyKernelOutput* out_inv)
{
    *out_for = PairKernelsOptimized(eos_pars, kernel_pars);

    out_inv->em[id_nue]  = out_for->em[id_anue];
    out_inv->em[id_anue] = out_for->em[id_nue];
    out_inv->em[id_nux]  = out_for->em[id_anux];
    out_inv->em[id_anux] = out_for->em[id_nux];

    out_inv->abs[id_nue]  = out_for->abs[id_anue];
    out_inv->abs[id_anue] = out_for->abs[id_nue];
    out_inv->abs[id_nux]  = out_for->abs[id_anux];
    out_inv->abs[id_anux] = out_for->abs[id_nux];

    return;
}

KOKKOS_INLINE_FUNCTION
void PairKernelsTable(const int n, double* nu_array,
                      GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput pair_1, pair_2;

    grey_pars->kernel_pars.pair_kernel_params.cos_theta = 1.;
    grey_pars->kernel_pars.pair_kernel_params.filter    = 0.;
    grey_pars->kernel_pars.pair_kernel_params.lmax      = 0;
    grey_pars->kernel_pars.pair_kernel_params.mu        = 1.;
    grey_pars->kernel_pars.pair_kernel_params.mu_prime  = 1.;

    for (int i = 0; i < n; i++)
    {
        grey_pars->kernel_pars.pair_kernel_params.omega = nu_array[i];

        for (int j = i; j < n; j++)
        {
            grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[j];

            PairKernelsM1Test(&grey_pars->eos_pars,
                              &grey_pars->kernel_pars.pair_kernel_params,
                              &pair_1, &pair_2);

            for (int idx = 0; idx < total_num_species; idx++)
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
