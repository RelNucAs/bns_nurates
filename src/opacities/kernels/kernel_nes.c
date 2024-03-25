//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nes.c
//  \brief contains kernels for inelastic neutrino scattering
//         on electrons and positrons
//
// Computation of inelastic neutrino-electron and neutrino-positron
// scattering using Eq. (43) from Mezzacappa & Bruenn, ApJ v.410, p.740 (1993)
// https://ui.adsabs.harvard.edu/abs/1993ApJ...410..740M/abstract

#include "math.h"

#include "kernels.h"
#include "bns_nurates.h"
#include "functions.h"
#include "functions.h"
#include "constants.h"

// Numerical constants
//---------------------------------------------------------------------------------------------------------------------
static const double kTaylorSeriesEpsilon = 1e-3;


// Physical constats
//---------------------------------------------------------------------------------------------------------------------
static const double kSinWeinbergThetaSquared = 0.2229;
static const double kNes =
    2. * kGf0 * kGf0 * kClight * kHbarClight * kHbarClight / (3. * kPi);
static const double kBPlus =
    (2. * kSinWeinbergThetaSquared + 1.) * (2. * kSinWeinbergThetaSquared + 1.);
static const double kBMinus =
    (2. * kSinWeinbergThetaSquared - 1.) * (2. * kSinWeinbergThetaSquared - 1.);
static const double kBZero =
    4. * kSinWeinbergThetaSquared * kSinWeinbergThetaSquared;

// Calculates in a safe way the exponential
double KernelExpFunc(double x)
{
    const double exp_  = SafeExp(-fabs(x));
    const int is_x_neg = signbit(x);

    return (is_x_neg - ! is_x_neg * exp_) / (1. - exp_ + (x == 0.));
}

// Saves the expression that are needed for the unapproximated integral
void ComputeFDIForInelastic(double w, double wp, double eta, double* fdi_diff_w,
                            double* fdi_diff_abs)
{
    double abs_val = fabs(w - wp);

    fdi_diff_w[0] = FDI_p1(eta - wp) - FDI_p1(eta - w);
    fdi_diff_w[1] = FDI_p2(eta - wp) - FDI_p2(eta - w);
    fdi_diff_w[2] = FDI_p3(eta - wp) - FDI_p3(eta - w);
    fdi_diff_w[3] = FDI_p4(eta - wp) - FDI_p4(eta - w);
    fdi_diff_w[4] = FDI_p5(eta - wp) - FDI_p5(eta - w);

    fdi_diff_abs[0] = FDI_p3(eta) - FDI_p3(eta - abs_val);
    fdi_diff_abs[1] = FDI_p4(eta) - FDI_p4(eta - abs_val);
    fdi_diff_abs[2] = FDI_p5(eta) - FDI_p5(eta - abs_val);
}

// Functions that calculates the kernel integral
//=========================================================================================================================================

// Not approximated out kernel integral
double MezzacappaIntOut(double w, double wp, double x, double y, int sign,
                        double b1, double b2, double* fdi_diff_w,
                        double* fdi_diff_abs)
{
    return (((b1 + b2) * (sign * fdi_diff_abs[2] - fdi_diff_w[4]) *
                 0.2 // All G5 terms

             - b1 * (w + wp) *
                   (

                       fdi_diff_w[3]

                       + 2. * ((w + wp) * fdi_diff_w[2] +
                               3. * w * wp * fdi_diff_w[1])

                           ) // G4(eta - wp) + G3(eta - wp) term

             + sign * ((b1 * x - b2 * y) * fdi_diff_abs[1] +
                       2. * (b1 * x * x + b2 * y * y) * fdi_diff_abs[0])) /
                (w * w * wp * wp)

            - 6. * b1 * fdi_diff_w[0]);
}

// Taylor expansion in the lowest energy of the function MezzacappaIntOut and
// MezzacappaIntIn
double MezzacappaIntOneEnergy(double x, double y, int sign, double b1,
                              double b2, const double* fdis)
{
    return -sign * y *
           (2. * (b1 + b2) * fdis[0]

            + (b1 * (y + 4. * x) + 3. * b2 * y) * fdis[1]

            + (b1 * (3. * y - 4. * x) + b2 * y) * fdis[2]

            + ((b1 + 6. * b2) * y * y * 0.2 + b1 * x * y + 2. * b1 * x * x) *
                  fdis[3]

            -
            ((6. * b1 + b2) * y * y * 0.2 - 3. * b1 * x * y + 2. * b1 * x * x) *
                fdis[4]) /
           (x * x);
}

// Taylor expansion in both energies of the function MezzacappaIntOut and
// MezzacappaIntIn
double MezzacappaIntTwoEnergies(double w, double wp, double x, double y,
                                double b1, double b2, const double* fdis)
{
    return y * (w - wp) *
           ((b1 + b2) *
                (4. * fdis[2]

                 + fdis[0] * (11. * y * y - 25. * x * y + 20. * x * x) / 30.)

            + (b1 - b2) * (2. * x - y) * fdis[1]) /
           (x * x);
}

// Functions for kernel calculation
//===================================================================================================

// Calculates and saves the neutrino electron scattering in and out kernel for
// every neutrino species
MyKernelOutput NESKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const double t = eos_params->temp, w = kernel_params->omega / t,
                 wp = kernel_params->omega_prime / t, x = fmax(w, wp),
                 y = fmin(w, wp), eta_e = eos_params->mu_e / t,
                 exp_factor           = KernelExpFunc(wp - w),
                 exp_factor_exchanged = KernelExpFunc(w - wp);

    MyKernelOutput output;

    if (y > eta_e * kTaylorSeriesEpsilon)
    {
        double fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_e, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] = kNes * t * t *
                             MezzacappaIntOut(w, wp, x, y, sign, kBPlus, kBZero,
                                              fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntOut(w, wp, x, y, sign, kBZero, kBPlus, fdi_diff_w,
                             fdi_diff_abs);
        output.abs[id_nux] = kNes * t * t *
                             MezzacappaIntOut(w, wp, x, y, sign, kBMinus,
                                              kBZero, fdi_diff_w, fdi_diff_abs);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntOut(w, wp, x, y, sign, kBZero, kBMinus, fdi_diff_w,
                             fdi_diff_abs);
    }
    else if (x > eta_e * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const double fdis[5] = {
            FDI_p2(eta_e - x) - FDI_p2(eta_e), FDI_p1(eta_e - x), FDI_p1(eta_e),
            FDI_0(eta_e - x) + y * FermiDistr(0., 1., eta_e - x) / 6.,
            FDI_0(eta_e) - y * FermiDistr(0., 1., eta_e) / 6.};

        output.abs[id_nue] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBPlus, kBZero, fdis);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBZero, kBPlus, fdis);
        output.abs[id_nux] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBMinus, kBZero, fdis);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBZero, kBMinus, fdis);
    }
    else
    {
        const double fdis[3] = {FermiDistr(0., 1., eta_e), FDI_0(eta_e),
                                FDI_p1(eta_e)};

        output.abs[id_nue] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBPlus, kBZero, fdis);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBZero, kBPlus, fdis);
        output.abs[id_nux] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBMinus, kBZero, fdis);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBZero, kBMinus, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates and saves the neutrino positron scattering in and out kernel for
// every neutrino species
MyKernelOutput NPSKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const double t = eos_params->temp, w = kernel_params->omega / t,
                 wp = kernel_params->omega_prime / t, x = fmax(w, wp),
                 y = fmin(w, wp), eta_p = -eos_params->mu_e / t,
                 exp_factor           = KernelExpFunc(wp - w),
                 exp_factor_exchanged = KernelExpFunc(w - wp);

    MyKernelOutput output;

    if (y > eta_p * kTaylorSeriesEpsilon)
    {
        double fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_p, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] = kNes * t * t *
                             MezzacappaIntOut(w, wp, x, y, sign, kBZero, kBPlus,
                                              fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntOut(w, wp, x, y, sign, kBPlus, kBZero, fdi_diff_w,
                             fdi_diff_abs);
        output.abs[id_nux] =
            kNes * t * t *
            MezzacappaIntOut(w, wp, x, y, sign, kBZero, kBMinus, fdi_diff_w,
                             fdi_diff_abs);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntOut(w, wp, x, y, sign, kBMinus, kBZero, fdi_diff_w,
                             fdi_diff_abs);
    }
    else if (x > eta_p * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const double fdis[5] = {
            FDI_p2(eta_p - x) - FDI_p2(eta_p), FDI_p1(eta_p - x), FDI_p1(eta_p),
            FDI_0(eta_p - x) + y * FermiDistr(0., 1., eta_p - x) / 6.,
            FDI_0(eta_p) - y * FermiDistr(0., 1., eta_p) / 6.};

        output.abs[id_nue] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBZero, kBPlus, fdis);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBPlus, kBZero, fdis);
        output.abs[id_nux] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBZero, kBMinus, fdis);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntOneEnergy(x, y, sign, kBMinus, kBZero, fdis);
    }
    else
    {
        const double fdis[3] = {FermiDistr(0., 1., eta_p), FDI_0(eta_p),
                                FDI_p1(eta_p)};

        output.abs[id_nue] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBZero, kBPlus, fdis);
        output.abs[id_anue] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBPlus, kBZero, fdis);
        output.abs[id_nux] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBZero, kBMinus, fdis);
        output.abs[id_anux] =
            kNes * t * t *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBMinus, kBZero, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates the full in and out kernels
MyKernelOutput InelasticScattKernels(InelasticScattKernelParams* kernel_params,
                                     MyEOSParams* eos_params)
{
    MyKernelOutput nes_kernel = NESKernels(kernel_params, eos_params);
    MyKernelOutput nps_kernel = NPSKernels(kernel_params, eos_params);

    MyKernelOutput tot_kernel = {0};

    for (int idx = 0; idx < total_num_species; idx++)
    {
        tot_kernel.em[idx]  = nes_kernel.em[idx] + nps_kernel.em[idx];
        tot_kernel.abs[idx] = nes_kernel.abs[idx] + nps_kernel.abs[idx];
    }

    return tot_kernel;
}

void InelasticKernelsTable(const int n, double* nu_array,
                           GreyOpacityParams* grey_pars, M1Matrix* out)
{
    MyKernelOutput inel_1, inel_2;

    InelasticScattKernelParams inelastic_pars =
        grey_pars->kernel_pars.inelastic_kernel_params;
    for (int i = 0; i < n; i++)
    {

        for (int j = i; j < n; j++)
        {

            // compute the pair kernels
            inelastic_pars.omega       = nu_array[i];
            inelastic_pars.omega_prime = nu_array[j];
            inel_1 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);

            inelastic_pars.omega       = nu_array[j];
            inelastic_pars.omega_prime = nu_array[i];
            inel_2 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);


            for (int idx = 0; idx < total_num_species; idx++)
            {
                out->m1_mat_em[idx][i][j] = inel_1.em[idx];
                out->m1_mat_em[idx][j][i] = inel_2.em[idx];

                out->m1_mat_ab[idx][i][j] = inel_1.abs[idx];
                out->m1_mat_ab[idx][j][i] = inel_2.abs[idx];
            }
        }
    }

    return;
}