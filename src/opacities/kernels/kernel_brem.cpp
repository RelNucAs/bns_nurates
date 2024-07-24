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


#include <math.h>
#include <assert.h>

#include "kernels.hpp"
#include "constants.hpp"
#include "bns_nurates.hpp"
#include "functions.hpp"

#define POW3(X) ((X) * (X) * (X))

// Constants for the bremsstrahlung reaction
#define kBremXmin 1.0E-10
#define kBremYmin 1.0E-10
#define kBremEtamin 1.0E-10
#define kCa -0.63                              // -1.26 / 2.
#define kSqrtPi 1.7724538509055159             // sqrt(kPi)
#define kPiSquared 9.869604401089358           // kPi * kPi;
#define kFourPiSquared 39.47841760435743       // 4. * kPiSquared
#define kPi2OneEighth 1.1538350678499893       // pow(kPi, 1. / 8.)
#define kPiHalfToFiveHalves 3.0924286813991433 // pow(0.5 * kPi, 2.5)
#define kC4BRT06                                                               \
    0.1409912109375 // (3*5*7*11)/2^{11} / 4 (C in BRT06 Eq.143 divided by four)
// #define kGfBrem 1.1663787E-11 // [MeV^-2]

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
double BremKernelS(double x, double y, double eta_star)
{
    static const double kFiveThirds = 5. / 3.;
    static const double kFiveSixths = 5. / 6.;

    // check non negativity of input quantities
    assert(x >= 0.);
    assert(y >= 0.);
    assert(eta_star >= 0.);

    // prevent singular behaviour
    x        = (x > kBremXmin) ? x : kBremXmin;
    y        = (y > kBremYmin) ? y : kBremYmin;
    eta_star = (eta_star > kBremEtamin) ? eta_star : kBremEtamin;

    // compute non-degenerate approximation, s_nd in Eqn. (45)

    const double s_nd_numerator =
        2. * kSqrtPi * pow(x + 2. - exp(-y / 12.), 1.5) *
        (x * x + 2. * x * y + kFiveThirds * y * y + 1.);
    const double s_nd_denominator = kSqrtPi + pow(kPi2OneEighth + x + y, 4.);
    const double s_nd             = s_nd_numerator / s_nd_denominator;

    // compute degenerate approximation, s_d in Eqn. (46)
    const double u     = sqrt(y / (2. * eta_star)) + 1.e-10;
    const double u2    = u * u;
    const double u_arg = u2 / (2. * sqrt(2. * u2 + 4.));
    double f_u = (1. - kFiveSixths * u * atan(2. / u) + u2 / (3. * (u2 + 4.)) +
                  atan(1. / u_arg) * u_arg / 3.);

    // @TODO: Leonardo check this! Doing this to prevent s_d from being a large
    // negative number
    if (fabs(f_u) < 1.e-14)
    {
        f_u = 1.e-14;
    }

    const double s_d = 3. * kPiHalfToFiveHalves * pow(eta_star, -2.5) *
                       (x * x + kFourPiSquared) * x * f_u /
                       (kFourPiSquared * (1. - SafeExp(-x)));

    /*
    if (s_d < 0.) {
      printf("s_D = %.5e\n"              , s_d);
      printf("one minus exp(-x) = %.5e\n", 1.-SafeExp(-x));
      printf("x = %.5e\n"                , x);
      printf("y = %.5e\n"                , y);
      printf("u = %.5e\n"                , u);
      printf("eta_star = %.5e\n"         , eta_star);
      printf("f_u = %.5e\n"              , f_u);
    } */

    const double pow_x_1_1 = pow(x, 1.1);

    // F, Eqn. (50)
    const double f_denominator = (3. + pow(x - 1.2, 2.) + pow(x, -4.)) *
                                 (1. + eta_star * eta_star) * (1. + pow(y, 4.));
    const double f_brem = 1. + 1. / f_denominator; // Eq.(50)

    // G, Eqn. (50)
    const double g_brem = 1. - 0.0044 * pow_x_1_1 * y /
                                   (0.8 + 0.06 * pow(y, 1.05)) *
                                   sqrt(eta_star) / (eta_star + 0.2);

    // h and C, Eqn. (50)
    const double h_brem = 0.1 * eta_star / (2.39 + 0.1 * pow(eta_star, 1.1));
    const double c_brem = 1.1 * pow_x_1_1 * h_brem /
                          (2.3 + h_brem * pow(x, 0.93) + 0.0001 * pow(x, 1.2)) *
                          30. / (30. + 0.005 * pow(x, 2.8));

    // p, Eqn. (50)
    const double p_brem = 0.67 + 0.18 * pow(y, 0.4);

    // interpolated formula for s in Eqn. (49)
    const double s_brem =
        pow(pow(s_nd, -p_brem) + pow(s_d, -p_brem), -1. / p_brem) * f_brem *
        (1. + c_brem * g_brem);

    assert(s_brem >= 0.);

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
double BremKernelG(double y, double eta_star)
{
    // check non negativity of input quantities
    assert(y >= 0.);
    assert(eta_star >= 0.);

    // prevent singular behavior
    y        = (y > kBremYmin) ? y : kBremYmin;
    eta_star = (eta_star > kBremEtamin) ? eta_star : kBremEtamin;

    // alpha_1, Eqn. (53)
    const double y2            = y * y;
    const double eta_star_inv  = 1. / eta_star;
    const double alpha_1_denom = 25. * y2 + 1.;

    const double alpha_1 =
        (0.5 + eta_star_inv) / (1. + eta_star_inv) * (1. / alpha_1_denom) +
        (0.5 + eta_star / 15.6) * 25. * y2 / alpha_1_denom;

    // alpha_2, Eqn. (53)
    const double alpha_2 =
        (0.63 + 0.04 * pow(eta_star, 1.45)) / (1. + 0.02 * pow(eta_star, 2.5));

    const double pow_eta_star_1_5 = pow(eta_star, 1.5);

    // alpha_3, Eqn. (53)
    const double alpha_3 =
        1.2 * SafeExp(0.6 * eta_star - 0.4 * pow_eta_star_1_5);

    // p_1, Eqn. (53)
    const double p_1 = (1.8 + 0.45 * eta_star) / (1. + 0.15 * pow_eta_star_1_5);

    // p_2, Eqn. (53)
    const double p_2 = 2.3 - 0.05 * eta_star / (1. + 0.025 * eta_star);

    // g, Eqn. (52)
    const double g =
        (alpha_1 + alpha_2 * pow(y, p_1)) /
        (1. + alpha_3 * pow(y, p_2) + alpha_2 * pow(y, p_1 + 2.) / 13.75);

    assert(g >= 0.);

    return g;
}

/* Compute the absorption kernels for a given NN Bremsstrahlung channel */
double BremSingleChannelAbsKernel(const double n_nuc, const double m_nuc,
                                  BremKernelParams* kernel_params,
                                  MyEOSParams* eos_params)
{
    static const double kHSquared = kH * kH; // [MeV^2 s^2]

    // EOS parameters
    const double temp = eos_params->temp; // temperature [MeV]

    // kernel parameters
    const double omega = kernel_params->omega; // neutrino energy [MeV]
    const double omega_prime =
        kernel_params->omega_prime; // primed neutrino energy [MeV]

    // dimensionless neutrino energy sum
    const double x = (omega + omega_prime) / temp;

    // temperature in units of 10 MeV
    const double temp_10 = temp / 10.;

    // nucleon effective degeneracy parameter, Eqn. (36) using baryon number
    // density instead of matter density
    const double eta_star = pow(3. * kPiSquared * n_nuc, 2. / 3.) *
                            (kHSquared * kMeV / kFourPiSquared) /
                            (2. * m_nuc * temp);

    // spin-fluctuation rate gamma, Eqn. (37)
    const double gamma = 1.63 * pow(eta_star, 3. / 2.) * temp_10;

    // pion mass parameter y, Eqn. (38)
    const double y = 1.94 / temp_10;

    // dimensionless fitting parameter s
    const double sb = BremKernelS(x, y, eta_star);

    // dimensionless fitting parameter g
    const double gb = BremKernelG(y, eta_star);

    // differential absorption kernel, Eqn. (35)
    const double s_kernel_abs =
        gamma / (x * x + pow(0.5 * gamma * gb, 2.)) * sb / temp;

    return s_kernel_abs;
}

/* Compute the angular independent part of the absorption kernels for the
 Bremsstrahlung reactions by summing the contributions of all NN channels */
double BremAllChannelsAbsKernel(BremKernelParams* kernel_params,
                                MyEOSParams* eos_params)
{
    static const double kTwentyeightThirds = 28. / 3.;

    static const double kMnGrams =
        kMn * kMeV / (kClight * kClight); // neutron mass in grams
    static const double kMpGrams =
        kMp * kMeV / (kClight * kClight); // proton mass in grams
    const double kMAvgGrams =
        sqrt(kMnGrams * kMpGrams); // geometric mean of nucleon masses in grams

    static const double kBremConst =
        kCa * kCa * kGf * kGf / kHbar; // kClight * pow(kHbar * kClight, 5.) *
                                       // kCa * kCa * kGfBrem * kGfBrem

    // EOS parameters
    const double nb = eos_params->nb; // baryon number density [cm^-3]
    const double xn = eos_params->yn; // neutron abundance/mass fraction
    const double xp = eos_params->yp; // proton abundance/mass fraction

    const double x_mean =
        sqrt(xn * xp); // geometric mean of nucleon abundances/mass fractions

    const double nn = nb * xn; // neutron number density [cm^-3]
    const double np = nb * xp; // protron number density [cm^-3]
    const double n_mean =
        nb * x_mean; // geometric mean of nucleon number densities [cm^-3]

    // compute single channel kernels
    const double s_abs_nn = BremSingleChannelAbsKernel(
        nn, kMnGrams, kernel_params, eos_params); // neutron-neutron
    const double s_abs_pp = BremSingleChannelAbsKernel(
        np, kMpGrams, kernel_params, eos_params); // const proton-proton
    const double s_abs_np = BremSingleChannelAbsKernel(
        n_mean, kMAvgGrams, kernel_params, eos_params); // neutron-proton

    // total absorption kernel
    double s_abs_tot = kBremConst * (nn * s_abs_nn + np * s_abs_pp +
                                     kTwentyeightThirds * n_mean * s_abs_np);

    // kernel correction due to medium dependence as in Fischer2016
    if (kernel_params->use_NN_medium_corr == true)
        s_abs_tot = s_abs_tot /
                    pow((1. + pow(nb / 0.15E+39, 0.3333333333333333) / 3.), 6.);

    return s_abs_tot;
}

/* Compute a specific Legendre coefficient in the expansion of production and
 * absorption kernels for the Bremsstrahlung reactions */
MyKernelOutput BremKernelsLegCoeff(BremKernelParams* kernel_params,
                                   MyEOSParams* eos_params)
{
    // kernel parameters
    const int l        = kernel_params->l;     // order of Legendre coefficient
    const double omega = kernel_params->omega; // neutrino energy [MeV]
    const double omega_prime =
        kernel_params->omega_prime; // primed neutrino energy [MeV]

    assert(l >= 0 && l <= 1);

    // EOS parameters
    const double temp = eos_params->temp; // temperature [MeV]

    // dimensionless neutrino energy sum
    const double x = (omega + omega_prime) / temp;

    // angular independent part of absorption kernel
    double s_abs = BremAllChannelsAbsKernel(kernel_params, eos_params);

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
        exit(EXIT_FAILURE);
    }

    // production kernel from detailed balance
    double s_em = s_abs * SafeExp(-x);

    MyKernelOutput brem_kernel;

    for (int idx = 0; idx < total_num_species; idx++)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}


void BremKernelsTable(const int n, double* nu_array,
                      GreyOpacityParams* grey_pars, M1Matrix* out)
{
    MyKernelOutput brem_ker;

    grey_pars->kernel_pars.brem_kernel_params.l = 0;
    grey_pars->kernel_pars.brem_kernel_params.use_NN_medium_corr =
        grey_pars->opacity_pars.use_NN_medium_corr;

    for (int i = 0; i < n; i++)
    {

        for (int j = i; j < n; j++)
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
double QBrem_BRT06(const double nb, const double temp, const double xn,
                   const double xp)
{
    const double rho = nb * kMb; // mass density [g cm-3]
    return 2.0778E+02 * 0.5 * kMeV * (xn * xn + xp * xp + 28. / 3. * xn * xp) *
           rho * rho * pow(temp, 5.5); // [MeV cm-3 s-1]
}

// Bremsstrahlung kernel from BRT06 Eq.(143) rewritten consistently
// to fit within the framework of the present library
MyKernelOutput BremKernelsBRT06(BremKernelParams* kernel_params,
                                MyEOSParams* eos_pars)
{
    static const double kHClight6FourPiSquared =
        kHClight * kHClight * kHClight * kHClight * kHClight * kHClight /
        (16. * kPi * kPi);
    const double omega       = kernel_params->omega;
    const double omega_prime = kernel_params->omega_prime;
    const double temp        = eos_pars->temp;

    const double x = 0.5 * (omega + omega_prime) / temp;
    const double q_nb =
        QBrem_BRT06(eos_pars->nb, temp, eos_pars->yn, eos_pars->yp);

    const double tmp = kHClight6FourPiSquared * kC4BRT06 *
                       (q_nb / pow(temp, 7.)) * bessk1(x) / x;
    const double s_em  = tmp * SafeExp(-x);
    const double s_abs = tmp * SafeExp(x);

    MyKernelOutput brem_kernel;
    for (int idx = 0; idx < total_num_species; idx++)
    {
        brem_kernel.abs[idx] = s_abs;
        brem_kernel.em[idx]  = s_em;
    }

    return brem_kernel;
}

void BremKernelsTableBRT06(const int n, double* nu_array,
                           GreyOpacityParams* grey_pars, M1Matrix* out)
{
    MyKernelOutput brem_ker;

    for (int i = 0; i < n; i++)
    {

        for (int j = i; j < n; j++)
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