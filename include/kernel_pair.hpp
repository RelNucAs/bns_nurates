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

#ifdef GSL_INCLUDES_H_
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_zeta.h>
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
        exit(EXIT_FAILURE);
    }
}

KOKKOS_INLINE_FUNCTION
double PairTFittedOld(int l, double alpha)
{
    double fit_poly = 1.;
    double fit_exp =
        exp(1.1599894050697828e-11 -
            1.0000000000000557 *
                alpha); // fitting coefficients are the same for l = 1,...,6

    if (alpha < 15.)
    {
        if (alpha == 0)
            return PairTWithAlpha0(l);

        switch (l)
        {

        case 1:
            fit_poly =
                0.6932349884688382 +
                alpha *
                    (0.19153844068873888 +
                     alpha *
                         (-0.022454086684801058 +
                          alpha *
                              (-0.01862035702348546 +
                               alpha *
                                   (0.010059951727924867 +
                                    alpha *
                                        (-0.0025591672605759157 +
                                         alpha *
                                             (0.00040809016119114534 +
                                              alpha *
                                                  (-4.382673963153738e-05 +
                                                   alpha *
                                                       (3.2238390304643073e-06 +
                                                        alpha *
                                                            (-1.6018407662684595e-07 +
                                                             alpha *
                                                                 (5.130720144162727e-09 +
                                                                  alpha *
                                                                      (-9.531900556136287e-11 +
                                                                       alpha *
                                                                           7.764460276911894e-13)))))))))));
            break;

        case 2:
            fit_poly =
                0.8224794727922443 +
                alpha *
                    (0.1290371349050982 +
                     alpha *
                         (-0.030657101372762428 +
                          alpha *
                              (-0.003453495292249957 +
                               alpha *
                                   (0.00426086049575675 +
                                    alpha *
                                        (-0.0013457908088575755 +
                                         alpha *
                                             (0.0002501877313487743 +
                                              alpha *
                                                  (-3.0919945500091964e-05 +
                                                   alpha *
                                                       (2.6176160114632282e-06 +
                                                        alpha *
                                                            (-1.5051104090763755e-07 +
                                                             alpha *
                                                                 (5.628852629867178e-09 +
                                                                  alpha *
                                                                      (-1.2360258213665597e-10 +
                                                                       alpha *
                                                                           1.2097121671548247e-12)))))))))));
            break;

        case 3:
            fit_poly =
                0.9015387115221923 +
                alpha *
                    (0.0791164948210426 +
                     alpha *
                         (-0.025154468942278033 +
                          alpha *
                              (0.0021057203236187075 +
                               alpha *
                                   (0.0011880768678479855 +
                                    alpha *
                                        (-0.0005242569477484338 +
                                         alpha *
                                             (0.00011078575689764255 +
                                              alpha *
                                                  (-1.48222315936738e-05 +
                                                   alpha *
                                                       (1.331164063502168e-06 +
                                                        alpha *
                                                            (-8.029349298289137e-08 +
                                                             alpha *
                                                                 (3.1269450241639527e-09 +
                                                                  alpha *
                                                                      (-7.111741520288681e-11 +
                                                                       alpha *
                                                                           7.17877896667645e-13)))))))))));
            break;

        case 4:
            fit_poly =
                0.9470282577605772 +
                alpha *
                    (0.045559097999768754 +
                     alpha *
                         (-0.016997198961210103 +
                          alpha *
                              (0.0030071993630479362 +
                               alpha *
                                   (2.1281035907218678e-05 +
                                    alpha *
                                        (-0.00014559133060632152 +
                                         alpha *
                                             (3.856825962408047e-05 +
                                              alpha *
                                                  (-5.677096233745809e-06 +
                                                   alpha *
                                                       (5.389942800937989e-07 +
                                                        alpha *
                                                            (-3.3760229096333e-08 +
                                                             alpha *
                                                                 (1.3519272756431448e-09 +
                                                                  alpha *
                                                                      (-3.1423838815286814e-11 +
                                                                       alpha *
                                                                           3.228347864267874e-13)))))))))));
            break;

        case 5:
            fit_poly =
                0.9721171063352136 +
                alpha *
                    (0.025128000871717175 +
                     alpha *
                         (-0.010335683014480727 +
                          alpha *
                              (0.0023768431000457074 +
                               alpha *
                                   (-0.0002672125800539407 +
                                    alpha *
                                        (-1.096586220062312e-05 +
                                         alpha *
                                             (9.471179226979555e-06 +
                                              alpha *
                                                  (-1.7228094292976135e-06 +
                                                   alpha *
                                                       (1.7941662304610306e-07 +
                                                        alpha *
                                                            (-1.1844514209510994e-08 +
                                                             alpha *
                                                                 (4.908468308315054e-10 +
                                                                  alpha *
                                                                      (-1.1689674516055724e-11 +
                                                                       alpha *
                                                                           1.223042342844013e-13)))))))))));
            break;

        case 6:
            fit_poly =
                0.9855501299373773 +
                alpha *
                    (0.013449009951258391 +
                     alpha *
                         (-0.005890098985143524 +
                          alpha *
                              (0.00154830636662188 +
                               alpha *
                                   (-0.00025382574940076774 +
                                    alpha *
                                        (2.257168344772742e-05 +
                                         alpha *
                                             (7.929667164380111e-08 +
                                              alpha *
                                                  (-3.0974403462775056e-07 +
                                                   alpha *
                                                       (4.3475802836454606e-08 +
                                                        alpha *
                                                            (-3.249472337669054e-09 +
                                                             alpha *
                                                                 (1.4415672453643435e-10 +
                                                                  alpha *
                                                                      (-3.584054805259039e-12 +
                                                                       alpha *
                                                                           3.862549379641218e-14)))))))))));
            break;

        default:
            printf(
                "PairTFitted (kernel_pair.c): l = %d must be within 1 and 6\n",
                l);
            exit(EXIT_FAILURE);
        }
    }

    return fit_poly * fit_exp;
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
        sum += PairTWithAlpha0(2) / tgamma(k); // l=0
        for (int l = 1; l <= (int)((k - 1.) / 2.); l++)
        {
            tmp = tmp / (eta * eta);
            sum += PairTWithAlpha0(2 * l + 2) * tmp / tgamma(k - 1 - 2 * l + 1);
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

/* Calculate G_n(a,b) from Eqn. (12) of Pons et. al. (1998)
 *
 * Inputs:
 *    n:            integer
 *    a,b,eta,y,z:  double
 *
 * Output:
 *    G_n(a,b,eta,y,z) = F_n(eta,b) - F_n(eta,a) - F_n(eta+y+z,b) +
 * F_n(eta+y+z,a)
 * */

/* void PairGOptimized(int n, double eta, double y, double z, double* g_out) */
/* { */
/*     /\* */
/*     const double pair_f_y_z_1 = PairFBackup(n, eta, y + z); */
/*     const double pair_f_y_z_2 = PairFBackup(n, eta + y + z, y + z); */

/*     const double pair_f_y_1 = PairFBackup(n, eta, y); */
/*     const double pair_f_y_2 = PairFBackup(n, eta + y + z, y); */

/*     const double pair_f_z_1 = PairFBackup(n, eta, z); */
/*     const double pair_f_z_2 = PairFBackup(n, eta + y + z, z); */
/*     *\/ */


/*     const double eta_y_z = eta + y + z; */

/*     const double pair_f_y_z_1 = PairFOptimized(n, eta, y + z); */
/*     const double pair_f_y_z_2 = PairFOptimized(n, eta_y_z, y + z); */

/*     const double pair_f_y_1 = PairFOptimized(n, eta, y); */
/*     const double pair_f_y_2 = PairFOptimized(n, eta_y_z, y); */

/*     const double pair_f_z_1 = PairFOptimized(n, eta, z); */
/*     const double pair_f_z_2 = PairFOptimized(n, eta_y_z, z); */


/*     double pair_f_min_1 = (y > z) ? pair_f_z_1 : pair_f_y_1; */
/*     double pair_f_max_1 = (y > z) ? pair_f_y_1 : pair_f_z_1; */

/*     double pair_f_min_2 = (y > z) ? pair_f_z_2 : pair_f_y_2; */
/*     double pair_f_max_2 = (y > z) ? pair_f_y_2 : pair_f_z_2; */

/*     g_out[0] = pair_f_y_z_1 - pair_f_y_1 - (pair_f_y_z_2 - pair_f_y_2); */
/*     g_out[1] = pair_f_y_z_1 - pair_f_z_1 - (pair_f_y_z_2 - pair_f_z_2); */
/*     g_out[2] = pair_f_min_1 - pair_f_min_2; */
/*     g_out[3] = pair_f_y_z_1 - pair_f_max_1 - (pair_f_y_z_2 - pair_f_max_2);
 */

/*     return; */
/* } */

/* void F_for_G(double eta, double y, double z, double* F_out) */
/* { */
/*     F_out[0 + 0 * 8] = FDI_0(eta); */
/*     F_out[1 + 0 * 8] = FDI_0(eta - y); */
/*     F_out[2 + 0 * 8] = FDI_0(eta - z); */
/*     F_out[3 + 0 * 8] = FDI_0(eta + y); */
/*     F_out[4 + 0 * 8] = FDI_0(eta + z); */
/*     F_out[5 + 0 * 8] = FDI_0(eta - y - z); */
/*     F_out[6 + 0 * 8] = (y <= z) ? F_out[1] : F_out[2]; */
/*     F_out[7 + 0 * 8] = (y <= z) ? F_out[4] : F_out[3]; */

/*     F_out[0 + 1 * 8] = FDI_p1(eta); */
/*     F_out[1 + 1 * 8] = FDI_p1(eta - y); */
/*     F_out[2 + 1 * 8] = FDI_p1(eta - z); */
/*     F_out[3 + 1 * 8] = FDI_p1(eta + y); */
/*     F_out[4 + 1 * 8] = FDI_p1(eta + z); */
/*     F_out[5 + 1 * 8] = FDI_p1(eta - y - z); */
/*     F_out[6 + 1 * 8] = (y <= z) ? F_out[1 + 1 * 8] : F_out[2 + 1 * 8]; */
/*     F_out[7 + 1 * 8] = (y <= z) ? F_out[4 + 1 * 8] : F_out[3 + 1 * 8]; */

/*     F_out[0 + 2 * 8] = FDI_p2(eta); */
/*     F_out[1 + 2 * 8] = FDI_p2(eta - y); */
/*     F_out[2 + 2 * 8] = FDI_p2(eta - z); */
/*     F_out[3 + 2 * 8] = FDI_p2(eta + y); */
/*     F_out[4 + 2 * 8] = FDI_p2(eta + z); */
/*     F_out[5 + 2 * 8] = FDI_p2(eta - y - z); */
/*     F_out[6 + 2 * 8] = (y <= z) ? F_out[1 + 2 * 8] : F_out[2 + 2 * 8]; */
/*     F_out[7 + 2 * 8] = (y <= z) ? F_out[4 + 2 * 8] : F_out[3 + 2 * 8]; */

/*     F_out[0 + 3 * 8] = FDI_p3(eta); */
/*     F_out[1 + 3 * 8] = FDI_p3(eta - y); */
/*     F_out[2 + 3 * 8] = FDI_p3(eta - z); */
/*     F_out[3 + 3 * 8] = FDI_p3(eta + y); */
/*     F_out[4 + 3 * 8] = FDI_p3(eta + z); */
/*     F_out[5 + 3 * 8] = FDI_p3(eta - y - z); */
/*     F_out[6 + 3 * 8] = (y <= z) ? F_out[1 + 3 * 8] : F_out[2 + 3 * 8]; */
/*     F_out[7 + 3 * 8] = (y <= z) ? F_out[4 + 3 * 8] : F_out[3 + 3 * 8]; */

/*     F_out[0 + 4 * 8] = FDI_p4(eta); */
/*     F_out[1 + 4 * 8] = FDI_p4(eta - y); */
/*     F_out[2 + 4 * 8] = FDI_p4(eta - z); */
/*     F_out[3 + 4 * 8] = FDI_p4(eta + y); */
/*     F_out[4 + 4 * 8] = FDI_p4(eta + z); */
/*     F_out[5 + 4 * 8] = FDI_p4(eta - y - z); */
/*     F_out[6 + 4 * 8] = (y <= z) ? F_out[1 + 4 * 8] : F_out[2 + 4 * 8]; */
/*     F_out[7 + 4 * 8] = (y <= z) ? F_out[4 + 4 * 8] : F_out[3 + 4 * 8]; */

/*     F_out[0 + 5 * 8] = FDI_p5(eta); */
/*     F_out[1 + 5 * 8] = FDI_p5(eta - y); */
/*     F_out[2 + 5 * 8] = FDI_p5(eta - z); */
/*     F_out[3 + 5 * 8] = FDI_p5(eta + y); */
/*     F_out[4 + 5 * 8] = FDI_p5(eta + z); */
/*     F_out[5 + 5 * 8] = FDI_p5(eta - y - z); */
/*     F_out[6 + 5 * 8] = (y <= z) ? F_out[1 + 5 * 8] : F_out[2 + 5 * 8]; */
/*     F_out[7 + 5 * 8] = (y <= z) ? F_out[4 + 5 * 8] : F_out[3 + 5 * 8]; */

/*     return; */
/* } */

/* void PairGOptimized(double eta, double y, double z, double* g_out) */
/* { */
/*     double m = MIN(y, z); */

/*     double F[48]; */

/*     F_for_G(eta, y, z, F); */

/*     g_out[0] = (1 * (POW0(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                      POW0(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))); */
/*     g_out[1] = (1 * (POW0(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                      POW0(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))); */
/*     g_out[2] = (y <= z) ? g_out[1] : g_out[0]; */
/*     g_out[3] = (-1 * (POW0(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))); */

/*     g_out[4] = (1 * (POW1(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                      POW1(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                (1 * (POW0(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) - */
/*                      POW0(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))); */
/*     g_out[5] = (1 * (POW1(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                      POW1(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                (1 * (POW0(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) - */
/*                      POW0(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))); */
/*     g_out[6] = (y <= z) ? g_out[5] : g_out[4]; */
/*     g_out[7] = (-1 * (POW1(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) + */
/*                (-1 * (POW0(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))); */

/*     g_out[8] = (1 * (POW2(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                      POW2(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                (2 * (POW1(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) - */
/*                      POW1(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                (1 * (POW0(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) - */
/*                      POW0(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))); */
/*     g_out[9] = (1 * (POW2(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                      POW2(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                (2 * (POW1(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) - */
/*                      POW1(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                (1 * (POW0(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) - */
/*                      POW0(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))); */
/*     g_out[10] = (y <= z) ? g_out[9] : g_out[8]; */
/*     g_out[11] = (-1 * (POW2(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) + */
/*                 (-2 * (POW1(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) + */
/*                 (-1 * (POW0(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))); */

/*     g_out[12] = (1 * (POW3(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                       POW3(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (3 * (POW2(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) - */
/*                       POW2(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (3 * (POW1(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) - */
/*                       POW1(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (1 * (POW0(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) - */
/*                       POW0(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))); */
/*     g_out[13] = (1 * (POW3(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                       POW3(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (3 * (POW2(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) - */
/*                       POW2(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (3 * (POW1(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) - */
/*                       POW1(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (1 * (POW0(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) - */
/*                       POW0(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))); */
/*     g_out[14] = (y <= z) ? g_out[13] : g_out[12]; */
/*     g_out[15] = (-1 * (POW3(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) + */
/*                 (-3 * (POW2(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) + */
/*                 (-3 * (POW1(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) + */
/*                 (-1 * (POW0(m) * (F[6 + 3 * 8] - F[7 + 3 * 8]))); */

/*     g_out[16] = (1 * (POW4(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                       POW4(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (4 * (POW3(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) - */
/*                       POW3(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (6 * (POW2(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) - */
/*                       POW2(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (4 * (POW1(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) - */
/*                       POW1(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) + */
/*                 (1 * (POW0(y) * (F[1 + 4 * 8] - F[4 + 4 * 8]) - */
/*                       POW0(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))); */
/*     g_out[17] = (1 * (POW4(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                       POW4(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (4 * (POW3(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) - */
/*                       POW3(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (6 * (POW2(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) - */
/*                       POW2(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (4 * (POW1(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) - */
/*                       POW1(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) + */
/*                 (1 * (POW0(z) * (F[2 + 4 * 8] - F[3 + 4 * 8]) - */
/*                       POW0(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))); */
/*     g_out[18] = (y <= z) ? g_out[17] : g_out[16]; */
/*     g_out[19] = (-1 * (POW4(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) + */
/*                 (-4 * (POW3(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) + */
/*                 (-6 * (POW2(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) + */
/*                 (-4 * (POW1(m) * (F[6 + 3 * 8] - F[7 + 3 * 8]))) + */
/*                 (-1 * (POW0(m) * (F[6 + 4 * 8] - F[7 + 4 * 8]))); */

/*     g_out[20] = (1 * (POW5(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) - */
/*                       POW5(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (5 * (POW4(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) - */
/*                       POW4(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (10 * (POW3(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) - */
/*                        POW3(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (10 * (POW2(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) - */
/*                        POW2(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) + */
/*                 (5 * (POW1(y) * (F[1 + 4 * 8] - F[4 + 4 * 8]) - */
/*                       POW1(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))) + */
/*                 (1 * (POW0(y) * (F[1 + 5 * 8] - F[4 + 5 * 8]) - */
/*                       POW0(y + z) * (F[5 + 5 * 8] - F[0 + 5 * 8]))); */
/*     g_out[21] = (1 * (POW5(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) - */
/*                       POW5(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) + */
/*                 (5 * (POW4(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) - */
/*                       POW4(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) + */
/*                 (10 * (POW3(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) - */
/*                        POW3(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) + */
/*                 (10 * (POW2(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) - */
/*                        POW2(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) + */
/*                 (5 * (POW1(z) * (F[2 + 4 * 8] - F[3 + 4 * 8]) - */
/*                       POW1(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))) + */
/*                 (1 * (POW0(z) * (F[2 + 5 * 8] - F[3 + 5 * 8]) - */
/*                       POW0(y + z) * (F[5 + 5 * 8] - F[0 + 5 * 8]))); */
/*     g_out[22] = (y <= z) ? g_out[21] : g_out[20]; */
/*     g_out[23] = (-1 * (POW5(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) + */
/*                 (-5 * (POW4(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) + */
/*                 (-10 * (POW3(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) + */
/*                 (-10 * (POW2(m) * (F[6 + 3 * 8] - F[7 + 3 * 8]))) + */
/*                 (-5 * (POW1(m) * (F[6 + 4 * 8] - F[7 + 4 * 8]))) + */
/*                 (-1 * (POW0(m) * (F[6 + 5 * 8] - F[7 + 5 * 8]))); */

/*     return; */
/* } */

/* Calculate Psi_l(y,z) from Eqn. (11) of Pons et. al.
 *
 * Inputs:
 *      l:    integer (mode number)
 *      y:    dimensionless energy
 *      z:    dimensionless energy
 *      eta:  electron degenracy parameter
 *
 * Output:
 *      Psi_l(y,z) = sum_{n=0..2} [c_{ln} G_n(y,y+z) + d_{ln} G_n(z,y+z)] +
 * sum_{n=3..(2l+5)} a_{ln} [G_n(0,min(y,z)) - G_n(max(y,z),y+z)]
 */
// @TODO: compute only a coeffs needed for the specific l value
KOKKOS_INLINE_FUNCTION
void PairPsiCoeffs(double y, double z, double a[4][12], double c[4][3],
                   double d[4][3])
{
    const double y2 = y * y;
    const double y3 = y * y2;
    const double y4 = y * y3;
    const double y5 = y * y4;

    const double z2 = z * z;
    const double z3 = z * z2;
    const double z4 = z * z3;
    const double z5 = z * z4;

    // @TODO: compute ratio between constant numbers (e.g 4/3) only once

    a[0][3] = 8. / (3. * y2);
    a[0][4] = -4. / (3. * y2 * z);
    a[0][5] = 4. / (15. * y2 * z2);
    c[0][0] = (4. * y / z2) * (2. * z2 / 3. + y * z + 2. * y2 / 5.);
    c[0][1] = -(4. * y / (3. * z2)) * (3. * y + 4. * z);
    c[0][2] = 8. * y / (3. * z2);
    d[0][0] = 4. * z3 / (15. * y2);
    d[0][1] = -4. * z2 / (3. * y2);
    d[0][2] = 8. * z / (3. * y2);

    // Coefficients for l = 1
    a[1][3] = 8. / (3. * y2);
    a[1][4] = -4. * (4. * y + 3. * z) / (3. * y3 * z);
    a[1][5] = 4. * (13. * y + 18. * z) / (15. * y3 * z2);
    a[1][6] = -4. * (y + 3. * z) / (5. * y3 * z3);
    a[1][7] = 16. / (35. * y3 * z3);
    c[1][0] = -(4. * y / z3) *
              (2. * y3 / 7. + 4. * y2 * z / 5. + 4. * y * z2 / 5. + z3 / 3.);
    c[1][1] = (4. * y / z3) * (4. * y2 / 5. + 7. * y * z / 5. + 2. * z2 / 3.);
    c[1][2] = -(4. * y / z3) * (3. * y / 5. + z / 3.);
    d[1][0] = -(4. * z3 / (105. * y3)) * (14. * y + 9. * z);
    d[1][1] = (4. * z2 / (5. * y3)) * (7. * y / 3. + 2. * z);
    d[1][2] = -(4. * z / y3) * (y / 3. + 3. * z / 5.);

    // Coefficients for l = 2
    a[2][3] = 8. / (3. * y2);
    a[2][4] = -4. * (10. * y + 9. * z) / (3. * y3 * z);
    a[2][5] = 4. * (73. * y2 + 126. * y * z + 36. * z * z) / (15. * y4 * z * z);
    a[2][6] = -12. * (y2 + 3. * y * z + 8. * z2 / 5.) / (y4 * z3);
    a[2][7] = 48. * (2. * y2 + 13. * y * z + 12. * z2) / (35. * y4 * z4);
    a[2][8] = -24. * (y + 2. * z) / (7. * y4 * z4);
    a[2][9] = 8. / (7. * y4 * z4);
    c[2][0] = 4. * y *
              (2. * y4 / 7. + 6. * y3 * z / 7 + 32. * y2 * z2 / 35. +
               2. * y * z3 / 5. + z4 / 15.) /
              z4;
    c[2][1] = -4. * y *
              (6. * y3 / 7. + 12. * y2 * z / 7. + y * z2 + 2. * z3 / 15.) / z4;
    c[2][2] = 8. * y * (z2 / 6. + 3 * y * z / 2. + 12 * y2 / 7.) / (5. * z4);
    d[2][0] = 4. * z3 * (16. * y2 + 27. * y * z + 12. * z2) / (105. * y4);
    d[2][1] = -4. * z2 * (18. * z2 / 35. + 6. * y * z / 7. + y2 / 3.) / y4;
    d[2][2] = 8. * z * (y2 / 6. + 3. * y * z / 2. + 12. * z2 / 7.) / (5. * y4);

    // Coefficients for l = 3
    a[3][3] = 8. / (3. * y2);
    a[3][4] = -4. * (19. * y + 18. * z) / (3. * y3 * z);
    a[3][5] = 4. * (253. * y2 + 468. * y * z + 180. * z2) / (15. * y4 * z2);
    a[3][6] = -8. * (149. * y3 + 447. * y2 * z + 330. * y * z2 + 50. * z3) /
              (15. * y5 * z3);
    a[3][7] = 8. *
              (116. * y3 + 2916. * y2 * z / 5. + 696. * y * z2 + 200. * z3) /
              (21. * y5 * z4);
    a[3][8] = -40. * (5. * y3 + 54. * y2 * z + 108. * y * z2 + 50. * z3) /
              (21. * y5 * z5);
    a[3][9] = 40. * (10. * y2 + 43. * y * z + 100. * z2 / 3.) / (21. * y5 * z5);
    a[3][10] = -40. * (3. * y + 5. * z) / (9. * y5 * z5);
    a[3][11] = 320. / (99. * y5 * z5);
    c[3][0]  = -4. * y2 *
              (10. * y4 / 33. + 20. * y3 * z / 21. + 68. * y2 * z2 / 63. +
               18. * y * z3 / 35. + 3. * z4 / 35.) /
              z5;
    c[3][1] = 4. * y2 *
              (20. * y3 / 21. + 130. * y2 * z / 63. + 48. * y * z2 / 35. +
               9. * z3 / 35.) /
              z5;
    c[3][2] =
        -4. * y2 * (50. * y2 / 63. + 6. * y * z / 7. + 6. * z2 / 35.) / z5;
    d[3][0] = -4. * z3 *
              (9. * y3 + 34. * y2 * z + 40. * y * z2 + 500. * z3 / 33.) /
              (105. * y5);
    d[3][1] = 4. * z2 *
              (3. * y3 + 24. * y2 * z + 130. * y * z2 / 3. + 200. * z3 / 9.) /
              (35. * y5);
    d[3][2] =
        -4. * z2 * (50. * z2 / 63. + 6. * z * y / 7. + 6. * y2 / 35.) / y5;

    return;
}

KOKKOS_INLINE_FUNCTION
void PairPsiOptimized(int l, double y, double z, double eta, double* psi_out)
{
    assert(l == 0);

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
        (POW2(kAlpha1[0]) * pair_psi[0] + POW2(kAlpha2[0]) * pair_psi[1]) /
        phi_denom;
    phi_out[1] =
        phi_prefactor *
        (POW2(kAlpha1[0]) * pair_psi[1] + POW2(kAlpha2[0]) * pair_psi[0]) /
        phi_denom;
    phi_out[2] =
        phi_prefactor *
        (POW2(kAlpha1[1]) * pair_psi[0] + POW2(kAlpha2[1]) * pair_psi[1]) /
        phi_denom;
    phi_out[3] =
        phi_prefactor *
        (POW2(kAlpha1[1]) * pair_psi[1] + POW2(kAlpha2[1]) * pair_psi[0]) /
        phi_denom;

    return;
}

KOKKOS_INLINE_FUNCTION
MyKernelOutput PairKernelsOptimized(MyEOSParams* eos_pars,
                                    PairKernelParams* kernel_pars)
{
    // EOS specific parameters
    double eta  = eos_pars->mu_e / eos_pars->temp;
    double temp = eos_pars->temp;

    // kernel specific parameters
    double omega       = kernel_pars->omega;
    double omega_prime = kernel_pars->omega_prime;

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

//=====================================
// Tested, but not optimized functions
//=====================================

#ifdef GSL_INCLUDES_H_
KOKKOS_INLINE_FUNCTION
double PairT(int l, double alpha, double tolerance)
{

    assert(alpha >= 0 && l >= 1);

    if (alpha == 0 && l != 1)
    {
        return pow(2., -l) * (pow(2., l) - 2.) *
               gsl_sf_zeta_int(l); // Computed in Mathematica
    }
    else if (alpha == 0 && l == 1)
    {
        return log(2.0);
    }
    else
    {
        double val    = fabs(42. * tolerance);
        double result = 0.;
        int n         = 1;
        while (tolerance < fabs(val))
        {
            val = pow(-1., n + 1) * exp(-n * alpha) / pow(n, l);
            result += val;
            n++;
        }

        // printf("l = %d, alpha = %.3e: n = %d\n", l, alpha, n-1);
        return result;
    }
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
double PairFBackup(int k, double eta, double x1)
{

    double result = 0.;
    double tol    = 1.0E-10;

    if (eta < 0.)
    {
        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum += PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
        }
        result = tgamma(k + 1) * (PairT(k + 1., -eta, tol) - sum);
    }
    else if (k != 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        for (int l = 0; l <= (int)((k - 1.) / 2.); l++)
        {
            sum += PairTWithAlpha0(2 * l + 2) * pow(eta, k - 1. - 2. * l) /
                   tgamma(k - 1 - 2 * l + 1);
        }

        double sum_2 = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum_2 +=
                PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 tgamma(k + 1) *
                     (2.0 * sum + pow(-1., k) * PairT(k + 1, eta, tol) - sum_2);
    }
    else if (k == 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum += PairT(k + 1 - l, x1 - eta, tol) * pow(x1, l) / tgamma(l + 1);
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * PairT(k + 1, eta, tol) - sum);
    }
    else
    {

        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum += pow(-1., k - l) * PairT(k + 1 - l, eta - x1, tol) *
                   pow(x1, l) / tgamma(l + 1);
        }
        result = pow(x1, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * PairT(k + 1, eta, tol) - sum);
    }

    return result;
}

KOKKOS_INLINE_FUNCTION
double PairF(int k, double eta, double x1)
{

    double result = 0.;
    // double tol = 1e-6;

    if (eta < 0.)
    {
        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum +=
                PairTFitted(k + 1 - l, x1 - eta) * pow(x1, l) / tgamma(l + 1);
        }
        result = tgamma(k + 1) * (PairTFitted(k + 1., -eta) - sum);
    }
    else if (k != 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        for (int l = 0; l <= (int)((k - 1.) / 2.); l++)
        {
            sum += PairTWithAlpha0(2 * l + 2) * pow(eta, k - 1. - 2. * l) /
                   tgamma(k - 1 - 2 * l + 1);
        }

        double sum_2 = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum_2 +=
                PairTFitted(k + 1 - l, x1 - eta) * pow(x1, l) / tgamma(l + 1);
        }

        result =
            pow(eta, k + 1.) / (k + 1.) +
            tgamma(k + 1) *
                (2.0 * sum + pow(-1., k) * PairTFitted(k + 1, eta) - sum_2);
    }
    else if (k == 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum +=
                PairTFitted(k + 1 - l, x1 - eta) * pow(x1, l) / tgamma(l + 1);
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * PairTFitted(k + 1, eta) - sum);
    }
    else
    {

        double sum = 0.;
        for (int l = 0; l <= k; l++)
        {
            sum += pow(-1., k - l) * PairTFitted(k + 1 - l, eta - x1) *
                   pow(x1, l) / tgamma(l + 1);
        }
        result = pow(x1, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * PairTFitted(k + 1, eta) - sum);
    }

    return result;
}


/* Calculate G_n(a,b) from Eqn. (12) of Pons et. al. (1998)
 *
 * Inputs:
 *    n:            integer
 *    a,b,eta,y,z:  double
 *
 * Output:
 *    G_n(a,b,eta,y,z) = F_n(eta,b) - F_n(eta,a) - F_n(eta+y+z,b) +
 * F_n(eta+y+z,a)
 * */
KOKKOS_INLINE_FUNCTION
double PairG(int n, double a, double b, double eta, double y, double z)
{

    double result = 0.;
    result = PairF(n, eta, b) - PairF(n, eta, a) - PairF(n, eta + y + z, b) +
             PairF(n, eta + y + z, a);

    return result;
}

/* Calculate Psi_l(y,z) from Eqn. (11) of Pons et. al.
 *
 * Inputs:
 *      l:    integer (mode number)
 *      y:    dimensionless energy
 *      z:    dimensionless energy
 *      eta:  electron degenracy parameter
 *
 * Output:
 *      Psi_l(y,z) = sum_{n=0..2} [c_{ln} G_n(y,y+z) + d_{ln} G_n(z,y+z)] +
 * sum_{n=3..(2l+5)} a_{ln} [G_n(0,min(y,z)) - G_n(max(y,z),y+z)]
 */
KOKKOS_INLINE_FUNCTION
double PairPsi(int l, double y, double z, double eta)
{

    assert(0 <= l && l <= 3);

    double a[4][12], c[4][3], d[4][3];

    a[0][3] = 8. / (3. * y * y);
    a[0][4] = -4. / (3. * y * y * z);
    a[0][5] = 4. / (15. * y * y * z * z);
    c[0][0] = (4. * y / (z * z)) * (2. * z * z / 3. + y * z + 2. * y * y / 5.);
    c[0][1] = -(4. * y / (3. * z * z)) * (3. * y + 4. * z);
    c[0][2] = 8. * y / (3. * z * z);
    d[0][0] = 4. * z * z * z / (15. * y * y);
    d[0][1] = -4. * z * z / (3. * y * y);
    d[0][2] = 8. * z / (3. * y * y);

    // Coefficients for l = 1
    a[1][3] = 8. / (3. * y * y);
    a[1][4] = -4. * (4. * y + 3. * z) / (3. * y * y * y * z);
    a[1][5] = 4. * (13. * y + 18. * z) / (15. * y * y * y * z * z);
    a[1][6] = -4. * (y + 3. * z) / (5. * y * y * y * z * z * z);
    a[1][7] = 16. / (35. * y * y * y * z * z * z);
    c[1][0] =
        -(4. * y / (z * z * z)) * (2. * y * y * y / 7. + 4. * y * y * z / 5. +
                                   4. * y * z * z / 5. + z * z * z / 3.);
    c[1][1] = (4. * y / (z * z * z)) *
              (4. * y * y / 5. + 7. * y * z / 5. + 2. * z * z / 3.);
    c[1][2] = -(4. * y / (z * z * z)) * (3. * y / 5. + z / 3.);
    d[1][0] = -(4. * z * z * z / (105. * y * y * y)) * (14. * y + 9. * z);
    d[1][1] = (4. * z * z / (5. * y * y * y)) * (7. * y / 3. + 2. * z);
    d[1][2] = -(4. * z / (y * y * y)) * (y / 3. + 3. * z / 5.);

    // Coefficients for l = 2
    a[2][3] = 8. / (3. * y * y);
    a[2][4] = -4. * (10. * y + 9. * z) / (3. * y * y * y * z);
    a[2][5] = 4. * (73. * y * y + 126. * y * z + 36. * z * z) /
              (15. * y * y * y * y * z * z);
    a[2][6] = -12. * (y * y + 3. * y * z + 8. * z * z / 5.) /
              (y * y * y * y * z * z * z);
    a[2][7] = 48. * (2. * y * y + 13. * y * z + 12. * z * z) /
              (35. * y * y * y * y * z * z * z * z);
    a[2][8] = -24. * (y + 2. * z) / (7. * y * y * y * y * z * z * z * z);
    a[2][9] = 8. / (7. * y * y * y * y * z * z * z * z);
    c[2][0] = 4. * y *
              (2. * y * y * y * y / 7. + 6. * y * y * y * z / 7 +
               32. * y * y * z * z / 35. + 2. * y * z * z * z / 5. +
               z * z * z * z / 15.) /
              (z * z * z * z);
    c[2][1] = -4. * y *
              (6. * y * y * y / 7. + 12. * y * y * z / 7. + y * z * z +
               2. * z * z * z / 15.) /
              (z * z * z * z);
    c[2][2] = 8. * y * (z * z / 6. + 3 * y * z / 2. + 12 * y * y / 7.) /
              (5. * z * z * z * z);
    d[2][0] = 4. * z * z * z * (16. * y * y + 27. * y * z + 12. * z * z) /
              (105. * y * y * y * y);
    d[2][1] = -4. * z * z * (18. * z * z / 35. + 6. * y * z / 7. + y * y / 3.) /
              (y * y * y * y);
    d[2][2] = 8. * z * (y * y / 6. + 3. * y * z / 2. + 12. * z * z / 7.) /
              (5. * y * y * y * y);

    // Coefficients for l = 3
    a[3][3] = 8. / (3. * y * y);
    a[3][4] = -4. * (19. * y + 18. * z) / (3. * y * y * y * z);
    a[3][5] = 4. * (253. * y * y + 468. * y * z + 180. * z * z) /
              (15. * y * y * y * y * z * z);
    a[3][6] = -8. *
              (149. * y * y * y + 447. * y * y * z + 330. * y * z * z +
               50. * z * z * z) /
              (15. * y * y * y * y * y * z * z * z);
    a[3][7] = 8. *
              (116. * y * y * y + 2916. * y * y * z / 5. + 696. * y * z * z +
               200. * z * z * z) /
              (21. * y * y * y * y * y * z * z * z * z);
    a[3][8] = -40. *
              (5. * y * y * y + 54. * y * y * z + 108. * y * z * z +
               50. * z * z * z) /
              (21. * y * y * y * y * y * z * z * z * z * z);
    a[3][9] = 40. * (10. * y * y + 43. * y * z + 100. * z * z / 3.) /
              (21. * y * y * y * y * y * z * z * z * z * z);
    a[3][10] =
        -40. * (3. * y + 5. * z) / (9. * y * y * y * y * y * z * z * z * z * z);
    a[3][11] = 320. / (99. * y * y * y * y * y * z * z * z * z * z);
    c[3][0]  = -4. * y * y *
              (10. * y * y * y * y / 33. + 20. * y * y * y * z / 21. +
               68. * y * y * z * z / 63. + 18. * y * z * z * z / 35. +
               3. * z * z * z * z / 35.) /
              (z * z * z * z * z);
    c[3][1] = 4. * y * y *
              (20. * y * y * y / 21. + 130. * y * y * z / 63. +
               48. * y * z * z / 35. + 9. * z * z * z / 35.) /
              (z * z * z * z * z);
    c[3][2] = -4. * y * y *
              (50. * y * y / 63. + 6. * y * z / 7. + 6. * z * z / 35.) /
              (z * z * z * z * z);
    d[3][0] = -4. * z * z * z *
              (9. * y * y * y + 34. * y * y * z + 40. * y * z * z +
               500. * z * z * z / 33.) /
              (105. * y * y * y * y * y);
    d[3][1] = 4. * z * z *
              (3. * y * y * y + 24. * y * y * z + 130. * y * z * z / 3. +
               200. * z * z * z / 9.) /
              (35. * y * y * y * y * y);
    d[3][2] = -4. * z * z *
              (50. * z * z / 63. + 6. * z * y / 7. + 6. * y * y / 35.) /
              (y * y * y * y * y);

    double result = 0.;

    for (int n = 0; n <= 2; n++)
    {
        result += c[l][n] * PairG(n, y, y + z, eta, y, z) +
                  d[l][n] * PairG(n, z, y + z, eta, y, z);
    }

    double min_yz = (y > z) ? z : y;
    double max_yz = (y > z) ? y : z;

    for (int n = 3; n <= 2 * l + 5; n++)
    {
        result += a[l][n] * (PairG(n, 0, min_yz, eta, y, z) -
                             PairG(n, max_yz, y + z, eta, y, z));
    }

    return result;
}

KOKKOS_INLINE_FUNCTION
double PairPhi(int l, double omega, double omega_prime, double eta, double temp,
               int e_x)
{

    assert(e_x >= 0 && e_x <= 1);

    double y = omega / temp;
    double z = omega_prime / temp;

    double psi_1 = PairPsi(l, y, z, eta);
    double psi_2 = PairPsi(l, z, y, eta);

    double result = kGSqr * temp * temp *
                    (kAlpha1[e_x] * kAlpha1[e_x] * psi_1 +
                     kAlpha2[e_x] * kAlpha2[e_x] * psi_2) /
                    (kPi * (1. - exp(y + z)));

    return result;
}

/* Calculates the production and absorption kernels for the pair process from
 * Eqns. (2) and (3) of Pons et. al. integrated over the phi variable
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega')
 * P_l(mu) P_l(mu_prime) 2 pi R^a_{TP}(omega,omega',cos theta) =
 * e^{(omega+omega')/T} R^p_{TP}
 *
 * l is an integer.
 */
KOKKOS_INLINE_FUNCTION
MyKernelQuantity PairKernelsPhiIntegrated(MyEOSParams* eos_pars,
                                          PairKernelParams* kernel_pars)
{

    // kernel specific parameters
    double omega       = kernel_pars->omega;
    double omega_prime = kernel_pars->omega_prime;
    double mu          = kernel_pars->mu;
    double mu_prime    = kernel_pars->mu_prime;
    int lmax           = kernel_pars->lmax;
    double filterpar   = kernel_pars->filter;

    // EOS specific parameters
    double eta  = eos_pars->mu_e / eos_pars->temp;
    double temp = eos_pars->temp;

    double pair_kernel_production_e = 0.;
    double pair_kernel_absorption_e = 0.;
    double pair_kernel_production_x = 0.;
    double pair_kernel_absorption_x = 0.;
    double pair_phi_e               = 0.;
    double pair_phi_x               = 0.;

    assert(lmax >= 0 && lmax <= 3);

    for (int l = 0; l <= lmax; l++)
    {
        double legendre_l_mu       = gsl_sf_legendre_Pl(l, mu);
        double legendre_l_mu_prime = gsl_sf_legendre_Pl(l, mu_prime);
        pair_phi_e = PairPhi(l, omega, omega_prime, eta, temp, 0);
        pair_phi_x = PairPhi(l, omega, omega_prime, eta, temp, 1);

        pair_kernel_production_e +=
            2. * kPi * (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) *
            (2. * l + 1.) * pair_phi_e * legendre_l_mu * legendre_l_mu_prime /
            2.;
        pair_kernel_production_x +=
            2. * kPi * (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) *
            (2. * l + 1.) * pair_phi_x * legendre_l_mu * legendre_l_mu_prime /
            2.;
    }

    pair_kernel_absorption_e =
        exp((omega + omega_prime) / temp) * pair_kernel_production_e;
    pair_kernel_absorption_x =
        exp((omega + omega_prime) / temp) * pair_kernel_production_x;

  MyKernelQuantity pair_kernel = {
      .em_e  = pair_kernel_production_e, .abs_e = pair_kernel_absorption_e,
      .em_x  = pair_kernel_production_x, .abs_x = pair_kernel_absorption_x};

    return pair_kernel;
}

KOKKOS_INLINE_FUNCTION
MyKernelQuantity PairKernelsPhiMuIntegrated(MyEOSParams* eos_pars,
                                            PairKernelParams* kernel_pars)
{

    // kernel specific parameters
    double omega       = kernel_pars->omega;
    double omega_prime = kernel_pars->omega_prime;

    // EOS specific parameters
    double eta  = eos_pars->mu_e / eos_pars->temp;
    double temp = eos_pars->temp;

    int l = 0;

    double pair_phi_e = PairPhi(l, omega, omega_prime, eta, temp, 0);
    double pair_phi_x = PairPhi(l, omega, omega_prime, eta, temp, 1);

    double pair_kernel_production_e = 2. * kPi * pair_phi_e;
    double pair_kernel_production_x = 2. * kPi * pair_phi_x;

    double pair_kernel_absorption_e =
        SafeExp((omega + omega_prime) / temp) * pair_kernel_production_e;
    double pair_kernel_absorption_x =
        SafeExp((omega + omega_prime) / temp) * pair_kernel_production_x;

    const double kPrefactor = 1. / (kClight * pow(kH * kClight, 3.));
    // const double kPrefactor = 1.;

    MyKernelQuantity pair_kernel = {
        .em_e  = kPrefactor * pair_kernel_production_e,
        .abs_e = kPrefactor * pair_kernel_absorption_e,
        .em_x  = kPrefactor * pair_kernel_production_x,
        .abs_x = kPrefactor * pair_kernel_absorption_x
    };

    return pair_kernel;
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

/* Calculates the production and absorption kernels for the pair process from
 * Eqns. (2) and (3) of Pons et. al.
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega')
 * P_l(cos theta) R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T}
 * R^p_{TP}(omega,omega',cos theta)
 *
 * l is an integer.
 */
KOKKOS_INLINE_FUNCTION
MyKernelQuantity PairKernels(MyEOSParams* eos_pars,
                             PairKernelParams* kernel_pars)
{

    // kernel specific parameters
    double omega       = kernel_pars->omega;
    double omega_prime = kernel_pars->omega_prime;
    double cos_theta   = kernel_pars->cos_theta;
    int lmax           = kernel_pars->lmax;
    double filterpar   = kernel_pars->filter;

    // EOS specific parameters
    double eta  = eos_pars->mu_e / eos_pars->temp;
    double temp = eos_pars->temp;

    double pair_kernel_production_e = 0.;
    double pair_kernel_absorption_e = 0.;
    double pair_kernel_production_x = 0.;
    double pair_kernel_absorption_x = 0.;
    double pair_phi_e               = 0.;
    double pair_phi_x               = 0.;

    assert(lmax >= 0 && lmax <= 3);

    for (int l = 0; l <= lmax; l++)
    {
        double legendre_l = gsl_sf_legendre_Pl(l, cos_theta);
        pair_phi_e        = PairPhi(l, omega, omega_prime, eta, temp, 0);
        pair_phi_x        = PairPhi(l, omega, omega_prime, eta, temp, 1);

        pair_kernel_production_e +=
            (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) *
            (2. * l + 1.) * pair_phi_e * legendre_l / 2.;
        pair_kernel_production_x +=
            (1. / (1. + filterpar * l * l * (l + 1) * (l + 1))) *
            (2. * l + 1.) * pair_phi_x * legendre_l / 2.;
    }

    pair_kernel_absorption_e =
        exp((omega + omega_prime) / temp) * pair_kernel_production_e;
    pair_kernel_absorption_x =
        exp((omega + omega_prime) / temp) * pair_kernel_production_x;

    MyKernelQuantity pair_kernel = {
        .em_e  = pair_kernel_production_e,
        .abs_e = pair_kernel_absorption_e,
        .em_x  = pair_kernel_production_x,
        .abs_x = pair_kernel_absorption_x
    };

    return pair_kernel;
}

/* Calculates the production and absorption kernels for the pair process for M1
 * (l = 0)
 *
 * R^p_{TP}(omega,omega',cos theta) = (1/2) Phi_l(omega,omega')
 * R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T}
 * R^p_{TP}(omega,omega',cos theta)
 *
 */
KOKKOS_INLINE_FUNCTION
MyKernelQuantity PairKernelsM1(MyEOSParams* eos_pars,
                               PairKernelParams* kernel_pars)
{
    // EOS specific parameters
    double eta  = eos_pars->mu_e / eos_pars->temp;
    double temp = eos_pars->temp;

    // kernel specific parameters
    double omega       = kernel_pars->omega;
    double omega_prime = kernel_pars->omega_prime;

    double pair_kernel_production_e = 0.;
    double pair_kernel_absorption_e = 0.;
    double pair_kernel_production_x = 0.;
    double pair_kernel_absorption_x = 0.;

    double pair_phi_e = 0.;
    double pair_phi_x = 0.;

    pair_phi_e = PairPhi(0, omega, omega_prime, eta, temp, 0);
    pair_phi_x = PairPhi(0, omega, omega_prime, eta, temp, 1);

    pair_kernel_production_e = 0.5 * pair_phi_e;
    pair_kernel_production_x = 0.5 * pair_phi_x;

    pair_kernel_absorption_e =
        SafeExp((omega + omega_prime) / temp) * pair_kernel_production_e;
    pair_kernel_absorption_x =
        SafeExp((omega + omega_prime) / temp) * pair_kernel_production_x;

    MyKernelQuantity pair_kernel = {
        .em_e  = pair_kernel_production_e,
        .abs_e = pair_kernel_absorption_e,
        .em_x  = pair_kernel_production_x,
        .abs_x = pair_kernel_absorption_x};

    return pair_kernel;
}

//=======================
// Deprecated functions
//=======================
KOKKOS_INLINE_FUNCTION
void TabulatePairTFunction(double xi, double xf, double x[dim_pair_t],
                           double t[6][dim_pair_t])
{
    const double dx = (xf - xi) / ((double)dim_pair_t);

    for (int i = 0; i < dim_pair_t; i++)
    {
        x[i] = xi + (double)i * dx;
        for (int j = 1; j < 7; j++)
        {
            t[j - 1][i] = PairT(j, x[i], tol_t);
        }
    }

    return;
}

KOKKOS_INLINE_FUNCTION
void PairTInterpolated(PairKernelParams* kernel_pars, double alpha, double* out)
{

    assert(alpha >= 0);

    const double dx =
        kernel_pars->alpha[1] - kernel_pars->alpha[0]; // uniform 1d grid
    int idx = floor((alpha - kernel_pars->alpha[0]) / dx);

    idx = (idx < 0) ? 0 : idx;
    idx = (idx > dim_pair_t - 2) ? dim_pair_t : idx;

    assert(alpha >= kernel_pars->alpha[idx]);

    const double r = (alpha - kernel_pars->alpha[idx]) /
                     (kernel_pars->alpha[idx + 1] - kernel_pars->alpha[idx]);

    for (int i = 0; i < 6; i++)
    {
        out[i] = kernel_pars->pair_t[i][idx] * (1. - r) +
                 kernel_pars->pair_t[i][idx + 1] * r;
    }

    return;
}

KOKKOS_INLINE_FUNCTION
double PairFInterpolated(PairKernelParams* kernel_pars, int k, double eta,
                         double x1)
{
    double pair_t_1[6], pair_t_2[6];

    PairTInterpolated(kernel_pars, fabs(eta), pair_t_1); // PairT(|eta|)

    double result = 0.;

    if (eta < 0.)
    {
        double sum = 0.;
        PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1 - eta)
        for (int l = 0; l <= k; l++)
        {
            sum += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
        }
        result = tgamma(k + 1) * (pair_t_1[k + 1 - 1] - sum);
    }
    else if (k != 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        PairTInterpolated(kernel_pars, 0., pair_t_2); // PairT(0.)
        for (int l = 0; l <= (int)((k - 1.) / 2.); l++)
        {
            sum += pair_t_2[2 * l + 2 - 1] * pow(eta, k - 1. - 2. * l) /
                   tgamma(k - 1 - 2 * l + 1);
        }

        double sum_2 = 0.;
        PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1-eta)
        for (int l = 0; l <= k; l++)
        {
            sum_2 += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 tgamma(k + 1) *
                     (2.0 * sum + pow(-1., k) * pair_t_1[k + 1 - 1] - sum_2);
    }
    else if (k == 0 && 0. <= eta && eta <= x1)
    {

        double sum = 0.;
        PairTInterpolated(kernel_pars, x1 - eta, pair_t_2); // PairT(x1-eta)
        for (int l = 0; l <= k; l++)
        {
            sum += pair_t_2[k + 1 - l - 1] * pow(x1, l) / tgamma(l + 1);
        }

        result = pow(eta, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * pair_t_1[k + 1 - 1] - sum);
    }
    else
    {

        double sum = 0.;
        PairTInterpolated(kernel_pars, eta - x1, pair_t_2); // PairT(eta-x1)
        for (int l = 0; l <= k; l++)
        {
            sum += pow(-1., k - l) * pair_t_2[k + 1 - l - 1] * pow(x1, l) /
                   tgamma(l + 1);
        }
        result = pow(x1, k + 1.) / (k + 1.) +
                 tgamma(k + 1) * (pow(-1., k) * pair_t_1[k + 1 - 1] - sum);
    }

    return result;
}
#endif // GSL_INCLUDES_H_
