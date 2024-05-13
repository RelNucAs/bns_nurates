//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.c
//  \brief contains pair kernels and associated helper functions

#include <math.h>
#include <assert.h>
#include <stddef.h>
#include "kernels.h"
#include "functions.h"
#include "constants.h"

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

void F_for_G(double eta, double y, double z, double* F_out)
{
    F_out[0 + 0 * 8] = FDI_0(eta);
    F_out[1 + 0 * 8] = FDI_0(eta - y);
    F_out[2 + 0 * 8] = FDI_0(eta - z);
    F_out[3 + 0 * 8] = FDI_0(eta + y);
    F_out[4 + 0 * 8] = FDI_0(eta + z);
    F_out[5 + 0 * 8] = FDI_0(eta - y - z);
    F_out[6 + 0 * 8] = (y <= z) ? F_out[1] : F_out[2];
    F_out[7 + 0 * 8] = (y <= z) ? F_out[4] : F_out[3];

    F_out[0 + 1 * 8] = FDI_p1(eta);
    F_out[1 + 1 * 8] = FDI_p1(eta - y);
    F_out[2 + 1 * 8] = FDI_p1(eta - z);
    F_out[3 + 1 * 8] = FDI_p1(eta + y);
    F_out[4 + 1 * 8] = FDI_p1(eta + z);
    F_out[5 + 1 * 8] = FDI_p1(eta - y - z);
    F_out[6 + 1 * 8] = (y <= z) ? F_out[1 + 1 * 8] : F_out[2 + 1 * 8];
    F_out[7 + 1 * 8] = (y <= z) ? F_out[4 + 1 * 8] : F_out[3 + 1 * 8];

    F_out[0 + 2 * 8] = FDI_p2(eta);
    F_out[1 + 2 * 8] = FDI_p2(eta - y);
    F_out[2 + 2 * 8] = FDI_p2(eta - z);
    F_out[3 + 2 * 8] = FDI_p2(eta + y);
    F_out[4 + 2 * 8] = FDI_p2(eta + z);
    F_out[5 + 2 * 8] = FDI_p2(eta - y - z);
    F_out[6 + 2 * 8] = (y <= z) ? F_out[1 + 2 * 8] : F_out[2 + 2 * 8];
    F_out[7 + 2 * 8] = (y <= z) ? F_out[4 + 2 * 8] : F_out[3 + 2 * 8];

    F_out[0 + 3 * 8] = FDI_p3(eta);
    F_out[1 + 3 * 8] = FDI_p3(eta - y);
    F_out[2 + 3 * 8] = FDI_p3(eta - z);
    F_out[3 + 3 * 8] = FDI_p3(eta + y);
    F_out[4 + 3 * 8] = FDI_p3(eta + z);
    F_out[5 + 3 * 8] = FDI_p3(eta - y - z);
    F_out[6 + 3 * 8] = (y <= z) ? F_out[1 + 3 * 8] : F_out[2 + 3 * 8];
    F_out[7 + 3 * 8] = (y <= z) ? F_out[4 + 3 * 8] : F_out[3 + 3 * 8];

    F_out[0 + 4 * 8] = FDI_p4(eta);
    F_out[1 + 4 * 8] = FDI_p4(eta - y);
    F_out[2 + 4 * 8] = FDI_p4(eta - z);
    F_out[3 + 4 * 8] = FDI_p4(eta + y);
    F_out[4 + 4 * 8] = FDI_p4(eta + z);
    F_out[5 + 4 * 8] = FDI_p4(eta - y - z);
    F_out[6 + 4 * 8] = (y <= z) ? F_out[1 + 4 * 8] : F_out[2 + 4 * 8];
    F_out[7 + 4 * 8] = (y <= z) ? F_out[4 + 4 * 8] : F_out[3 + 4 * 8];

    F_out[0 + 5 * 8] = FDI_p5(eta);
    F_out[1 + 5 * 8] = FDI_p5(eta - y);
    F_out[2 + 5 * 8] = FDI_p5(eta - z);
    F_out[3 + 5 * 8] = FDI_p5(eta + y);
    F_out[4 + 5 * 8] = FDI_p5(eta + z);
    F_out[5 + 5 * 8] = FDI_p5(eta - y - z);
    F_out[6 + 5 * 8] = (y <= z) ? F_out[1 + 5 * 8] : F_out[2 + 5 * 8];
    F_out[7 + 5 * 8] = (y <= z) ? F_out[4 + 5 * 8] : F_out[3 + 5 * 8];

    return;
}

void PairGOptimized(double eta, double y, double z, double* g_out)
{
    double m = MIN(y, z);

    double F[48];

    F_for_G(eta, y, z, F);

    g_out[0] = (1 * (POW0(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                     POW0(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8])));
    g_out[1] = (1 * (POW0(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                     POW0(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8])));
    g_out[2] = (y <= z) ? g_out[1] : g_out[0];
    g_out[3] = (-1 * (POW0(m) * (F[6 + 0 * 8] - F[7 + 0 * 8])));

    g_out[4] = (1 * (POW1(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                     POW1(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
               (1 * (POW0(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) -
                     POW0(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8])));
    g_out[5] = (1 * (POW1(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                     POW1(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
               (1 * (POW0(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) -
                     POW0(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8])));
    g_out[6] = (y <= z) ? g_out[5] : g_out[4];
    g_out[7] = (-1 * (POW1(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) +
               (-1 * (POW0(m) * (F[6 + 1 * 8] - F[7 + 1 * 8])));

    g_out[8] = (1 * (POW2(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                     POW2(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
               (2 * (POW1(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) -
                     POW1(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
               (1 * (POW0(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) -
                     POW0(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8])));
    g_out[9] = (1 * (POW2(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                     POW2(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
               (2 * (POW1(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) -
                     POW1(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
               (1 * (POW0(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) -
                     POW0(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8])));
    g_out[10] = (y <= z) ? g_out[9] : g_out[8];
    g_out[11] = (-1 * (POW2(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) +
                (-2 * (POW1(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) +
                (-1 * (POW0(m) * (F[6 + 2 * 8] - F[7 + 2 * 8])));

    g_out[12] = (1 * (POW3(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                      POW3(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (3 * (POW2(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) -
                      POW2(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (3 * (POW1(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) -
                      POW1(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (1 * (POW0(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) -
                      POW0(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8])));
    g_out[13] = (1 * (POW3(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                      POW3(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (3 * (POW2(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) -
                      POW2(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (3 * (POW1(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) -
                      POW1(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (1 * (POW0(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) -
                      POW0(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8])));
    g_out[14] = (y <= z) ? g_out[13] : g_out[12];
    g_out[15] = (-1 * (POW3(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) +
                (-3 * (POW2(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) +
                (-3 * (POW1(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) +
                (-1 * (POW0(m) * (F[6 + 3 * 8] - F[7 + 3 * 8])));

    g_out[16] = (1 * (POW4(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                      POW4(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (4 * (POW3(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) -
                      POW3(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (6 * (POW2(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) -
                      POW2(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (4 * (POW1(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) -
                      POW1(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) +
                (1 * (POW0(y) * (F[1 + 4 * 8] - F[4 + 4 * 8]) -
                      POW0(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8])));
    g_out[17] = (1 * (POW4(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                      POW4(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (4 * (POW3(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) -
                      POW3(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (6 * (POW2(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) -
                      POW2(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (4 * (POW1(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) -
                      POW1(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) +
                (1 * (POW0(z) * (F[2 + 4 * 8] - F[3 + 4 * 8]) -
                      POW0(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8])));
    g_out[18] = (y <= z) ? g_out[17] : g_out[16];
    g_out[19] = (-1 * (POW4(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) +
                (-4 * (POW3(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) +
                (-6 * (POW2(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) +
                (-4 * (POW1(m) * (F[6 + 3 * 8] - F[7 + 3 * 8]))) +
                (-1 * (POW0(m) * (F[6 + 4 * 8] - F[7 + 4 * 8])));

    g_out[20] = (1 * (POW5(y) * (F[1 + 0 * 8] - F[4 + 0 * 8]) -
                      POW5(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (5 * (POW4(y) * (F[1 + 1 * 8] - F[4 + 1 * 8]) -
                      POW4(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (10 * (POW3(y) * (F[1 + 2 * 8] - F[4 + 2 * 8]) -
                       POW3(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (10 * (POW2(y) * (F[1 + 3 * 8] - F[4 + 3 * 8]) -
                       POW2(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) +
                (5 * (POW1(y) * (F[1 + 4 * 8] - F[4 + 4 * 8]) -
                      POW1(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))) +
                (1 * (POW0(y) * (F[1 + 5 * 8] - F[4 + 5 * 8]) -
                      POW0(y + z) * (F[5 + 5 * 8] - F[0 + 5 * 8])));
    g_out[21] = (1 * (POW5(z) * (F[2 + 0 * 8] - F[3 + 0 * 8]) -
                      POW5(y + z) * (F[5 + 0 * 8] - F[0 + 0 * 8]))) +
                (5 * (POW4(z) * (F[2 + 1 * 8] - F[3 + 1 * 8]) -
                      POW4(y + z) * (F[5 + 1 * 8] - F[0 + 1 * 8]))) +
                (10 * (POW3(z) * (F[2 + 2 * 8] - F[3 + 2 * 8]) -
                       POW3(y + z) * (F[5 + 2 * 8] - F[0 + 2 * 8]))) +
                (10 * (POW2(z) * (F[2 + 3 * 8] - F[3 + 3 * 8]) -
                       POW2(y + z) * (F[5 + 3 * 8] - F[0 + 3 * 8]))) +
                (5 * (POW1(z) * (F[2 + 4 * 8] - F[3 + 4 * 8]) -
                      POW1(y + z) * (F[5 + 4 * 8] - F[0 + 4 * 8]))) +
                (1 * (POW0(z) * (F[2 + 5 * 8] - F[3 + 5 * 8]) -
                      POW0(y + z) * (F[5 + 5 * 8] - F[0 + 5 * 8])));
    g_out[22] = (y <= z) ? g_out[21] : g_out[20];
    g_out[23] = (-1 * (POW5(m) * (F[6 + 0 * 8] - F[7 + 0 * 8]))) +
                (-5 * (POW4(m) * (F[6 + 1 * 8] - F[7 + 1 * 8]))) +
                (-10 * (POW3(m) * (F[6 + 2 * 8] - F[7 + 2 * 8]))) +
                (-10 * (POW2(m) * (F[6 + 3 * 8] - F[7 + 3 * 8]))) +
                (-5 * (POW1(m) * (F[6 + 4 * 8] - F[7 + 4 * 8]))) +
                (-1 * (POW0(m) * (F[6 + 5 * 8] - F[7 + 5 * 8])));

    return;
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
// @TODO: compute only a coeffs needed for the specific l value
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

void PairPsiOptimized(int l, double y, double z, double eta, double* psi_out)
{
    assert(l == 0);

    double G[24];
    double a_yz[4][12], c_yz[4][3], d_yz[4][3];
    double a_zy[4][12], c_zy[4][3], d_zy[4][3];

    PairPsiCoeffs(y, z, a_yz, c_yz, d_yz);
    PairPsiCoeffs(z, y, a_zy, c_zy, d_zy);

    PairGOptimized(eta, y, z, G);

    psi_out[0] = c_yz[0][0] * G[0] + d_yz[0][0] * G[1] + c_yz[0][1] * G[4] +
                 d_yz[0][1] * G[5] + c_yz[0][2] * G[8] + d_yz[0][2] * G[9] +
                 a_yz[0][3] * (G[15] - G[14]) + a_yz[0][4] * (G[19] - G[18]) +
                 a_yz[0][5] * (G[23] - G[22]); // Psi(y,z)
    psi_out[1] = c_zy[0][0] * G[1] + d_zy[0][0] * G[0] + c_zy[0][1] * G[5] +
                 d_zy[0][1] * G[4] + c_zy[0][2] * G[9] + d_zy[0][2] * G[8] +
                 a_zy[0][3] * (G[15] - G[14]) + a_zy[0][4] * (G[19] - G[18]) +
                 a_zy[0][5] * (G[23] - G[22]); // Psi(z,y)

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
void PairPhiOptimized(const double omega, const double omega_prime, int l,
                      double eta, double temp, double* phi_out)
{
    static const double kPairPhi = kGSqr / kPi;

    const double y = omega / temp;
    const double z = omega_prime / temp;

    const double phi_prefactor = kPairPhi * temp * temp;
    const double phi_denom     = 1. - exp(y + z);

    double pair_psi[2] = {0.};

    PairPsiOptimized(l, y, z, eta, pair_psi);

    phi_out[0] = phi_prefactor *
                 (kAlpha1[0] * kAlpha1[0] * pair_psi[0] +
                  kAlpha2[0] * kAlpha2[0] * pair_psi[1]) /
                 phi_denom;
    phi_out[1] = phi_prefactor *
                 (kAlpha1[0] * kAlpha1[0] * pair_psi[1] +
                  kAlpha2[0] * kAlpha2[0] * pair_psi[0]) /
                 phi_denom;
    phi_out[2] = phi_prefactor *
                 (kAlpha1[1] * kAlpha1[1] * pair_psi[0] +
                  kAlpha2[1] * kAlpha2[1] * pair_psi[1]) /
                 phi_denom;
    phi_out[3] = phi_prefactor *
                 (kAlpha1[1] * kAlpha1[1] * pair_psi[1] +
                  kAlpha2[1] * kAlpha2[1] * pair_psi[0]) /
                 phi_denom;

    return;
}


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

void PairKernelsTable(const int n, double* nu_array,
                      GreyOpacityParams* grey_pars, M1Matrix* out)
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
double PairG(int n, double a, double b, double eta, double y, double z)
{

    double result = 0.;
    result = PairF(n, eta, b) - PairF(n, eta, a) - PairF(n, eta + y + z, b) +
             PairF(n, eta + y + z, a);

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

    MyKernelQuantity pair_kernel = {.abs_e = pair_kernel_absorption_e,
                                    .em_e  = pair_kernel_production_e,
                                    .abs_x = pair_kernel_absorption_x,
                                    .em_x  = pair_kernel_production_x};

    return pair_kernel;
}

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
        .abs_e = kPrefactor * pair_kernel_absorption_e,
        .em_e  = kPrefactor * pair_kernel_production_e,
        .abs_x = kPrefactor * pair_kernel_absorption_x,
        .em_x  = kPrefactor * pair_kernel_production_x};

    return pair_kernel;
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
 * Eqns. (2) and (3) of Pons et. al.
 *
 * R^p_{TP}(omega,omega',cos theta) = sum_{l} ((2l+1)/2) Phi_l(omega,omega')
 * P_l(cos theta) R^a_{TP}(omega,omega',cos theta) = e^{(omega+omega')/T}
 * R^p_{TP}(omega,omega',cos theta)
 *
 * l is an integer.
 */
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

    MyKernelQuantity pair_kernel = {.abs_e = pair_kernel_absorption_e,
                                    .em_e  = pair_kernel_production_e,
                                    .abs_x = pair_kernel_absorption_x,
                                    .em_x  = pair_kernel_production_x};

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

    MyKernelQuantity pair_kernel = {.abs_e = pair_kernel_absorption_e,
                                    .em_e  = pair_kernel_production_e,
                                    .abs_x = pair_kernel_absorption_x,
                                    .em_x  = pair_kernel_production_x};

    return pair_kernel;
}

//=======================
// Deprecated functions
//=======================

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
