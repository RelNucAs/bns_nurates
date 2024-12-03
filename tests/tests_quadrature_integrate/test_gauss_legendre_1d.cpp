//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_gauss_legendre_1d.c
//  \brief Test 1d Gauss-Legendre routines with known values

#include <stdio.h>
#include <math.h>
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

int main()
{

    double points_2[2] = {-(1. / 3.) * sqrt(3.), (1. / 3.) * sqrt(3.)};
    double w_2[2]      = {1., 1.};
    double points_3[3] = {-(1. / 5.) * sqrt(15.), 0, (1. / 5.) * sqrt(15.)};
    double w_3[3]      = {5. / 9., 8. / 9., 5. / 9.};
    double points_4[4] = {-(1. / 35.) * sqrt(525. + 70. * sqrt(30.)),
                          -(1. / 35.) * sqrt(525. - 70. * sqrt(30.)),
                          (1. / 35.) * sqrt(525. - 70. * sqrt(30.)),
                          (1. / 35.) * sqrt(525. + 70. * sqrt(30.))};
    double w_4[4]      = {
        (1. / 36.) * (18. - sqrt(30.)), (1. / 36.) * (18. + sqrt(30.)),
        (1. / 36.) * (18. + sqrt(30.)), (1. / 36.) * (18. - sqrt(30.))};
    double max_error = -42.;

    int n     = 2;
    double x1 = -1.;
    double x2 = 1.;

    MyQuadrature quad = quadrature_default;
    quad.dim          = 1;
    quad.type         = kGauleg;
    quad.x1           = x1;
    quad.x2           = x2;
    quad.nx           = n;
    GaussLegendre(&quad);

    printf("Generating Gauss-Legendre quadratures ... \n");
    printf("x1 = %f \t x2 = %f \t nx = %d\n", quad.x1, quad.x2, quad.nx);
    printf("points \t\t answer \t weights \t answer\n");
    for (int i = 0; i < quad.nx; i++)
    {
        printf("%f \t %f \t %f \t %f\n", quad.points[i], points_2[i], quad.w[i],
               w_2[i]);
        if (fabs(quad.points[i] - points_2[i]) > max_error)
        {
            max_error = fabs(quad.points[i] - points_2[i]);
        }
        if (fabs(quad.w[i] - w_2[i]) > max_error)
        {
            max_error = fabs(quad.w[i] - w_2[i]);
        }
    }

    n       = 3;
    quad.nx = n;

    GaussLegendre(&quad);

    printf("Generating Gauss-Legendre quadratures ... \n");
    printf("x1 = %f \t x2 = %f \t nx = %d\n", quad.x1, quad.x2, quad.nx);
    printf("points \t\t answer \t weights \t answer\n");
    for (int i = 0; i < quad.nx; i++)
    {
        printf("%f \t %f \t %f \t %f\n", quad.points[i], points_3[i], quad.w[i],
               w_3[i]);
        if (fabs(quad.points[i] - points_3[i]) > max_error)
        {
            max_error = fabs(quad.points[i] - points_3[i]);
        }
        if (fabs(quad.w[i] - w_3[i]) > max_error)
        {
            max_error = fabs(quad.w[i] - w_3[i]);
        }
    }

    n       = 4;
    quad.nx = n;

    GaussLegendre(&quad);

    printf("Generating Gauss-Legendre quadratures ... \n");
    printf("x1 = %f \t x2 = %f \t nx = %d\n", quad.x1, quad.x2, quad.nx);
    printf("points \t\t answer \t weights \t answer\n");
    for (int i = 0; i < quad.nx; i++)
    {
        printf("%f \t %f \t %f \t %f\n", quad.points[i], points_4[i], quad.w[i],
               w_4[i]);
        if (fabs(quad.points[i] - points_4[i]) > max_error)
        {
            max_error = fabs(quad.points[i] - points_4[i]);
        }
        if (fabs(quad.w[i] - w_4[i]) > max_error)
        {
            max_error = fabs(quad.w[i] - w_4[i]);
        }
    }

    if (max_error < 1e-12)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}