//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_gauss_laguerre_1d.c
//  \brief Test 1d Gauss-Laguerre routine with known values

#include <stdio.h>
#include <math.h>
#include "../../include/bns_nurates.hpp"
#include "../../include/integration.hpp"

int main() {

  double max_error = -42.;

  double points_2[2] = {0.585786437626904951, 3.41421356237309505};
  double w_2[2] = {0.853553390593273762, 0.146446609406726238};
  double points_3[3] = {0.415774556783479083, 2.29428036027904172, 6.2899450829374792};
  double w_3[3] = {0.711093009929173015, 0.278517733569240849, 0.0103892565015861358};
  double points_4[4] = {0.322547689619392312, 1.74576110115834658, 4.53662029692112798, 9.39507091230113313};
  double w_4[4] = {0.603154104341633602, 0.357418692437799687, 0.0388879085150053843, 5.3929470556132745E-4};

  printf("========================================= \n");
  printf("Generating Gauss-Laguerre quadratures ... \n");
  printf("========================================= \n");

  int n = 2;
  double alpha = 0.;

  MyQuadrature quad = quadrature_default;
  quad.dim =1;
  quad.type = kGaulag;
  quad.nx = n;
  quad.alpha = alpha;
  GaussLaguerre(&quad);

  printf("x1 = 0 \t\t x2 = inf \t nx = %d\n", quad.nx);
  printf("points \t\t answer \t weights \t answer\n");
  for (int i = 0; i < quad.nx; i++) {
    printf("%f \t %f \t %f \t %f\n", quad.points[i], points_2[i], quad.w[i], w_2[i]);
    if (fabs(quad.points[i] - points_2[i]) > max_error) {
      max_error = fabs(quad.points[i] - points_2[i]);
    }
    if (fabs(quad.w[i] - w_2[i]) > max_error) {
      max_error = fabs(quad.w[i] - w_2[i]);
    }
  }

  quad.nx = 3;
  GaussLaguerre(&quad);

  printf("\n");
  printf("x1 = 0 \t\t x2 = inf \t nx = %d\n", quad.nx);
  printf("points \t\t answer \t weights \t answer\n");
  for (int i = 0; i < quad.nx; i++) {
    printf("%f \t %f \t %f \t %f\n", quad.points[i], points_3[i], quad.w[i], w_3[i]);
    if (fabs(quad.points[i] - points_3[i]) > max_error) {
      max_error = fabs(quad.points[i] - points_3[i]);
    }
    if (fabs(quad.w[i] - w_3[i]) > max_error) {
      max_error = fabs(quad.w[i] - w_3[i]);
    }
  }

  quad.nx = 4;
  GaussLaguerre(&quad);

  printf("\n");
  printf("x1 = 0 \t\t x2 = inf \t nx = %d\n", quad.nx);
  printf("points \t\t answer \t weights \t answer\n");
  for (int i = 0; i < quad.nx; i++) {
    printf("%f \t %f \t %f \t %f\n", quad.points[i], points_4[i], quad.w[i], w_4[i]);
    if (fabs(quad.points[i] - points_4[i]) > max_error) {
      max_error = fabs(quad.points[i] - points_4[i]);
    }
    if (fabs(quad.w[i] - w_4[i]) > max_error) {
      max_error = fabs(quad.w[i] - w_4[i]);
    }
  }

  printf("\nMax error: %e\n", max_error);

  if (max_error < 1e-5) {
    return 0;
  } else {
    return 1;
  }
}