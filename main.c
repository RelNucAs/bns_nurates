#include <math.h>
#include <stdio.h>
#include "src/constants.h"
#include "src/bns_nurates.h"
#include "src/integration/integration.h"
#include "tests/tests.h"
#include "src/opacities/kernels/kernels.h"
#include "src/functions/functions.h"
#include "src/opacities/opacities.h"

double arctan(double *x, void *params) {
  return atan(x[0]);
}

double arctan_der(double *x, void *params) {
  return 1. / (1. + x[0] * x[0]);
}

int main() {

  printf("Hello, World!\n");

  double nb = 3.730000e+14 / kMb;
  double temp = 1.204000e+01;
  double lep_mass = kMe;
  double yp = 3.134000e-01;
  double yn = 6.866000e-01;
  double mu_l = 2.495089e+02;
  double mu_hat = 6.015553e+01;
  double omega = 10.08;


  OpacityParams opacity_pars = {.omega = omega,
                                .use_dU = false,
                                .use_WM_ab = false};

  MyEOSParams EOS_pars = {.nb = nb,
                          .temp = temp,
                          .yp = yp,
                          .yn = yn,
                          .mu_e = mu_l,
                          .mu_mu = 0.,
                          .mu_p = 0.,
                          .mu_n = mu_hat + kQ,
                          .dU = 0.};

  //ElasticScattParams ker_pars = {.omega = omega,
  //                               .mu = 0.,
  //                               .mu_prime = 1.,
  //                               .use_WM_sc = false};

  MyOpacity test = AbsOpacity(&opacity_pars, &EOS_pars);
  printf("nu_abs_em output:\n");
  printf("%.5e, %.5e\n", test.em_nue , test.ab_nue /kClight);
  printf("%.5e, %.5e\n", test.em_anue, test.ab_anue/kClight);
  printf("%.5e, %.5e\n", test.em_nux , test.ab_nux /kClight);

  printf("\n");

  SourceCoeffs coeff_test = NuclAbsEmissionCoeffs(&EOS_pars);
  printf("emission coefficients output:\n");
  printf("R_nue = %.5e, R_anue = %.5e, R_nux = %.5e\n", coeff_test.R_nue, coeff_test.R_anue , coeff_test.R_nux);
  printf("Q_nue = %.5e, Q_anue = %.5e, Q_nux = %.5e\n", coeff_test.Q_nue, coeff_test.Q_anue , coeff_test.Q_nux);
  printf("\n");

  //MyKernel test_2 =  nu_N_scatt_kern_tot(&ker_pars, &EOS_pars);

  //printf("nu_N_scatt_kern_tot output:\n");
  //printf("%.5e, %.5e\n", test_2.absorption_e, test_2.production_e);
  //printf("%.5e, %.5e\n", test_2.absorption_x, test_2.production_x);

  //TestQuadratureInputOutput("/var/home/maitraya/Documents/bns_nurates/src/quadratures/");

  //TestGaussLegendreQuadratureMultiD();
  //TestQuadratureWithGSL();
  TestIntegrationMultiD();
  //PrintGaussLegendreQuadrature(10, -1., 1.);
  //printf("\n");
  //PrintGaussLaguerreQuadrature(5, 0.);
  //TestQuadratureWithGSL();
  //TestPairT();
  //printf("\n");
  //TestPairF();
  //printf("\n");
  //TestPairG();
  //printf("\n");
  //TestPairPsi();
  //printf("\n");
  //TestPairPhi();
  //printf("\n");
  //TestPairKernels();
  //printf("\n");
  //TestPairOpacities();
  //printf("\n");
  //TestBremKernelS("/var/home/maitraya/Documents/");

  //BremKernelS(-0.5, 0., 0.31);
  //TestBremKernelG("/var/home/maitraya/Documents/");
  //MyQuadrature quad;
  //quad.n = 3;
  //quad.dim = 1;
  //quad.type = kGaulag;
  //quad.x1 = -1.;
  //quad.x2 = 1.;
  //quad.alpha = 0.;

  //SaveQuadrature("/var/home/maitraya/Documents/bns_nurates/src/quadratures/", &quad);

  MyFunction f_NR;
  f_NR.dim = 1;
  f_NR.function = &arctan;
  
  MyFunction df_NR;
  df_NR.dim = 1;
  df_NR.function = &arctan_der;

  double test_NR = MNewt1d(3.14, -5., +5., 3.14/4., &f_NR, &df_NR);
  printf("Result of 1D Newton-Raphson is %.5lf\n", test_NR);

  return 0;
}
