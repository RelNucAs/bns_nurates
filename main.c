#include <math.h>
#include <stdio.h>
#include "include/constants.h"
//#include "src/bns_nurates.h"
//#include "src/integration/integration.h"
#include "include/kernels.h"
//#include "src/functions/functions.h"
#include "src/opacities/opacities.h"

int main() {

  double nb = 3.7300E+14 / kMb;
  double temp = 1.2040E+01;
  double omega = 10.;
  double omega_prime = 20.;

  double x = (omega + omega_prime) / (temp);
  double y = 19.4 / temp;
  double eta_star = pow(3. * kPi * kPi * nb, 2. / 3.) * (kH * kH * kMeV / (4. * kPi * kPi)) / (2. * kMb * temp);

  double out = BremKernelS(x, y, eta_star);

  printf("Result: %.5e\n", out);


  // activate only iso-energetic scattering
  GreyOpacityParams my_grey_opacity_params;

  my_grey_opacity_params.eos_pars.mu_e = 2.4951E+02;
  my_grey_opacity_params.eos_pars.mu_p = 0.;
  my_grey_opacity_params.eos_pars.mu_n = 6.0156E+01 + kQ;
  my_grey_opacity_params.eos_pars.temp = temp;
  my_grey_opacity_params.eos_pars.yp = 3.1340E-01;
  my_grey_opacity_params.eos_pars.yn = 6.8660E-01;
  my_grey_opacity_params.eos_pars.nb = 3.7300E+14 / kMu;


  my_grey_opacity_params.opacity_pars.use_WM_sc = false;
  my_grey_opacity_params.opacity_pars.use_WM_ab = false;
  my_grey_opacity_params.opacity_pars.use_dU = false;

  // only enable Fermi distribution function for the three neutrino species
  my_grey_opacity_params.distr_pars.temp_t[0] = temp;
  my_grey_opacity_params.distr_pars.temp_t[1] = temp;
  my_grey_opacity_params.distr_pars.temp_t[2] = temp;
  my_grey_opacity_params.distr_pars.w_t[0] = 1.;
  my_grey_opacity_params.distr_pars.w_t[1] = 1.;
  my_grey_opacity_params.distr_pars.w_t[2] = 1.;
  my_grey_opacity_params.distr_pars.eta_t[0] = (my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p) / temp;
  my_grey_opacity_params.distr_pars.eta_t[1] = -(my_grey_opacity_params.eos_pars.mu_e - my_grey_opacity_params.eos_pars.mu_n + my_grey_opacity_params.eos_pars.mu_p) / temp;
  my_grey_opacity_params.distr_pars.eta_t[2] = 0.;


  my_grey_opacity_params.m1_pars.chi[0] = 1.;
  my_grey_opacity_params.m1_pars.chi[1] = 1.;
  my_grey_opacity_params.m1_pars.chi[2] = 1.;


  // disable thin distribution function for the time being
  my_grey_opacity_params.distr_pars.w_f[0] = 0.;
  my_grey_opacity_params.distr_pars.w_f[1] = 0.;
  my_grey_opacity_params.distr_pars.w_f[2] = 0.;
  my_grey_opacity_params.distr_pars.c_f[0] = 1.;
  my_grey_opacity_params.distr_pars.c_f[1] = 1.;
  my_grey_opacity_params.distr_pars.c_f[2] = 1.;
  my_grey_opacity_params.distr_pars.temp_f[0] = temp;
  my_grey_opacity_params.distr_pars.temp_f[1] = temp;
  my_grey_opacity_params.distr_pars.temp_f[2] = temp;

  M1Opacities out2;//BetaNumberIntegral(&my_grey_opacity_params);
  M1Opacities out3;//BetaEnergyIntegral(&my_grey_opacity_params);

  printf("eta_0_nue  = %.5e\n", out2.eta_0[id_nue]);
  printf("eta_0_anue = %.5e\n", out2.eta_0[id_anue]);
  printf("eta_0_nux  = %.5e\n", out2.eta_0[id_nux]);
  printf("kappa_0_a_nue  = %.5e\n", out2.kappa_0_a[id_nue]);
  printf("kappa_0_a_anue = %.5e\n", out2.kappa_0_a[id_anue]);
  printf("kappa_0_a_nux  = %.5e\n", out2.kappa_0_a[id_nux]);

  printf("\n");
  
  printf("eta_nue  = %.5e\n", out3.eta[id_nue]);
  printf("eta_anue = %.5e\n", out3.eta[id_anue]);
  printf("eta_nux  = %.5e\n", out3.eta[id_nux]);
  printf("kappa_a_nue  = %.5e\n", out3.kappa_a[id_nue]);
  printf("kappa_a_anue = %.5e\n", out3.kappa_a[id_anue]);
  printf("kappa_a_nux  = %.5e\n", out3.kappa_a[id_nux]);

  return 0;
}