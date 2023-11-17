#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/constants.h"
#include "include/bns_nurates.h"
#include "include/integration.h"
#include "include/distribution.h"
#include "include/constants.h"
#include "include/functions.h"
#include "include/kernels.h"
#include "include/functions.h"
#include "src/opacities/opacities.h"


int main() {
  // Test inelastic scattering
  printf("======================================== \n");
  printf("Testing inelastic scattering reaction... \n");
  printf("======================================== \n");
  printf("\n");

  const double dens =  1.; // [g cm^-3];
  const double temp = 10.; // [MeV]
  const double ye = 0.1;
  
  // EOS parameters
  MyEOSParams eos_pars;
  eos_pars.nb = dens / kMu;         // [cm-3]
  eos_pars.temp = temp;             // [MeV]
  eos_pars.ye = ye;
  eos_pars.yp = ye;
  eos_pars.yn = 6.8660E-01;
  eos_pars.mu_e = 2.4951E+02;       // [MeV]
  eos_pars.mu_p = 0.;               // [MeV]
  eos_pars.mu_n = 6.0156E+01 + kQ;  // [MeV]

  /////////////
  // Kernels //
  /////////////

  printf("rho  = %.5e g cm^3\n", dens);
  printf("temp = %.5e MeV   \n", temp);
  printf("ye   = %.5e       \n", ye);
  printf("\n");

  // Kernel parameters
  InelasticScattKernelParams inelastic_pars;
  inelastic_pars.omega = 10.;       // [MeV]
  inelastic_pars.omega_prime = 20.; // [MeV]

  printf("E_nu       = %.5lf MeV\n", inelastic_pars.omega);
  printf("E_nu_prime = %.5lf MeV\n", inelastic_pars.omega_prime);
  printf("\n");

  // Compute kernels
  MyKernelOutput nes_kernels = NESKernels(&inelastic_pars, &eos_pars);
  MyKernelOutput nps_kernels = NPSKernels(&inelastic_pars, &eos_pars);

  // Print kernel output
  printf("NES kernels [cm^3 s^-1]:\n");

  printf("nue:  R_in = %.5e, R_out = %.5e\n", nes_kernels.em[id_nue] , nes_kernels.abs[id_nue]);
  printf("anue: R_in = %.5e, R_out = %.5e\n", nes_kernels.em[id_anue], nes_kernels.abs[id_anue]);
  printf("nux:  R_in = %.5e, R_out = %.5e\n", nes_kernels.em[id_nux],  nes_kernels.abs[id_nux]);
  printf("\n");

  printf("NPS kernels [cm^3 s^-1]:\n");

  printf("nue:  R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_nue] , nps_kernels.abs[id_nue]);
  printf("anue: R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_anue], nps_kernels.abs[id_anue]);
  printf("nux:  R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_nux],  nps_kernels.abs[id_nux]);
  printf("\n");

  ///////////////////////////////////////////
  // Emissivity and inverse mean free path //
  ///////////////////////////////////////////

  // Generate quadratures
  MyQuadrature my_quadrature_1d = {.nx = 60, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&my_quadrature_1d);

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = opacity_params_default_none;

  // Distribution parameters
  NuDistributionParams distr_pars = NuEquilibriumParams(&eos_pars); // consider neutrino distribution function at equilibrium
  
  // Opacity flags
  OpacityFlags opacity_flags = opacity_flags_default_none;
  opacity_flags.use_inelastic_scatt = 1; // activate only inelastic scattering
  
  // Kernel parameters
  MyKernelParams kernel_pars;
  kernel_pars.inelastic_kernel_params = inelastic_pars;

  // Neutrino energy
  double e_nu = 10.; // [MeV]

  // M1 Grey Coefficient parameters  
  GreyOpacityParams grey_pars = {.distr_pars = distr_pars,
                                 .eos_pars = eos_pars,
                                 .opacity_pars = opacity_pars,
                                 .kernel_pars = kernel_pars,
                                 .opacity_flags = opacity_flags};

  // Compute emissivity and inverse mean free path
  SpectralOpacities inelastic_output = ComputeSpectralOpacitiesNotStimulated(e_nu, &my_quadrature_1d, &grey_pars);

  // Print emissivity and inverse mean free path
  printf("Emissivity [s-1] and inverse mean free path [cm-1] (E = %.2lf MeV):\n", e_nu);

  printf("nue:  eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_nue], inelastic_output.kappa[id_nue]);
  printf("anue: eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_anue], inelastic_output.kappa[id_anue]);
  printf("nux:  eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_nux], inelastic_output.kappa[id_nux]);
  printf("\n");

  ///////////////////////////
  // 1D CCSN input profile //
  ///////////////////////////
 
  printf("Emissivity [s-1] and inverse mean free path along 1D CCSN input profile:\n");
  printf("\n");

  char filepath[300] = {'\0'};
  char filedir[300] = SOURCE_DIR;
  char outname[200] = "/tests/tests_opacities_m1/nurates_CCSN/nurates_1.008E+01.txt";

  strcat(filepath, filedir);
  strcat(filepath, outname);

  //printf("Data_directory: %s\n", filepath);

  FILE *fptr;
  fptr = fopen(filepath, "r");
  if (fptr == NULL) {
    printf("%s: The file %s does not exist!\n", __FILE__, filepath);
    exit(1);
  }

  // store columns here
  int num_data = 102;
  int zone[num_data];
  double r[num_data];
  double rho[num_data];
  double T[num_data];
  double Ye[num_data];
  double mu_e[num_data];
  double mu_hat[num_data];
  double Yh[num_data];
  double Ya[num_data];
  double Yp[num_data];
  double Yn[num_data];
  double em_nue[num_data];
  double l_nue_inv[num_data];
  double em_anue[num_data];
  double l_anue_inv[num_data];

  // read in the data file
  int i = 0;
  char line[1000];
  while (fgets(line, sizeof(line), fptr) != NULL) {
    if (line[1] == '#' && i == 0) {
      sscanf(line + 14, "%lf\n", &e_nu);
      continue;
    } else if (line[1] == '#' && i != 0) {
      continue;
    }
    
    sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &zone[i], &r[i], &rho[i], &T[i], &Ye[i], &mu_e[i], &mu_hat[i], &Yh[i], &Ya[i], &Yp[i], &Yn[i], &em_nue[i], &l_nue_inv[i], &em_anue[i], &l_anue_inv[i]);
    i++;
  }

  // print input data
  if (false) {
    printf("\n");
    printf("Printing out input table ... %s\n", outname);
    printf("\n");

    printf("Input data:\n");
    printf("E_nu: %lf\n", e_nu);
    for (int i = 0; i < num_data; i++) {
      printf("%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
             zone[i], r[i], rho[i], T[i], Ye[i], mu_e[i], mu_hat[i], Yh[i], Ya[i], Yp[i], Yn[i], em_nue[i], l_nue_inv[i], em_anue[i], l_anue_inv[i]);
    }

    printf("\n");
    printf("End of input table.\n");
    printf("\n");
  }

  printf("E_nu = %.5lf MeV\n", e_nu);
  printf("\n");

  printf("Generated table:\n");
  printf("r j-nue j-anue j-nux kappa-nue kappa-anue kappa-nux\n");

  for (int i = 0; i < 102; i++) {

    // populate EOS parameters from table
    grey_pars.eos_pars.mu_e = mu_e[i];
    grey_pars.eos_pars.mu_p = 0.;
    grey_pars.eos_pars.mu_n = mu_hat[i] + kQ;
    grey_pars.eos_pars.temp = T[i];
    grey_pars.eos_pars.yp = Yp[i];
    grey_pars.eos_pars.yn = Yn[i];
    grey_pars.eos_pars.ye = Ye[i];
    grey_pars.eos_pars.nb = rho[i] / kMb;
 
    // Neutrino distribution function at equilibrium
    grey_pars.distr_pars = NuEquilibriumParams(&grey_pars.eos_pars); 
    
    inelastic_output = ComputeSpectralOpacitiesNotStimulated(e_nu, &my_quadrature_1d, &grey_pars);

    printf("%e %e %e %e %e %e %e\n",
            r[i], 
            inelastic_output.j[id_nue], inelastic_output.j[id_anue], inelastic_output.j[id_nux],
            inelastic_output.kappa[id_nue], inelastic_output.kappa[id_anue], inelastic_output.kappa[id_nux]
          );
    
  }

  return 0;
}