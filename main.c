#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_zeta.h>

//#include <hdf5.h>
//#include <hdf5_hl.h>
#include <assert.h>

#include "include/constants.h"
#include "include/bns_nurates.h"
#include "include/integration.h"
#include "include/distribution.h"
#include "include/constants.h"
#include "include/functions.h"
#include "include/kernels.h"
#include "include/functions.h"
#include "include/opacities.h"

#define n_points 4
#define n_bins 128

void save_rates_to_file(const int id_point, const char* reac_type, const int n_leg, MyEOSParams *eos_pars, OpacityFlags *opacity_flags, const double rho, const double Abar, const double Zbar, const double mu_hat, const bool neglect_blocking) {
  const char folder[200] = "/home/leonardo/Desktop/PhD_work/BNS_muons/bns_refactor/bns_nurates/tests/tests_code_comparison/";
  char filename[300], full_path[300];
  char file_type[100];

  if (strcmp(reac_type, "pair") || strcmp(reac_type, "brem") ||strcmp(reac_type,"all")) {
    if (neglect_blocking) {
      if (kirchoff_flag) return;
      strcpy(file_type, "_noblock");
    } else {
      if (kirchoff_flag) {
        strcpy(file_type, "_kirchoff");
      } else {
        strcpy(file_type, "");
      }
    }
  } else {
    strcpy(file_type, "");
  }
  
  sprintf(filename, "bns_nurates_data_%s%s_point_%d_%d_bins.txt", reac_type, file_type, id_point, n_bins);

  strcpy(full_path, folder); 
  strcat(full_path, filename);
  
  FILE *fp = fopen(full_path, "w");
  if (fp == NULL) {
    printf("Error opening the file %s", filename);
    return;
  }

  fprintf(fp, "############################################################\n");
  fprintf(fp, "# bns_nurates vs WeakRates_stand_alone vs NuLib comparison #\n");
  fprintf(fp, "############################################################\n");
  fprintf(fp, "\n");


  fprintf(fp, "Generating quadratures with 2 * %d points ...\n", n_leg);
  MyQuadrature quad_1d = {.nx = n_leg, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&quad_1d);
  MyQuadrature quad_2d = {.nx = n_leg, .ny = n_leg, .dim = 2, .type = kGauleg, .x1 = 0., .x2 = 1., .y1 = 0., .y2 = 1.};
  GaussLegendreMultiD(&quad_2d);
  fprintf(fp, "Quadratures generated.\n");
  fprintf(fp, "\n");

  // Energy bins in MeV (from NuLib)
  const double e_bins[n_bins] =  {
    5.00000000e-01,
    1.50000000e+00,
    2.51335229e+00,
    3.55376572e+00,
    4.62196295e+00,
    5.71868593e+00,
    6.84469644e+00,
    8.00077657e+00,
    9.18772934e+00,
    1.04063792e+01,
    1.16575725e+01,
    1.29421785e+01,
    1.42610893e+01,
    1.56152211e+01,
    1.70055143e+01,
    1.84329348e+01,
    1.98984739e+01,
    2.14031496e+01,
    2.29480071e+01,
    2.45341193e+01,
    2.61625880e+01,
    2.78345442e+01,
    2.95511493e+01,
    3.13135956e+01,
    3.31231073e+01,
    3.49809413e+01,
    3.68883879e+01,
    3.88467721e+01,
    4.08574541e+01,
    4.29218304e+01,
    4.50413351e+01,
    4.72174403e+01,
    4.94516574e+01,
    5.17455383e+01,
    5.41006764e+01,
    5.65187074e+01,
    5.90013109e+01,
    6.15502113e+01,
    6.41671790e+01,
    6.68540317e+01,
    6.96126357e+01,
    7.24449070e+01,
    7.53528129e+01,
    7.83383732e+01,
    8.14036617e+01,
    8.45508073e+01,
    8.77819961e+01,
    9.10994725e+01,
    9.45055407e+01,
    9.80025665e+01,
    1.01592979e+02,
    1.05279272e+02,
    1.09064005e+02,
    1.12949809e+02,
    1.16939381e+02,
    1.21035493e+02,
    1.25240990e+02,
    1.29558793e+02,
    1.33991901e+02,
    1.38543393e+02,
    1.43216431e+02,
    1.48014260e+02,
    1.52940213e+02,
    1.57997712e+02,
    1.63190270e+02,
    1.68521492e+02,
    1.73995082e+02,
    1.79614843e+02,
    1.85384676e+02,
    1.91308591e+02,
    1.97390701e+02,
    2.03635231e+02,
    2.10046519e+02,
    2.16629017e+02,
    2.23387299e+02,
    2.30326057e+02,
    2.37450112e+02,
    2.44764412e+02,
    2.52274037e+02,
    2.59984204e+02,
    2.67900267e+02,
    2.76027725e+02,
    2.84372224e+02,
    2.92939559e+02,
    3.01735681e+02,
    3.10766700e+02,
    3.20038888e+02,
    3.29558686e+02,
    3.39332706e+02,
    3.49367737e+02,
    3.59670750e+02,
    3.70248900e+02,
    3.81109535e+02,
    3.92260199e+02,
    4.03708637e+02,
    4.15462800e+02,
    4.27530853e+02,
    4.39921178e+02,
    4.52642382e+02,
    4.65703300e+02,
    4.79113004e+02,
    4.92880809e+02,
    5.07016277e+02,
    5.21529227e+02,
    5.36429739e+02,
    5.51728163e+02,
    5.67435125e+02,
    5.83561535e+02,
    6.00118593e+02,
    6.17117801e+02,
    6.34570965e+02,
    6.52490209e+02,
    6.70887978e+02,
    6.89777052e+02,
    7.09170551e+02,
    7.29081945e+02,
    7.49525064e+02,
    7.70514109e+02,
    7.92063656e+02,
    8.14188675e+02,
    8.36904534e+02,
    8.60227009e+02,
    8.84172302e+02,
    9.08757043e+02,
    9.33998310e+02,
    9.59913634e+02,
    9.86521015e+02,
    1.01383894e+03
  };

    
  // Print opacity flags
  fprintf(fp, "List of included reactions:\n");
  fprintf(fp, "Neutrino abs on nucleons - e+/e- capture    : %d\n", opacity_flags->use_abs_em);
  fprintf(fp, "Elastic scattering on nucleons              : %d\n", opacity_flags->use_iso);
  fprintf(fp, "Nucleon-nucleon bremsstrahlung (and inverse): %d\n", opacity_flags->use_brem);
  fprintf(fp, "e+ e- annihilation (and inverse)            : %d\n", opacity_flags->use_pair);
  fprintf(fp, "Inelastic scattering on e+/e-               : %d\n", opacity_flags->use_inelastic_scatt);
  fprintf(fp, "\n");

   // Print properties of thermodynamic point
  fprintf(fp, "Thermodynamic point for comparison\n");
  fprintf(fp, "rho  = %lf g/cm^3\n", rho);
  fprintf(fp, "nb   = %lf 1/cm^3\n", eos_pars->nb);
  fprintf(fp, "T    = %lf MeV   \n", eos_pars->temp);
  fprintf(fp, "Ye   = %lf       \n", eos_pars->ye);
  fprintf(fp, "Yn   = %lf       \n", eos_pars->yn);
  fprintf(fp, "Yp   = %lf       \n", eos_pars->yp);
  fprintf(fp, "Abar = %lf       \n", Abar);
  fprintf(fp, "Zbar = %lf       \n", Zbar);
  fprintf(fp, "mu_e = %lf MeV   \n", eos_pars->mu_e);
  fprintf(fp, "mu_n = %lf MeV   \n", eos_pars->mu_n);
  fprintf(fp, "mu_p = %lf MeV   \n", eos_pars->mu_p);
  fprintf(fp, "mu_hat = %lf MeV \n", mu_hat);
  fprintf(fp, "\n");

  // Neutrino distribution function parameters
  NuDistributionParams distr_pars = NuEquilibriumParams(eos_pars); // Fermi-Dirac equilibrium distribution function

  // M1 quantities
  M1Quantities m1_pars;
  MyQuadratureIntegrand Jvals = NuEnergy(&distr_pars);
  MyQuadratureIntegrand nvals = NuNumber(&distr_pars);
  for (int idx = 0; idx < total_num_species; idx++) {
    m1_pars.n[idx] = nvals.integrand[idx];
    m1_pars.J[idx] = Jvals.integrand[idx];
  }
  
  // Print M1 quantities
  fprintf(fp, "Neutrino number and energy density:\n");
  fprintf(fp, "n_nue  = %.8e 1/cm^3, J_nue  = %.8e MeV/cm^3\n", m1_pars.n[id_nue] , m1_pars.J[id_nue] );
  fprintf(fp, "n_anue = %.8e 1/cm^3, J_anue = %.8e MeV/cm^3\n", m1_pars.n[id_anue], m1_pars.J[id_anue]);
  fprintf(fp, "n_nux  = %.8e 1/cm^3, J_nux  = %.8e MeV/cm^3\n", m1_pars.n[id_nux] , m1_pars.J[id_nux] );
  fprintf(fp, "n_anux = %.8e 1/cm^3, J_anux = %.8e MeV/cm^3\n", m1_pars.n[id_anux], m1_pars.J[id_anux] );
  fprintf(fp, "\n");

  // Grey opacity parameters
  GreyOpacityParams grey_pars = {.eos_pars = *eos_pars, .opacity_flags = *opacity_flags, .opacity_pars = opacity_params_default_none, .distr_pars = distr_pars, .m1_pars =  m1_pars};
  grey_pars.opacity_pars.neglect_blocking = neglect_blocking;
  
  // Compute spectral emissivities and opacities for NuLib energies
  fprintf(fp, "# Energy [MeV], j_nue/anue/nux/anux [s-1], NON STIMULATED op_ab_nue/anue/nux/anux [cm-1], STIMULATED op_ab_nue/anue/nux/anux [cm-1], op_sc_nue/anue/nux/anux [cm-1]\n");
  for (int i = 0; i < n_bins; i++) {
    SpectralOpacities spec_opacities = ComputeSpectralOpacitiesNotStimulatedAbs(e_bins[i], &quad_1d, &grey_pars);
    fprintf(fp, "%lf ", e_bins[i]);
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fp, "%.10e ", spec_opacities.j[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fp, "%.10e ", spec_opacities.kappa[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fp, "%.10e ", spec_opacities.j[idx] / kClight + spec_opacities.kappa[idx]);
    }
    for (int idx = 0; idx < total_num_species; idx++) {
      fprintf(fp, "%.10e ", spec_opacities.kappa_s[idx]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  // Compute M1 source coefficients
  M1Opacities coeffs = ComputeM1Opacities(&quad_1d, &quad_2d, &grey_pars);
  M1Opacities stim_coeffs = ComputeM1OpacitiesNotStimulated(&quad_1d, &quad_2d, &grey_pars);
  
  // Print M1 source coefficients
  fprintf(fp, "Emission coefficients (no stimulated absorption):\n");
  fprintf(fp, "eta_0_nue  = %.5e cm-3 s-1    \n", stim_coeffs.eta_0[id_nue]);
  fprintf(fp, "eta_0_anue = %.5e cm-3 s-1    \n", stim_coeffs.eta_0[id_anue]);
  fprintf(fp, "eta_0_nux  = %.5e cm-3 s-1    \n", stim_coeffs.eta_0[id_nux]);
  fprintf(fp, "eta_0_anux = %.5e cm-3 s-1    \n", stim_coeffs.eta_0[id_anux]);
  fprintf(fp, "eta_nue    = %.5e MeV cm-3 s-1\n", stim_coeffs.eta[id_nue]);
  fprintf(fp, "eta_anue   = %.5e MeV cm-3 s-1\n", stim_coeffs.eta[id_anue]);
  fprintf(fp, "eta_nux    = %.5e MeV cm-3 s-1\n", stim_coeffs.eta[id_nux]);
  fprintf(fp, "eta_anux   = %.5e MeV cm-3 s-1\n", stim_coeffs.eta[id_anux]);
  fprintf(fp, "\n");

  fprintf(fp, "Emission coefficients (stimulated absorption):\n");
  fprintf(fp, "eta_0_nue  = %.5e cm-3 s-1    \n", coeffs.eta_0[id_nue]);
  fprintf(fp, "eta_0_anue = %.5e cm-3 s-1    \n", coeffs.eta_0[id_anue]);
  fprintf(fp, "eta_0_nux  = %.5e cm-3 s-1    \n", coeffs.eta_0[id_nux]);
  fprintf(fp, "eta_0_anux = %.5e cm-3 s-1    \n", coeffs.eta_0[id_anux]);
  fprintf(fp, "eta_nue    = %.5e MeV cm-3 s-1\n", coeffs.eta[id_nue]);
  fprintf(fp, "eta_anue   = %.5e MeV cm-3 s-1\n", coeffs.eta[id_anue]);
  fprintf(fp, "eta_nux    = %.5e MeV cm-3 s-1\n", coeffs.eta[id_nux]);
  fprintf(fp, "eta_anux   = %.5e MeV cm-3 s-1\n", coeffs.eta[id_anux]);
  fprintf(fp, "\n");

  fprintf(fp, "Absorption coefficients (no stimulated absorption):\n");
  fprintf(fp, "kappa_0_a_nue  = %.5e cm-1\n", stim_coeffs.kappa_0_a[id_nue]);
  fprintf(fp, "kappa_0_a_anue = %.5e cm-1\n", stim_coeffs.kappa_0_a[id_anue]);
  fprintf(fp, "kappa_0_a_nux  = %.5e cm-1\n", stim_coeffs.kappa_0_a[id_nux]);
  fprintf(fp, "kappa_0_a_anux = %.5e cm-1\n", stim_coeffs.kappa_0_a[id_anux]);
  fprintf(fp, "kappa_a_nue    = %.5e cm-1\n", stim_coeffs.kappa_a[id_nue]);
  fprintf(fp, "kappa_a_anue   = %.5e cm-1\n", stim_coeffs.kappa_a[id_anue]);
  fprintf(fp, "kappa_a_nux    = %.5e cm-1\n", stim_coeffs.kappa_a[id_nux]);
  fprintf(fp, "kappa_a_anux   = %.5e cm-1\n", stim_coeffs.kappa_a[id_anux]);
  fprintf(fp, "\n");

  fprintf(fp, "Absorption coefficients (stimulated absorption):\n");
  fprintf(fp, "kappa_0_a_nue  = %.5e cm-1\n", coeffs.kappa_0_a[id_nue]);
  fprintf(fp, "kappa_0_a_anue = %.5e cm-1\n", coeffs.kappa_0_a[id_anue]);
  fprintf(fp, "kappa_0_a_nux  = %.5e cm-1\n", coeffs.kappa_0_a[id_nux]);
  fprintf(fp, "kappa_0_a_anux = %.5e cm-1\n", coeffs.kappa_0_a[id_anux]);
  fprintf(fp, "kappa_a_nue    = %.5e cm-1\n", coeffs.kappa_a[id_nue]);
  fprintf(fp, "kappa_a_anue   = %.5e cm-1\n", coeffs.kappa_a[id_anue]);
  fprintf(fp, "kappa_a_nux    = %.5e cm-1\n", coeffs.kappa_a[id_nux]);
  fprintf(fp, "kappa_a_anux   = %.5e cm-1\n", coeffs.kappa_a[id_anux]);
  fprintf(fp, "\n");

  fprintf(fp, "Scattering coefficients:\n");
  fprintf(fp, "kappa_s_nue  = %.5e cm-1\n", coeffs.kappa_s[id_nue]);
  fprintf(fp, "kappa_s_anue = %.5e cm-1\n", coeffs.kappa_s[id_anue]);
  fprintf(fp, "kappa_s_nux  = %.5e cm-1\n", coeffs.kappa_s[id_nux]);
  fprintf(fp, "kappa_s_anux = %.5e cm-1\n", coeffs.kappa_s[id_anux]);
  fprintf(fp, "\n");

  fclose(fp);

  return;
}

int main() {

  // Thermodynamic points for comparison
  const double rho_array[n_points]   = {1435091918127104.0, 993097437675520.00, 511017121480704.00, 23989367865344.000}; // [g cm-3]
  const double temp_array[n_points]  = {51.150005340576172, 43.034099578857422, 31.396757125854492, 5.8661994934082031}; // [MeV]
  const double ye_array[n_points]    = {0.21549025177955627, 0.18022947013378143, 0.12433240562677383, 3.5000000149011612E-02};
  const double xn_array[n_points]    = {0.78450974821565389, 0.81977052986906962, 0.87566759439393282, 0.90764789847055527};
  const double xp_array[n_points]    = {0.21549025178434600, 0.18022947013093041, 0.12433240560606711, 7.3319472295880960E-03};
  const double abar_array[n_points]  = {1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 7.6262107448084606};
  const double zbar_array[n_points]  = {1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 2.0767784835527388};
  const double mue_array[n_points]   = {324.00790839794519, 269.64449212325241, 189.96101970245701, 46.198652099231772};
  const double mun_array[n_points]   = {585.93269122709501, 306.62799804149404, 87.380083431825298, 8.6121695893764922};
  const double mup_array[n_points]   = {419.35214965870006, 148.84071174828097, -65.656622311432869, -32.763517246657820};
  const double muhat_array[n_points] = {166.58054156839492, 157.78728629321293, 153.03670574325818, 41.375686836034305};

  // Reference baryonic mass
  const double m_ref =  931.494061; // [MeV]
  double nb;

  // Opacity flags
  OpacityFlags opacity_flags = {.use_abs_em = 0, .use_iso = 1, .use_brem = 0, .use_pair = 0, .use_inelastic_scatt = 0};
  const char* reac_type = "iso";

  // Number of quadrature points ==> 2 * n
  const int n_leg = 20;
  
  for (int id_p = 1; id_p <= 4; id_p++){

    nb = rho_array[id_p-1] * kClight * kClight / (m_ref * kMeV); // [cm^-3]

    // EOS parameters
    MyEOSParams eos_pars = {.nb = nb,
                            .temp = temp_array[id_p-1],
                            .ye = ye_array[id_p-1],
                            .yn = xn_array[id_p-1],
                            .yp = xp_array[id_p-1],
                            .mu_e = mue_array[id_p-1],
                            .mu_n = mun_array[id_p-1],
                            .mu_p = mup_array[id_p-1]};
  
    if (strcmp(reac_type, "pair") || strcmp(reac_type, "brem") || strcmp(reac_type, "all")) {
      save_rates_to_file(id_p, reac_type, n_leg, &eos_pars, &opacity_flags, rho_array[id_p-1], abar_array[id_p-1], zbar_array[id_p-1], muhat_array[id_p-1], false);
      save_rates_to_file(id_p, reac_type, n_leg, &eos_pars, &opacity_flags, rho_array[id_p-1], abar_array[id_p-1], zbar_array[id_p-1], muhat_array[id_p-1], true);
    } else {
      save_rates_to_file(id_p, reac_type, n_leg, &eos_pars, &opacity_flags, rho_array[id_p-1], abar_array[id_p-1], zbar_array[id_p-1], muhat_array[id_p-1], false);
    }

  }

  return 0;
}

/*
  // Thermodynamic point for comparison
  const double rho  = 1435091918127104.0;    // [g/cm^3]
  const double temp = 51.150005340576172;     // [MeV]
  const double ye   = 0.21549025177955627;
  const double xn   = 0.78450974821565389;
  const double xp   = 0.21549025178434600;
  const double Abar = 1.0000000000000000;
  const double Zbar = 1.0000000000000000;
  const double mu_e = 324.00790839794519;    // [MeV]
  const double mu_n = 585.93269122709501;    // [MeV]
  const double mu_p = 419.35214965870006;    // [MeV]
  const double mu_hat = 166.58054156839492;  // [MeV]


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
  printf("nux:  R_in = %.5e, R_out = %.5e\n", nes_kernels.em[id_nux] , nes_kernels.abs[id_nux]);
  printf("anux: R_in = %.5e, R_out = %.5e\n", nes_kernels.em[id_anux], nes_kernels.abs[id_anux]);
  printf("\n");

  printf("NPS kernels [cm^3 s^-1]:\n");

  printf("nue:  R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_nue] , nps_kernels.abs[id_nue]);
  printf("anue: R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_anue], nps_kernels.abs[id_anue]);
  printf("nux:  R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_nux],  nps_kernels.abs[id_nux]);
  printf("anux: R_in = %.5e, R_out = %.5e\n", nps_kernels.em[id_anux], nps_kernels.abs[id_anux]);
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
  double e_nu = 110.; // [MeV]

  // M1 Grey Coefficient parameters  
  GreyOpacityParams grey_pars = {.distr_pars = distr_pars,
                                 .eos_pars = eos_pars,
                                 .opacity_pars = opacity_pars,
                                 .kernel_pars = kernel_pars,
                                 .opacity_flags = opacity_flags};

  // Compute emissivity and inverse mean free path
  SpectralOpacities inelastic_output = ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, &my_quadrature_1d, &grey_pars);

  // Print emissivity and inverse mean free path
  printf("Emissivity [s-1] and inverse mean free path [cm-1] (E = %.2lf MeV):\n", e_nu);

  printf("nue:  eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_nue] , inelastic_output.kappa[id_nue]);
  printf("anue: eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_anue], inelastic_output.kappa[id_anue]);
  printf("nux:  eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_nux] , inelastic_output.kappa[id_nux]);
  printf("anux: eta = %.5e, kappa = %.5e\n", inelastic_output.j[id_anux], inelastic_output.kappa[id_anux]);
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
  printf("r j-nue j-anue j-nux j-anux kappa-nue kappa-anue kappa-nux kappa-anux\n");

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
    
    inelastic_output = ComputeSpectralOpacitiesNotStimulatedAbs(e_nu, &my_quadrature_1d, &grey_pars);

    //printf("%e %e %e %e %e %e %e %e %e\n",
    //        r[i], 
    //        inelastic_output.j[id_nue], inelastic_output.j[id_anue], inelastic_output.j[id_nux], inelastic_output.j[id_anux],
    //        inelastic_output.kappa[id_nue], inelastic_output.kappa[id_anue], inelastic_output.kappa[id_nux], inelastic_output.kappa[id_anux]
    //      );
    
  }

  static const int n = 10000;

  double alpha_min = 0.;
  double alpha_max = 700.;

  double pair_t_sum[6] = {0.};
  double pair_t_fit[6] = {0.};
  double rel_diff[6] = {0.};

  double dx = (alpha_max - alpha_min) / (double) n;
  double alpha;

  for (int idx=0; idx<n; idx++) {
    alpha = alpha_min + idx * dx;

    for (int j=0; j<6; j++) {
      pair_t_sum[j] = PairT(j+1, alpha, 1.0E-06);
      pair_t_fit[j] = PairTFitted(j+1, alpha);
      rel_diff[j] = fabs(pair_t_sum[j] - pair_t_fit[j]) / pair_t_sum[j];
    }
    //printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", alpha, rel_diff[0], rel_diff[1], rel_diff[2], rel_diff[3], rel_diff[4], rel_diff[5]);
  }

  //for (int l=2; l<=6; l++) {
    //printf("l = %d: zeta = %.30e\n", l, gsl_sf_zeta_int(l));
  //}
  
  double x1 = 10.;
  double x2 = 20.;

  PairKernelParams pair_pars = {.omega = x1, .omega_prime=x2};


  printf("Pair Kernel: %.5e\n", PairKernelsM1Optimized(&eos_pars,&pair_pars).abs_e);
  
  printf("%.5e\n", FDI_p2(15.67));
 





  MyQuadrature quad_GLQ_64 = {.dim = 1, .nx = 64, .alpha = 0};
  GaussLaguerre(&quad_GLQ_64);

  MyKernelOutput pair_out;

  double e_prime;

  PairKernelParams pair_pars = grey_pars.kernel_pars.pair_kernel_params;
  pair_pars.omega = e_bins[0];

  for (int i = 0; i < quad_GLQ_64.nx; i++) {
  //for (int i = 0; i < 1; i++) {
    //printf("%.10lf  %.10lf\n", quad_GLQ_64.points[i], quad_GLQ_64.w[i]);
    e_prime = quad_GLQ_64.points[i];
    pair_pars.omega_prime = e_prime * eos_pars.temp;
    pair_out = PairKernelsOptimized(&eos_pars, &pair_pars);
    // printf("%.10lf  %.10e\n", e_prime, pair_out.em[id_nue]); //, pair_out.em[id_anue], pair_out.em[id_nux], pair_out.em[id_anux]);
    //printf("%.10e ", pair_out.em[id_nue]);
    printf("%.10e ", e_prime);
    //exit(0);
  }





  const double e_array[18] = {1.        ,   3.        ,   5.23824436,   8.00973783,
        11.44152399,  15.69091386,  20.95269691,  27.46807011,
        35.53569425,  45.52538436,  57.89506183,  73.21174523,
        92.17754211, 115.66183319, 144.74112423, 180.74839245,
       225.33418934, 280.54230101}; // [MeV]




  const double eta_array[61] = {0.1       ,   0.11220185,   0.12589254,   0.14125375,
         0.15848932,   0.17782794,   0.19952623,   0.22387211,
         0.25118864,   0.28183829,   0.31622777,   0.35481339,
         0.39810717,   0.44668359,   0.50118723,   0.56234133,
         0.63095734,   0.70794578,   0.79432823,   0.89125094,
         1.        ,   1.12201845,   1.25892541,   1.41253754,
         1.58489319,   1.77827941,   1.99526231,   2.23872114,
         2.51188643,   2.81838293,   3.16227766,   3.54813389,
         3.98107171,   4.46683592,   5.01187234,   5.62341325,
         6.30957344,   7.07945784,   7.94328235,   8.91250938,
        10.        ,  11.22018454,  12.58925412,  14.12537545,
        15.84893192,  17.7827941 ,  19.95262315,  22.38721139,
        25.11886432,  28.18382931,  31.6227766 ,  35.48133892,
        39.81071706,  44.66835922,  50.11872336,  56.23413252,
        63.09573445,  70.79457844,  79.43282347,  89.12509381,
       100.};

  const double T_array[65] = {5.00000000e-02, 5.66630600e-02, 6.42140473e-02, 7.27712882e-02,
       8.24688774e-02, 9.34587789e-02, 1.05913208e-01, 1.20027329e-01,
       1.36022315e-01, 1.54148812e-01, 1.74690867e-01, 1.97970381e-01,
       2.24352152e-01, 2.54249589e-01, 2.88131194e-01, 3.26527902e-01,
       3.70041402e-01, 4.19353563e-01, 4.75237122e-01, 5.38567791e-01,
       6.10337981e-01, 6.91672352e-01, 7.83845439e-01, 8.88301622e-01,
       1.00667776e+00, 1.14082885e+00, 1.29285707e+00, 1.46514475e+00,
       1.66039170e+00, 1.88165749e+00, 2.13240942e+00, 2.41657686e+00,
       2.73861279e+00, 3.10356361e+00, 3.51714822e+00, 3.98584761e+00,
       4.51700644e+00, 5.11894814e+00, 5.80110531e+00, 6.57416756e+00,
       7.45024901e+00, 8.44307813e+00, 9.56821285e+00, 1.08432844e+01,
       1.22882734e+01, 1.39258235e+01, 1.57815954e+01, 1.78846698e+01,
       2.02680023e+01, 2.29689406e+01, 2.60298092e+01, 2.94985728e+01,
       3.34295880e+01, 3.78844550e+01, 4.29329829e+01, 4.86542836e+01,
       5.51380118e+01, 6.24857694e+01, 7.08126980e+01, 8.02492830e+01,
       9.09433987e+01, 1.03062625e+02, 1.16796874e+02, 1.32361366e+02,
       1.50000000e+02};

  double table_out[65][61][18][4][18][2] = {0};

  MyKernelOutput pair_out;

  PairKernelParams pair_pars = grey_pars.kernel_pars.pair_kernel_params;

  //printf("[");
  for (int i = 0; i < 65; i++) {
    eos_pars.temp = T_array[i];
    //printf("[");
    for (int j = 0; j < 61; j++) {
      eos_pars.mu_e = eos_pars.temp * eta_array[j];
      //printf("[");
      for (int k = 0; k < 18; k++) {
        pair_pars.omega = e_array[k];
        //printf("[");
        for (int l = 0; l < 4; l++) {
          //printf("[");
          for (int m = 0; m < 18; m++) {
            pair_pars.omega_prime = e_array[m];
            pair_out = PairKernelsOptimized(&eos_pars, &pair_pars);
            table_out[i][j][k][l][m][0] = pair_out.em[l];
            table_out[i][j][k][l][m][1] = pair_out.abs[l];
            //printf("[%.10e, %.10e]", pair_out.em[l], pair_out.abs[l]);
            //if (m != 17) { printf(", "); };
          }
          //printf("]");
          //if (l != 3) { printf(", "); };
        }
        //printf("]");
        //if (k != 17) { 
          //printf(", ");
        //} else {
          //printf("\n");
        //}
      }
      //printf("]");
      //if (j != 60) { 
        //printf(", ");
      //} else {
        //printf("\n");
      //}
    }
    //printf("]");
    //if (i != 64) {
      //printf(", ");
    //} else {
      //printf("\n");
    //}
  }
  //printf("]");
*/

 /* now the various steps involved in preparing an hdf5 file */
  //hid_t file_id;
  /* at this time our basic file just has the values strung out,
     so the hdf5 rank is 1 */
  //hsize_t dims[6] = {65, 61, 18, 4, 18, 2};
  //herr_t status;

  /* create a HDF5 file */
  //file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  /* create and write a double type dataset named "/datacube" */
  //status = H5LTmake_dataset(file_id, "/datacube", 3, dims,
                           //H5T_NATIVE_DOUBLE, &table_out);
  //assert(status != -1);
  /* add some hdf5 attributes with the metadata we need */
  //H5LTset_attribute_int(file_id, "/datacube", "nx", &nx, 1);
  //H5LTset_attribute_int(file_id, "/datacube", "nx", &nx, 1);
  //H5LTset_attribute_int(file_id, "/datacube", "ny", &ny, 1);
  //H5LTset_attribute_string(file_id, "/datacube", "x_units", "m");
  //H5LTset_attribute_string(file_id, "/datacube", "y_units", "kg");
  //H5LTset_attribute_string(file_id, "/datacube", "z_units", "degK");

  //status = H5Fclose (file_id); assert(status != -1);
  //assert(status != -1);



/*
  static const int n_alpha = 10000;

  double alpha_min = 1.0E-05;
  double alpha_max = 30.;

  double dalpha = (log10(alpha_max) - log10(alpha_min)) / (double) n_alpha;
  double alpha;

  for (int idx=0; idx<n_alpha; idx++) {
    alpha = pow(10., log10(alpha_min) + idx * dalpha);
    printf("%.15e ", alpha);
    for (int j=0; j<6; j++) {
      printf("%.15e ", PairT(j+1, alpha, 1.0E-15));
    }
    printf("\n");
  }


  alpha_min = 30.;
  alpha_max = 750.;

  dalpha = (alpha_max - alpha_min) / (double) n_alpha;

  for (int idx=0; idx<n_alpha; idx++) {
    alpha = alpha_min + idx * dalpha;
    printf("%.15e ", alpha);
    for (int j=0; j<6; j++) {
      printf("%.15e ", PairT(j+1, alpha, 1.0E-15));
    }
    printf("\n");
  }


  SpectralOpacities brem_out;
  double omega;

  for (int i = 0; i < 18; i++) {
    omega = e_array[i];
    brem_out = ComputeSpectralOpacitiesNotStimulatedAbs(omega, &quad_1d, &grey_pars);
    printf("%.15e, ", brem_out.j[0]);   
  }
*/

