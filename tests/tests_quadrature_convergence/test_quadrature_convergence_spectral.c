//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_quadrature_convergence_spectral.c
//  \brief Test convergence of spectral rates with number of quadrature points

#include <stdio.h>
#include <string.h>

#include "bns_nurates.h"
#include "constants.h"
#include "opacities.h"
#include "integration.h"
#include "distribution.h"


#define n_bins 128
#define n_points 4

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

void ComputeSpectralRatesGivenQuadrature(const int n_leg, const int id_point, const char *reac_type, MyEOSParams *eos_pars, OpacityFlags *opacity_flags) {

  const char folder[300] = "/home/leonardo/Desktop/PhD_work/BNS_muons/bns_refactor/bns_nurates/tests/tests_quadrature_convergence/output/spectral/";
  char filename[100], full_path[300];

  sprintf(filename, "bns_nurates_data_%s_point_%d_quad_%d.txt", reac_type, id_point, 2 * n_leg);
  printf("%s\n",filename);
  strcpy(full_path, folder); 
  strcat(full_path, filename);
  
  FILE *fp = fopen(full_path, "w");
  if (fp == NULL) {
    printf("Error opening the file %s\n", filename);
    return;
  }

  printf("Computing spectral rates for n_leg = 2 * %d....\n", n_leg);

  MyQuadrature quad_1d = {.nx = n_leg, .dim = 1, .type = kGauleg, .x1 = 0., .x2 = 1.};
  GaussLegendreMultiD(&quad_1d);

  fprintf(fp, "# Convergence test of spectral rates with number of quadrature points: n = 2 * %d\n", n_leg);
  fprintf(fp, "#\n");

  // Print opacity flags
  fprintf(fp, "# List of included reactions:\n");
  fprintf(fp, "# Neutrino abs on nucleons - e+/e- capture    : %d\n", opacity_flags->use_abs_em);
  fprintf(fp, "# Elastic scattering on nucleons              : %d\n", opacity_flags->use_iso);
  fprintf(fp, "# Nucleon-nucleon bremsstrahlung (and inverse): %d\n", opacity_flags->use_brem);
  fprintf(fp, "# e+ e- annihilation (and inverse)            : %d\n", opacity_flags->use_pair);
  fprintf(fp, "# Inelastic scattering on e+/e-               : %d\n", opacity_flags->use_inelastic_scatt);
  fprintf(fp, "#\n");

   // Print properties of thermodynamic point
  fprintf(fp, "# Thermodynamic point for comparison\n");
  fprintf(fp, "# nb   = %lf 1/cm^3\n", eos_pars->nb);
  fprintf(fp, "# T    = %lf MeV   \n", eos_pars->temp);
  fprintf(fp, "# Ye   = %lf       \n", eos_pars->ye);
  fprintf(fp, "# Yn   = %lf       \n", eos_pars->yn);
  fprintf(fp, "# Yp   = %lf       \n", eos_pars->yp);
  fprintf(fp, "# mu_e = %lf MeV   \n", eos_pars->mu_e);
  fprintf(fp, "# mu_n = %lf MeV   \n", eos_pars->mu_n);
  fprintf(fp, "# mu_p = %lf MeV   \n", eos_pars->mu_p);
  fprintf(fp, "#\n");

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
  
  // Grey opacity parameters
  GreyOpacityParams grey_pars = {.eos_pars = *eos_pars, .opacity_flags = *opacity_flags, .opacity_pars = opacity_params_default_none, .distr_pars = distr_pars, .m1_pars =  m1_pars};
  
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

  fclose(fp);
  
  printf("Done!\n");
  printf("\n");

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

  // Number of quadrature points
  int n_leg[6] = {5, 10, 20, 30, 40, 50};

  // Opacity flags
  OpacityFlags opacity_flags = {.use_abs_em = 1, .use_iso = 1, .use_brem = 1, .use_pair = 1, .use_inelastic_scatt = 0};
  const char* reac_type = "noinel";
 
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
  
    for (int i = 0; i < 6; i++) {
      ComputeSpectralRatesGivenQuadrature(n_leg[i], id_p, reac_type, &eos_pars, &opacity_flags);
    }
    
  }

  return 0;
}