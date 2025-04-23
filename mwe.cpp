#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bns_nurates.hpp"
#include "m1_opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

int main(int argc, char* argv[])
{
    // Fix the neutrino energy for computation of spectral rates
    const double nu_energy = 10.; // [MeV]
    
    // Input thermodynamic quantities (corresponding to point A in Chiesa+25 PRD)
    // N.B.: chemical potentials include the rest mass contribution
    const double nb = 4.208366627847035e+38;  // Baryon number density [cm-3]           
    const double T = 12.406403541564941;      // Temperature [MeV]
    const double ye = 0.07158458232879639;    // Electron fraction
    const double mu_e = 187.1814489;          // Electron chemical potential [MeV]
    const double mu_p = 1011.01797737;        // Proton chemical potential [MeV]
    const double mu_n = 1221.59013681;        // Neutron chemical potential [MeV]
    const double dU = 18.92714728;            // Nucleon interaction potential difference (Un-Up) [MeV]
    const double mp_eff = 278.87162217;       // Proton effective mass [MeV]
    const double mn_eff = 280.16495513;       // Neutron effective mass [MeV]
    const double dm = mn_eff - mp_eff;        // Nucleon effective mass difference [MeV]

    // Input gray neutrino quantities (library supports 4 neutrino species)
    // N.B.: These will be used in PART 2 for reconstructing the neutrino distribution functions
    //       and as normalization factors for energy-averaged opacities
    const double n_nue  = 3.739749408027436e+33;   // Electron neutrino number density [cm-3]
    const double n_anue = 1.2174961961689319e+35;  // Electron antineutrino number density [cm-3]
    const double n_nux  = 2.2438496448164613e+34;  // Heavy-type neutrino number density [cm-3]
    const double n_anux = 2.2438496448164613e+34;  // Heavy-type antineutrino number density [cm-3]
    const double j_nue  = 1.246583136009145e+35;   // Electron neutrino energy density [MeV cm-3]
    const double j_anue = 5.360307484839323e+36;   // Electron antineutrino number density [MeV cm-3]
    const double j_nux  = 8.726081952064015e+35;   // Heavy-type neutrino energy density [MeV cm-3]
    const double j_anux = 8.726081952064015e+35;   // Heavy-type antineutrino energy density [MeV cm-3]
    const double chi_nue   = 1. / 3.;              // Electron neutrino Eddington factor
    const double chi_anue  = 1. / 3.;              // Electron antineutrino Eddington factor
    const double chi_nux   = 1. / 3.;              // Heavy-type neutrino Eddington factor
    const double chi_anux  = 1. / 3.;              // Heavy-type antineutrino Eddington factor

    // Create a quadrature struct and populate it with data relative to
    // Gauss-Legendre quadrature
    // N.B.: only the 'nx' member can be modified, all the others should
    //       be always set as the following
    MyQuadrature my_quadrature;
    my_quadrature.nx   = 6;  // number of quadrature points   
    my_quadrature.dim  = 1;
    my_quadrature.type = kGauleg;
    my_quadrature.x1   = 0.;
    my_quadrature.x2   = 1.;
    GaussLegendre(&my_quadrature);

    // Create two structs that will hold the final result of the
    // computation of spectral and gray rates, respectively
    SpectralOpacities spectral_rates;
    M1Opacities gray_rates;

    // Create an opacity params structure, to activate/deactivate specific
    // reactions or corrections and pass physical parameters
    GreyOpacityParams my_grey_opacity_params = {0};
    my_grey_opacity_params.opacity_flags     = opacity_flags_default_none;
    my_grey_opacity_params.opacity_pars      = opacity_params_default_none;

    // Select active reactions
    my_grey_opacity_params.opacity_flags.use_abs_em =
        1; // Activate beta processes
    my_grey_opacity_params.opacity_flags.use_brem =
        1; // Activate Bremsstrahlung
    my_grey_opacity_params.opacity_flags.use_pair =
        1; // Activate pair processes
    my_grey_opacity_params.opacity_flags.use_iso =
        1; // Activate scattering on nucleons
    my_grey_opacity_params.opacity_flags.use_inelastic_scatt =
        1; // Activate scattering on leptons

    // Select corrections to rates
    my_grey_opacity_params.opacity_pars.use_dm_eff =
        0; // Do not use effective mass correction to beta processes
    my_grey_opacity_params.opacity_pars.use_dU =
        1; // Use effective potential correction to beta processes
    my_grey_opacity_params.opacity_pars.use_decay =
        1; // Include (inverse) nucleon decays to beta processes
    my_grey_opacity_params.opacity_pars.use_WM_ab =
        1; // Activate beta-processes weak magnetism correction
    my_grey_opacity_params.opacity_pars.use_WM_sc =
        1; // Activate isoenergetic scattering weak magnetism correction
    my_grey_opacity_params.opacity_pars.use_NN_medium_corr =
        1; // Activate NN bremsstrahlung medium correction as in Fischer+16

    // Pass thermodynamic conditions
    my_grey_opacity_params.eos_pars.nb =
        nb * 1e-21; // Set baryon number density (in baryon/nm^3)
    my_grey_opacity_params.eos_pars.temp = T; // Set temperature (in MeV)
    my_grey_opacity_params.eos_pars.yp   = ye;  // Set proton fraction
    my_grey_opacity_params.eos_pars.yn   = 1. - ye; // Set neutron fraction
    my_grey_opacity_params.eos_pars.mu_e =
        mu_e; // Set electron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.mu_p =
        mu_p; // Set proton chemical potential (in MeV) (NOTE: reactions only
              // depend on the difference between mu_n and mu_p)
    my_grey_opacity_params.eos_pars.mu_n =
        mu_n; // Set neutron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.dU =
        dU; // Set effective potential difference (in MeV)
    my_grey_opacity_params.eos_pars.dm_eff =
        dm; // Set effective mass difference (in MeV)

    printf("\nInput parameters\n");
    printf("----------------\n");
    printf("Neutrino energy (for spectral rates)               : %13.6e (MeV)\n", nu_energy);
    printf("Baryon number density                              : %13.6e (cm^-3)\n", nb);
    printf("Temperature                                        : %13.6e (MeV)\n", T);
    printf("Electron fraction                                  : %13.3f\n", ye);
    printf("Relativistic electron chemical potential           : %13.6e (MeV)\n", mu_e);
    printf("Relativistic proton chemical potential             : %13.6e (MeV)\n", mu_p);
    printf("Relativistic neutron chemical potential            : %13.6e (MeV)\n", mu_n);
    printf("Effective nucleon mass difference (n - p)          : %13.6e (MeV)\n", dm);
    printf("Effective nucleon potential difference (n - p)     : %13.6e (MeV)\n", dU);
    printf("Electron neutrinos number density 'n'              : %13.6e (cm^-3)\n", n_nue);
    printf("Electron antineutrinos number density 'n'          : %13.6e (cm^-3)\n", n_anue);
    printf("Heavy-type neutrinos number density 'n'            : %13.6e (cm^-3)\n", n_nux);
    printf("Heavy-type antineutrinos number density 'n'        : %13.6e (cm^-3)\n", n_anux);
    printf("Electron neutrinos energy density 'J'              : %13.6e (MeV cm^-3)\n", j_nue);
    printf("Electron antineutrinos energy density 'J'          : %13.6e (MeV cm^-3)\n", j_anue);
    printf("Heavy-type neutrinos energy density 'J'            : %13.6e (MeV cm^-3)\n", j_nux);
    printf("Heavy-type antineutrinos energy density 'J'        : %13.6e (MeV cm^-3)\n", j_anux);
    printf("Electron neutrinos Eddington parameter 'chi'       : %13.11f\n", chi_nue);
    printf("Electron antineutrinos Eddington parameter 'chi'   : %13.11f\n", chi_anue);
    printf("Heavy-type neutrinos Eddington parameter 'chi'     : %13.11f\n", chi_nux);
    printf("Heavy-type antineutrinos Eddington parameter 'chi' : %13.11f\n\n", chi_anux);

    ////////////////////////////////////////////////////////////
    // PART 1: compute rates assuming equilibrium with matter //
    ////////////////////////////////////////////////////////////

    // Compute neutrino distribution parameters assuming equilibrium
    my_grey_opacity_params.distr_pars =
        NuEquilibriumParams(&my_grey_opacity_params.eos_pars);

    // Compute gray neutrino number and energy densities assuming equilibrium
    // N.B.: they are required for the normalization factor of energy-averaged opacities,
    //       they do not enter in the computation of spectral rates
    ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                         &my_grey_opacity_params.distr_pars,
                         &my_grey_opacity_params.m1_pars);
    my_grey_opacity_params.m1_pars.chi[id_nue]  = 0.333333333333333333333333333;
    my_grey_opacity_params.m1_pars.chi[id_anue] = 0.333333333333333333333333333;
    my_grey_opacity_params.m1_pars.chi[id_nux]  = 0.333333333333333333333333333;
    my_grey_opacity_params.m1_pars.chi[id_anux] = 0.333333333333333333333333333;

    // Restore to input units (from MeV nm^-3 to g s^-2 nm^-1)
    my_grey_opacity_params.m1_pars.J[id_nue] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_anue] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_nux] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_anux] *= kBS_MeV;

    // Restore units of cm^-3 and MeV cm^-3 for printing
    const double n_nue_eq = my_grey_opacity_params.m1_pars.n[id_nue] * 1e21;
    const double n_anue_eq = my_grey_opacity_params.m1_pars.n[id_anue] * 1e21;
    const double n_nux_eq = my_grey_opacity_params.m1_pars.n[id_nux] * 1e21;
    const double n_anux_eq = my_grey_opacity_params.m1_pars.n[id_anux] * 1e21;

    const double J_nue_eq = my_grey_opacity_params.m1_pars.J[id_nue] * 1e21 / kBS_MeV;
    const double J_anue_eq = my_grey_opacity_params.m1_pars.J[id_anue] * 1e21 / kBS_MeV;
    const double J_nux_eq = my_grey_opacity_params.m1_pars.J[id_nux] * 1e21 / kBS_MeV;
    const double J_anux_eq = my_grey_opacity_params.m1_pars.J[id_anux] * 1e21 / kBS_MeV;

    printf("Reconstructed neutrino densities assuming equilibrium\n");
    printf("-----------------------------------------------------\n");
    printf("    nue            anue            nux            anux\n");
    printf("n  %13.6e  %13.6e   %13.6e  %13.6e     (cm^-3)\n", n_nue_eq, n_anue_eq, n_nux_eq, n_anux_eq);
    printf("J  %13.6e  %13.6e   %13.6e  %13.6e (MeV cm^-3)\n", J_nue_eq, J_anue_eq, J_nux_eq, J_anux_eq);
    printf("chi %12.10f   %12.10f    %12.10f   %12.10f\n\n", my_grey_opacity_params.m1_pars.chi[id_nue], my_grey_opacity_params.m1_pars.chi[id_anue], my_grey_opacity_params.m1_pars.chi[id_nux], my_grey_opacity_params.m1_pars.chi[id_anux]);

    // Compute and output spectral emissivities and inverse mean free paths (not in the stimulated
    // absorption formalism)
    spectral_rates = ComputeSpectralOpacitiesNotStimulatedAbs(nu_energy, &my_quadrature, &my_grey_opacity_params);
    
    // The numerical factors restore usual units (see output)
    // N.B.: 'j' and 'kappa' are the spectral emissivities and inverse mean free paths for the
    //       sum of the active inelastic reactions, 'j_s' and 'kappa_s' are the equivalent for
    //       the elastic scattering off nucleons
    printf("Spectral rates assuming equilibrium\n");
    printf("------------------------------\n");
    printf("     j             j_s           kappa         kappa_s\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_nue], spectral_rates.j_s[id_nue],
           spectral_rates.kappa[id_nue] * 1e7, spectral_rates.kappa_s[id_nue] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_anue], spectral_rates.j_s[id_anue],
           spectral_rates.kappa[id_anue] * 1e7, spectral_rates.kappa_s[id_anue] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_nux], spectral_rates.j_s[id_nux],
           spectral_rates.kappa[id_nux] * 1e7, spectral_rates.kappa_s[id_nux] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e\n\n",
           spectral_rates.j[id_anux], spectral_rates.j_s[id_anux],
           spectral_rates.kappa[id_anux] * 1e7, spectral_rates.kappa_s[id_anux] * 1e7);

    // Compute and output gray emissivities and opacities (Eqs. (19)-(23) in Chiesa+25 PRD)
    gray_rates = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    // The numerical factors restore usual units (see output)
    printf("Gray rates assuming equilibrium\n");
    printf("------------------------------\n");
    printf("     eta0          eta1          kappa0        kappa1        scat1\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_nue] * 1e21, gray_rates.eta[id_nue] * 1e21,
           gray_rates.kappa_0_a[id_nue] * 1e7, gray_rates.kappa_a[id_nue] * 1e7,
           gray_rates.kappa_s[id_nue] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_anue] * 1e21, gray_rates.eta[id_anue] * 1e21,
           gray_rates.kappa_0_a[id_anue] * 1e7, gray_rates.kappa_a[id_anue] * 1e7,
           gray_rates.kappa_s[id_anue] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_nux] * 1e21, gray_rates.eta[id_nux] * 1e21,
           gray_rates.kappa_0_a[id_nux] * 1e7, gray_rates.kappa_a[id_nux] * 1e7,
           gray_rates.kappa_s[id_nux] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n\n",
           gray_rates.eta_0[id_anux] * 1e21, gray_rates.eta[id_anux] * 1e21,
           gray_rates.kappa_0_a[id_anux] * 1e7, gray_rates.kappa_a[id_anux] * 1e7,
           gray_rates.kappa_s[id_anux] * 1e7);

    ////////////////////////////////////////////////////////////////////
    // PART 2: compute rates reconstructing the neutrino distribution //
    //         function from neutrino number and energy densities     //
    ////////////////////////////////////////////////////////////////////

    // Set neutrino number and energy densities
    my_grey_opacity_params.m1_pars.n[id_nue]    = n_nue * 1e-21; // nm^-3
    my_grey_opacity_params.m1_pars.J[id_nue]    = j_nue * 1e-21 * kBS_MeV; // g s^-2 nm^-1
    my_grey_opacity_params.m1_pars.chi[id_nue]  = chi_nue;
    my_grey_opacity_params.m1_pars.n[id_anue]   = n_anue * 1e-21; // nm^-3
    my_grey_opacity_params.m1_pars.J[id_anue]   = j_anue * 1e-21 * kBS_MeV; // g s^-2 nm^-1
    my_grey_opacity_params.m1_pars.chi[id_anue] = chi_anue;
    my_grey_opacity_params.m1_pars.n[id_nux]    = n_nux * 1e-21; // nm^-3
    my_grey_opacity_params.m1_pars.J[id_nux]    = j_nux * 1e-21 * kBS_MeV; // g s^-2 nm^-1
    my_grey_opacity_params.m1_pars.chi[id_nux]  = chi_nux;
    my_grey_opacity_params.m1_pars.n[id_anux]   = n_nux * 1e-21; // nm^-3
    my_grey_opacity_params.m1_pars.J[id_anux]   = j_nux * 1e-21 * kBS_MeV; // g s^-2 nm^-1
    my_grey_opacity_params.m1_pars.chi[id_anux] = chi_nux;

    // Compute neutrino distribution parameters from neutrino number/energy
    // densities
    my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
        &my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

    // Compute and output spectral emissivities and inverse mean free paths (not in the stimulated
    // absorption formalism)
    spectral_rates = ComputeSpectralOpacitiesNotStimulatedAbs(nu_energy, &my_quadrature, &my_grey_opacity_params);

    // The numerical factors restore usual units (see output)
    printf("Spectral rates reconstructing distribution function\n");
    printf("------------------------------\n");
    printf("     j             j_s           kappa         kappa_s\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_nue], spectral_rates.j_s[id_nue],
           spectral_rates.kappa[id_nue] * 1e7, spectral_rates.kappa_s[id_nue] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_anue], spectral_rates.j_s[id_anue],
           spectral_rates.kappa[id_anue] * 1e7, spectral_rates.kappa_s[id_anue] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e\n",
           spectral_rates.j[id_nux], spectral_rates.j_s[id_nux],
           spectral_rates.kappa[id_nux] * 1e7, spectral_rates.kappa_s[id_nux] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e\n\n",
           spectral_rates.j[id_anux], spectral_rates.j_s[id_anux],
           spectral_rates.kappa[id_anux] * 1e7, spectral_rates.kappa_s[id_anux] * 1e7);

    // Compute and output gray emissivities and opacities (Eqs. (19)-(23) in Chiesa+25 PRD)
    gray_rates = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    // The numerical factors restore usual units (see output)
    printf("Gray rates reconstructing distribution function\n");
    printf("----------------------------------------------\n");
    printf("     eta0          eta1          kappa0        kappa1        scat1\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_nue] * 1e21, gray_rates.eta[id_nue] * 1e21,
           gray_rates.kappa_0_a[id_nue] * 1e7, gray_rates.kappa_a[id_nue] * 1e7,
           gray_rates.kappa_s[id_nue] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_anue] * 1e21, gray_rates.eta[id_anue] * 1e21,
           gray_rates.kappa_0_a[id_anue] * 1e7, gray_rates.kappa_a[id_anue] * 1e7,
           gray_rates.kappa_s[id_anue] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           gray_rates.eta_0[id_nux] * 1e21, gray_rates.eta[id_nux] * 1e21,
           gray_rates.kappa_0_a[id_nux] * 1e7, gray_rates.kappa_a[id_nux] * 1e7,
           gray_rates.kappa_s[id_nux] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n\n\n",
           gray_rates.eta_0[id_anux] * 1e21, gray_rates.eta[id_anux] * 1e21,
           gray_rates.kappa_0_a[id_anux] * 1e7, gray_rates.kappa_a[id_anux] * 1e7,
           gray_rates.kappa_s[id_anux] * 1e7);

    printf("Units\n"
           "-----\n"
           "Spectral emissivity 'j'/'j_s'   :           s^-1\n"
           "Spectral imfp 'kappa'/'kappa_s' :          cm^-1\n"
           "Gray number emissivity 'eta0'   :     cm^-3 s^-1\n"
           "Gray energy emissivity 'eta1'   : MeV cm^-3 s^-1\n"
           "Gray number opacity 'kappa0'    :     cm^-1 s^-1\n"
           "Gray energy opacity 'kappa1'    : MeV cm^-1 s^-1\n"
           "Gray scattering opacity 'scat1' : MeV cm^-1 s^-1\n");

    return 0;
}
