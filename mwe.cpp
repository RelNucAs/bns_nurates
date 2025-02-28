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
    // Constants
    const double nb = 4.155277120030484e+38;
    const double T = 12.39;
    const double ye = 0.07;
    const double mu_e = 1.87e2;
    const double mu_p = 0.;
    const double mu_n = 2.1e2;
    const double dm = 0.;
    const double dU = 1.9e1;
    const double n_nue  = 3.739749408027436e+33;
    const double n_anue = 1.2174961961689319e+35;
    const double n_nux  = 2.2438496448164613e+34;
    const double j_nue  = 1.246583136009145e+35;
    const double j_anue = 5.360307484839323e+36;
    const double j_nux  = 8.726081952064015e+35;
    const double chi_nue  = 1. / 3.;
    const double chi_anue = 1. / 3.;
    const double chi_nux  = 1. / 3.;

    // Create a quadrature struct and populate it with data relative to
    // Gauss-Legendre quadrature
    MyQuadrature my_quadrature;
    my_quadrature.nx   = 6;
    my_quadrature.dim  = 1;
    my_quadrature.type = kGauleg;
    my_quadrature.x1   = 0.;
    my_quadrature.x2   = 1.;
    GaussLegendre(&my_quadrature);

    // Create an opacities struct that will hold the final result of the
    // computation
    M1Opacities opacities;

    // Create an opacity params structure, to activate/deactivate specific
    // reactions or corrections and pass physical parameters
    GreyOpacityParams my_grey_opacity_params = {0};
    my_grey_opacity_params.opacity_flags     = opacity_flags_default_none;
    my_grey_opacity_params.opacity_pars      = opacity_params_default_none;
    my_grey_opacity_params.opacity_pars.use_dm_eff =
        0; // Do not use effective mass correction to beta-processes
    my_grey_opacity_params.opacity_pars.use_dU =
        1; // Use effective potential correction to beta-processes
    my_grey_opacity_params.opacity_pars.use_decay =
        1; // Activate beta-decay contribution
    my_grey_opacity_params.opacity_pars.use_WM_ab =
        1; // Activate beta-processes weak magnetism correction
    my_grey_opacity_params.opacity_pars.use_WM_sc =
        1; // Activate scattering weak magnetism correction
    my_grey_opacity_params.opacity_flags.use_abs_em =
        1; // Activate beta-processes
    my_grey_opacity_params.opacity_flags.use_brem =
        1; // Activate Bremsstrahlung
    my_grey_opacity_params.opacity_flags.use_pair =
        1; // Activate pair processes
    my_grey_opacity_params.opacity_flags.use_iso =
        1; // Activate scattering on nucleons
    my_grey_opacity_params.opacity_flags.use_inelastic_scatt =
        1; // Activate scattering on leptons
    my_grey_opacity_params.eos_pars.mu_e =
        mu_e; // Set electron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.mu_p =
        mu_p; // Set proton chemical potential (in MeV) (NOTE: reactions only
            // depend on the difference between mu_n and mu_p, hence mu_p can be
            // set to 0)
    my_grey_opacity_params.eos_pars.mu_n =
        mu_n; // Set neutron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.temp = T; // Set temperature (in MeV)
    my_grey_opacity_params.eos_pars.yp   = ye;  // Set proton fraction
    my_grey_opacity_params.eos_pars.yn   = 1. - ye; // Set neutron fraction
    my_grey_opacity_params.eos_pars.nb =
        nb * 1e-21; // Set baryon number density (in baryon/nm^3)
    my_grey_opacity_params.eos_pars.dU =
        dU; // Set effective potential difference (in MeV)
    my_grey_opacity_params.eos_pars.dm_eff =
        dm; // Set effective mass difference (in MeV)

    printf("\nInput parameters\n");
    printf("----------------\n");
    printf("Baryon number density                           : %13.6e (cm^-3)\n", nb);
    printf("Temperature                                     : %13.6e (MeV)\n", T);
    printf("Electron fraction                               : %13.3f\n", ye);
    printf("Relativistic electron chemical potential        : %13.6e (MeV)\n", mu_e);
    printf("Relativistic proton chemical potential          : %13.6e (MeV)\n", mu_p);
    printf("Relativistic neutron chemical potential         : %13.6e (MeV)\n", mu_n);
    printf("Effective nucleon mass difference (n - p)       : %13.6e (MeV)\n", dm);
    printf("Effective nucleon potential difference (n - p)  : %13.6e (MeV)\n", dU);
    printf("Electron neutrinos number density 'n'           : %13.6e (cm^-3)\n", n_nue);
    printf("Electron antineutrinos number density 'n'       : %13.6e (cm^-3)\n", n_anue);
    printf("Heavy neutrinos number density 'n'              : %13.6e (cm^-3)\n", n_nux);
    printf("Electron neutrinos energy density 'J'           : %13.6e (MeV cm^-3)\n", j_nue);
    printf("Electron antineutrinos energy density 'J'       : %13.6e (MeV cm^-3)\n", j_anue);
    printf("Heavy neutrinos energy density 'J'              : %13.6e (MeV cm^-3)\n", j_nux);
    printf("Electron neutrinos Eddington parameter 'chi'    : %13.11f\n", chi_nue);
    printf("Electron antineutrinos Eddington parameter 'chi': %13.11f\n", chi_anue);
    printf("Heavy neutrinos Eddington parameter 'chi'       : %13.11f\n\n", chi_nux);

    ////////////////////////////////////////////////////////////////
    // PART 1: compute opacities assuming equilibrium with matter //
    ////////////////////////////////////////////////////////////////

    // Compute neutrino distribution parameters assuming equilibrium
    my_grey_opacity_params.distr_pars =
        NuEquilibriumParams(&my_grey_opacity_params.eos_pars);

    // Compute neutrino number and energy densities assuming equilibrium
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
    printf("n  %13.6e  %13.6e   %13.6e  %13.6e\n", n_nue_eq, n_anue_eq, n_nux_eq, n_anux_eq);
    printf("J  %13.6e  %13.6e   %13.6e  %13.6e\n", J_nue_eq, J_anue_eq, J_nux_eq, J_anux_eq);
    printf("chi %12.10f   %12.10f    %12.10f   %12.10f\n\n", 1./3., 1./3., 1./3., 1./3.);

    // Compute and output opacities
    opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    // The numerical factors restore usual units (see output)
    printf("Opacities assuming equilibrium\n");
    printf("------------------------------\n");
    printf("     eta0          eta1          kappa0        kappa1        scat1\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
           opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
           opacities.kappa_s[0] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[1] * 1e21, opacities.eta[1] * 1e21,
           opacities.kappa_0_a[1] * 1e7, opacities.kappa_a[1] * 1e7,
           opacities.kappa_s[1] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
           opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
           opacities.kappa_s[2] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n\n",
           opacities.eta_0[3] * 1e21, opacities.eta[3] * 1e21,
           opacities.kappa_0_a[3] * 1e7, opacities.kappa_a[3] * 1e7,
           opacities.kappa_s[3] * 1e7);

    ////////////////////////////////////////////////////////////////////////
    // PART 2: compute opacities reconstructing the neutrino distribution //
    //         function from neutrino number and energy densities         //
    ////////////////////////////////////////////////////////////////////////

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

    // Compute and output opacities
    opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    // The numerical factors restore usual units (see output)
    printf("Opacities reconstructing distribution function\n");
    printf("----------------------------------------------\n");
    printf("     eta0          eta1          kappa0        kappa1        scat1\n");
    printf(" nue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
           opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
           opacities.kappa_s[0] * 1e7);
    printf("anue %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[1] * 1e21, opacities.eta[1] * 1e21,
           opacities.kappa_0_a[1] * 1e7, opacities.kappa_a[1] * 1e7,
           opacities.kappa_s[1] * 1e7);
    printf(" nux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n",
           opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
           opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
           opacities.kappa_s[2] * 1e7);
    printf("anux %-13.6e %-13.6e %-13.6e %-13.6e %-13.6e\n\n\n",
           opacities.eta_0[3] * 1e21, opacities.eta[3] * 1e21,
           opacities.kappa_0_a[3] * 1e7, opacities.kappa_a[3] * 1e7,
           opacities.kappa_s[3] * 1e7);

    printf("Units\n"
           "-----\n"
           "Number emissivity 'eta0'  :     cm^-3 s^-1\n"
           "Energy emissivity 'eta1'  : MeV cm^-3 s^-1\n"
           "Number opacity 'kappa0'   :     cm^-1 s^-1\n"
           "Energy opacity 'kappa1'   : MeV cm^-1 s^-1\n"
           "Scattering opacity 'scat1': MeV cm^-1 s^-1\n");

    return 0;
}
