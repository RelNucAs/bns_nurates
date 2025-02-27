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

    // Create a quadrature struct and populate it with data relative to
    // Gauss-Legendre quadrature
    MyQuadrature my_quadrature;
    my_quadrature.nx   = 20;
    my_quadrature.dim  = 1;
    my_quadrature.type = kGauleg;
    my_quadrature.x1   = 0.;
    my_quadrature.x2   = 1.;
    GaussLegendre(&my_quadrature);

    // Create an opacities struct that will hold the final result of the
    // computation
    M1Opacities opacities;
    double cactus2cgsRho = 6.1762691458861632e+17;
    double cactus2cgsEps = 8.987551787368178e+20;

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
        1.87e2; // Set electron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.mu_p =
        0.; // Set proton chemical potential (in MeV) (NOTE: reactions only
            // depend on the difference between mu_n and mu_p, hence mu_p can be
            // set to 0)
    my_grey_opacity_params.eos_pars.mu_n =
        2.10e2; // Set neutron chemical potential (in MeV)
    my_grey_opacity_params.eos_pars.temp = 12.39; // Set temperature (in MeV)
    my_grey_opacity_params.eos_pars.yp   = 0.07;  // Set proton fraction
    my_grey_opacity_params.eos_pars.yn   = 1. - 0.07; // Set neutron fraction
    my_grey_opacity_params.eos_pars.nb =
        4.155277120030484e+17; // Set baryon number density (in baryon/nm^3)
    my_grey_opacity_params.eos_pars.dU =
        19.; // Set effective potential difference (in MeV)

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
    // Restore correct units
    my_grey_opacity_params.m1_pars.J[id_nue] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_anue] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_nux] *= kBS_MeV;
    my_grey_opacity_params.m1_pars.J[id_anux] *= kBS_MeV;

    // Compute and output opacities
    opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    printf("Opacities assuming equilibrium:\n");
    printf("eta0_nue: %.17e, eta1_nue: %.17e, kappa0_nue: %.17e, kappa1_nue: "
           "%.17e, scat1_nue: %.17e\n",
           opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
           opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
           opacities.kappa_s[0] * 1e7);
    printf("eta0_anue: %.17e, eta1_anue: %.17e, kappa0_anue: %.17e, "
           "kappa1_anue: %.17e, scat1_anue: %.17e\n",
           opacities.eta_0[1] * 1e21, opacities.eta[1] * 1e21,
           opacities.kappa_0_a[1] * 1e7, opacities.kappa_a[1] * 1e7,
           opacities.kappa_s[1] * 1e7);
    printf("eta0_nux: %.17e, eta1_nux: %.17e, kappa0_nux: %.17e, kappa1_nux: "
           "%.17e, scat1_nux: %.17e\n",
           opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
           opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
           opacities.kappa_s[2] * 1e7);
    printf("eta0_anux: %.17e, eta1_anux: %.17e, kappa0_anux: %.17e, "
           "kappa1_anux: %.17e, scat1_anux: %.17e\n",
           opacities.eta_0[3] * 1e21, opacities.eta[3] * 1e21,
           opacities.kappa_0_a[3] * 1e7, opacities.kappa_a[3] * 1e7,
           opacities.kappa_s[3] * 1e7);

    ////////////////////////////////////////////////////////////////////////
    // PART 2: compute opacities reconstructing the neutrino distribution //
    //         function from neutrino number and energy densities         //
    ////////////////////////////////////////////////////////////////////////

    double n_nue  = 3739749408027.4355;
    double n_anue = 121749619616893.19;
    double n_nux  = 22438496448164.613;
    double j_nue  = 124658313600914.5 * kBS_MeV;
    double j_anue = 5360307484839324.0 * kBS_MeV;
    double j_nux  = 872608195206401.5 * kBS_MeV;

    // Set neutrino number and energy densities
    my_grey_opacity_params.m1_pars.n[id_nue]    = n_nue;
    my_grey_opacity_params.m1_pars.J[id_nue]    = j_nue;
    my_grey_opacity_params.m1_pars.chi[id_nue]  = 0.34;
    my_grey_opacity_params.m1_pars.n[id_anue]   = n_anue;
    my_grey_opacity_params.m1_pars.J[id_anue]   = j_anue;
    my_grey_opacity_params.m1_pars.chi[id_anue] = 0.34;
    my_grey_opacity_params.m1_pars.n[id_nux]    = n_nux;
    my_grey_opacity_params.m1_pars.J[id_nux]    = j_nux;
    my_grey_opacity_params.m1_pars.chi[id_nux]  = 0.34;
    my_grey_opacity_params.m1_pars.n[id_anux]   = n_nux;
    my_grey_opacity_params.m1_pars.J[id_anux]   = j_nux;
    my_grey_opacity_params.m1_pars.chi[id_anux] = 0.34;

    // Compute neutrino distribution parameters from neutrino number/energy
    // densities
    my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
        &my_grey_opacity_params.m1_pars, &my_grey_opacity_params.eos_pars);

    // Compute and output opacities
    opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                   &my_grey_opacity_params);

    printf("Opacities reconstructing distribution function:\n");
    printf("eta0_nue: %.17e, eta1_nue: %.17e, kappa0_nue: %.17e, kappa1_nue: "
           "%.17e, scat1_nue: %.17e\n",
           opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
           opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
           opacities.kappa_s[0] * 1e7);
    printf("eta0_anue: %.17e, eta1_anue: %.17e, kappa0_anue: %.17e, "
           "kappa1_anue: %.17e, scat1_anue: %.17e\n",
           opacities.eta_0[1] * 1e21, opacities.eta[1] * 1e21,
           opacities.kappa_0_a[1] * 1e7, opacities.kappa_a[1] * 1e7,
           opacities.kappa_s[1] * 1e7);
    printf("eta0_nux: %.17e, eta1_nux: %.17e, kappa0_nux: %.17e, kappa1_nux: "
           "%.17e, scat1_nux: %.17e\n",
           opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
           opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
           opacities.kappa_s[2] * 1e7);
    printf("eta0_anux: %.17e, eta1_anux: %.17e, kappa0_anux: %.17e, "
           "kappa1_anux: %.17e, scat1_anux: %.17e\n",
           opacities.eta_0[3] * 1e21, opacities.eta[3] * 1e21,
           opacities.kappa_0_a[3] * 1e7, opacities.kappa_a[3] * 1e7,
           opacities.kappa_s[3] * 1e7);

    return 0;
}
