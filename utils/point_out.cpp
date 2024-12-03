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
    double cactus2cgsRho = 6.1762691458861632e+17;
    double cactus2cgsEps = 8.987551787368178e+20;

    double cgs2cactusRho    = 1.619100425158886e-18;
    double cgs2cactusPress  = 1.8014921788094724e-39;
    double cgs2cactusEps    = 1.112650056053618e-21;
    double cgs2cactusMass   = 5.0278543128934301e-34;
    double cgs2cactusEnergy = 5.5942423830703013e-55;
    double cgs2cactusTime   = 203012.91587112966e0;
    double cgs2cactusLength = 6.7717819596091924e-06;

    double mev_to_erg = 1.60217733e-6;
    double erg_to_mev = 6.24150636e5;

    // char* name[6] = {"A", "B", "C", "D", "E", "F"};
    double rho[6]    = {6.9e14, 9.81e13, 9.87e12, 1e12, 1e11, 1.01e10};
    double T[6]      = {12.39, 16.63, 8.74, 6.61, 3.60, 2.17};
    double ye[6]     = {0.07, 0.06, 0.06, 0.11, 0.14, 0.19};
    double ynue[6]   = {0.09e-4, 0.36e-3, 0.82e-3, 0.72e-2, 10.92e-3, 12.56e-3};
    double yanue[6]  = {2.93e-4, 2.21e-3, 2.01e-3, 0.45e-2, 7.88e-3, 13e-3};
    double ynux[6]   = {0.54e-4, 0.92e-3, 0.84e-3, 0.15e-2, 0.29e-3, 0.36e-3};
    double znue[6]   = {0.03e-2, 0.18e-1, 0.22e-1, 1.51e-1, 12.48e-2, 11.89e-2};
    double zanue[6]  = {1.29e-2, 1.23e-1, 0.57e-1, 0.93e-1, 10.21e-2, 15.56e-2};
    double znux[6]   = {.21e-2, 0.48e-1, 0.23e-1, 0.39e-1, 0.83e-2, 0.96e-2};
    double chinue[6] = {1. / 3., 1. / 3., 1. / 3., 1. / 3., 0.34, 0.53};
    double chianue[6] = {1. / 3., 1. / 3., 1. / 3., 1. / 3., 0.37, 0.57};
    double chinux[6]  = {1. / 3., 1. / 3., 1. / 3., 1. / 3., 0.36, 0.59};
    // double chinue[6]  = {0.33, 0.33, 0.33, 0.33, 0.34, 0.53};
    // double chianue[6] = {0.33, 0.33, 0.33, 0.33, 0.37, 0.57};
    // double chinux[6]  = {0.33, 0.33, 0.33, 0.33, 0.36, 0.59};
    double deltau[6] = {1.9e1, 3.48e1, 4.98, 0.47, 0.45e-1, 0.72e-3};
    double mn[6]     = {2.81e2, 7.36e2, 9.16e2, 9.37e2, 9.39e2, 9.40e2};
    double mp[6]     = {2.80e2, 7.35e2, 9.15e2, 9.36e2, 9.38e2, 9.38e2};
    double hatmu[6]  = {2.10e2, 9.9e1, 4.11e1, 1.76e1, 8.44, 4.57};
    double mue[6]    = {1.87e2, 8.23e1, 3.65e1, 1.94e1, 9.05, 4.17};

    int i = 0, r;

    double nb;

    MyQuadrature my_quadrature;

    M1Opacities opacities;

    double n_nue, n_anue, n_nux, j_nue, j_anue, j_nux;

    FILE* outfile1;
    FILE* outfile2;

    my_quadrature.nx   = 20;
    my_quadrature.dim  = 1;
    my_quadrature.type = kGauleg;
    my_quadrature.x1   = 0.;
    my_quadrature.x2   = 1.;
    GaussLegendre(&my_quadrature);

    // opacity params structure
    GreyOpacityParams my_grey_opacity_params = {0};

    // reaction flags
    my_grey_opacity_params.opacity_flags = opacity_flags_default_none;
    my_grey_opacity_params.opacity_pars  = opacity_params_default_none;
    my_grey_opacity_params.opacity_pars.use_dm_eff = 0;

    for (int CO = 0; CO < 6; ++CO)
    {

#define STRh(X) #X
#define STR(X) STRh(X)
        char buf1[100];
        char buf2[100];

        switch (CO)
        {
        case 0:
            my_grey_opacity_params.opacity_pars.use_decay = 0;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 0;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 0;
            my_grey_opacity_params.opacity_pars.use_dU    = 0;
            sprintf(buf1, "nm_eq_none_%s.dat", STR(REAL_TYPE));
            sprintf(buf1, "nm_eq_none_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_none_%s.dat", STR(REAL_TYPE));
            break;
        case 1:
            my_grey_opacity_params.opacity_pars.use_decay = 1;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 0;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 0;
            my_grey_opacity_params.opacity_pars.use_dU    = 0;
            sprintf(buf1, "nm_eq_decay_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_decay_%s.dat", STR(REAL_TYPE));
            break;
        case 2:
            my_grey_opacity_params.opacity_pars.use_decay = 0;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 1;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 0;
            my_grey_opacity_params.opacity_pars.use_dU    = 0;
            sprintf(buf1, "nm_eq_WMab_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_WMab_%s.dat", STR(REAL_TYPE));
            break;
        case 3:
            my_grey_opacity_params.opacity_pars.use_decay = 0;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 0;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 1;
            my_grey_opacity_params.opacity_pars.use_dU    = 0;
            sprintf(buf1, "nm_eq_WMsc_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_WMsc_%s.dat", STR(REAL_TYPE));
            break;
        case 4:
            my_grey_opacity_params.opacity_pars.use_decay = 0;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 0;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 0;
            my_grey_opacity_params.opacity_pars.use_dU    = 1;
            sprintf(buf1, "nm_eq_dU_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_dU_%s.dat", STR(REAL_TYPE));
            break;
        case 5:
            my_grey_opacity_params.opacity_pars.use_decay = 1;
            my_grey_opacity_params.opacity_pars.use_WM_ab = 1;
            my_grey_opacity_params.opacity_pars.use_WM_sc = 1;
            my_grey_opacity_params.opacity_pars.use_dU    = 1;
            sprintf(buf1, "nm_eq_all_%s.dat", STR(REAL_TYPE));
            sprintf(buf2, "nm_m1_all_%s.dat", STR(REAL_TYPE));
            break;
        }

        outfile1 = fopen(buf1, "w");
        outfile2 = fopen(buf2, "w");

        for (r = 0; r < 6; ++r)
        {

            my_grey_opacity_params.opacity_flags.use_abs_em          = 0;
            my_grey_opacity_params.opacity_flags.use_brem            = 0;
            my_grey_opacity_params.opacity_flags.use_pair            = 0;
            my_grey_opacity_params.opacity_flags.use_iso             = 0;
            my_grey_opacity_params.opacity_flags.use_inelastic_scatt = 0;

            switch (r)
            {
            case 0:
                my_grey_opacity_params.opacity_flags.use_abs_em = 1;
                break;
            case 1:
                my_grey_opacity_params.opacity_flags.use_brem = 1;
                break;
            case 2:
                my_grey_opacity_params.opacity_flags.use_pair = 1;
                break;
            case 3:
                my_grey_opacity_params.opacity_flags.use_iso = 1;
                break;
            case 4:
                my_grey_opacity_params.opacity_flags.use_inelastic_scatt = 1;
                break;
            case 5:
                my_grey_opacity_params.opacity_flags.use_abs_em          = 1;
                my_grey_opacity_params.opacity_flags.use_brem            = 1;
                my_grey_opacity_params.opacity_flags.use_pair            = 1;
                my_grey_opacity_params.opacity_flags.use_iso             = 1;
                my_grey_opacity_params.opacity_flags.use_inelastic_scatt = 1;
                break;
            }

            for (i = 0; i < 6; ++i)
            {
                nb = rho[i] / 1.66053906892e-24 * 1e-21;

                // populate EOS quantities
                my_grey_opacity_params.eos_pars.mu_e = mue[i];
                my_grey_opacity_params.eos_pars.mu_p = 0;
                my_grey_opacity_params.eos_pars.mu_n = hatmu[i];
                my_grey_opacity_params.eos_pars.temp = T[i];
                my_grey_opacity_params.eos_pars.yp   = ye[i];
                my_grey_opacity_params.eos_pars.yn   = 1. - ye[i];
                my_grey_opacity_params.eos_pars.nb   = nb;
                my_grey_opacity_params.eos_pars.dU   = deltau[i];

                my_grey_opacity_params.distr_pars =
                    NuEquilibriumParams(&my_grey_opacity_params.eos_pars);
                ComputeM1DensitiesEq(&my_grey_opacity_params.eos_pars,
                                     &my_grey_opacity_params.distr_pars,
                                     &my_grey_opacity_params.m1_pars);
                my_grey_opacity_params.m1_pars.chi[id_nue] =
                    0.333333333333333333333333333;
                my_grey_opacity_params.m1_pars.chi[id_anue] =
                    0.333333333333333333333333333;
                my_grey_opacity_params.m1_pars.chi[id_nux] =
                    0.333333333333333333333333333;
                my_grey_opacity_params.m1_pars.chi[id_anux] =
                    0.333333333333333333333333333;
                my_grey_opacity_params.m1_pars.J[id_nue] *= kBS_MeV;
                my_grey_opacity_params.m1_pars.J[id_anue] *= kBS_MeV;
                my_grey_opacity_params.m1_pars.J[id_nux] *= kBS_MeV;
                my_grey_opacity_params.m1_pars.J[id_anux] *= kBS_MeV;

                opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                               &my_grey_opacity_params);

                fprintf(outfile1,
                        "%.17e %.17e %.17e %.17e %.17e "
                        "%.17e %.17e %.17e %.17e %.17e "
                        "%.17e %.17e %.17e %.17e %.17e\n",
                        opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
                        opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
                        opacities.kappa_s[0] * 1e7, opacities.eta_0[1] * 1e21,
                        opacities.eta[1] * 1e21, opacities.kappa_0_a[1] * 1e7,
                        opacities.kappa_a[1] * 1e7, opacities.kappa_s[1] * 1e7,
                        opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
                        opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
                        opacities.kappa_s[2] * 1e7);

                n_nue  = ynue[i] * nb;
                n_anue = yanue[i] * nb;
                n_nux  = ynux[i] * nb;
                j_nue  = znue[i] * nb * kBS_MeV;
                j_anue = zanue[i] * nb * kBS_MeV;
                j_nux  = znux[i] * nb * kBS_MeV;

                my_grey_opacity_params.m1_pars.n[id_nue]    = n_nue;
                my_grey_opacity_params.m1_pars.J[id_nue]    = j_nue;
                my_grey_opacity_params.m1_pars.chi[id_nue]  = chinue[i];
                my_grey_opacity_params.m1_pars.n[id_anue]   = n_anue;
                my_grey_opacity_params.m1_pars.J[id_anue]   = j_anue;
                my_grey_opacity_params.m1_pars.chi[id_anue] = chianue[i];
                my_grey_opacity_params.m1_pars.n[id_nux]    = n_nux;
                my_grey_opacity_params.m1_pars.J[id_nux]    = j_nux;
                my_grey_opacity_params.m1_pars.chi[id_nux]  = chinux[i];
                my_grey_opacity_params.m1_pars.n[id_anux]   = n_nux;
                my_grey_opacity_params.m1_pars.J[id_anux]   = j_nux;
                my_grey_opacity_params.m1_pars.chi[id_anux] = chinux[i];
                my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
                    &my_grey_opacity_params.m1_pars,
                    &my_grey_opacity_params.eos_pars);

                opacities = ComputeM1Opacities(&my_quadrature, &my_quadrature,
                                               &my_grey_opacity_params);

                fprintf(outfile2,
                        "%.17e %.17e %.17e %.17e %.17e "
                        "%.17e %.17e %.17e %.17e %.17e "
                        "%.17e %.17e %.17e %.17e %.17e\n",
                        opacities.eta_0[0] * 1e21, opacities.eta[0] * 1e21,
                        opacities.kappa_0_a[0] * 1e7, opacities.kappa_a[0] * 1e7,
                        opacities.kappa_s[0] * 1e7, opacities.eta_0[1] * 1e21,
                        opacities.eta[1] * 1e21, opacities.kappa_0_a[1] * 1e7,
                        opacities.kappa_a[1] * 1e7, opacities.kappa_s[1] * 1e7,
                        opacities.eta_0[2] * 1e21, opacities.eta[2] * 1e21,
                        opacities.kappa_0_a[2] * 1e7, opacities.kappa_a[2] * 1e7,
                        opacities.kappa_s[2] * 1e7);
            }
        }

        fflush(outfile1);
        fflush(outfile2);
        fclose(outfile1);
        fclose(outfile2);
    }
    return 0;
}
