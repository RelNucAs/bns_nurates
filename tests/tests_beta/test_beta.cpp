//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_beta.c
//  \brief test routines for beta processes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Kokkos_Core.hpp>

#include "constants.hpp"
#include "weak_magnetism.hpp"
#include "formulas_beta_iso.hpp"

void generate_comparison_data(const bool use_dU)
{

    using DevExeSpace  = Kokkos::DefaultExecutionSpace;
    using DevMemSpace  = Kokkos::DefaultExecutionSpace::memory_space;
    using HostMemSpace = Kokkos::HostSpace;
    using LayoutWrapper =
        Kokkos::LayoutRight; // increments last index fastest views defined like

    // Define path of input file
    char filepath[300] = {'\0'};
    char filedir[300]  = SOURCE_DIR;
    char outname[200]  = {'\0'};

    if (use_dU == false)
    {
        strcpy(outname, "/inputs/nurates_CCSN/nurates_1.008E+01.txt");
    }
    else
    {
        strcpy(outname, "/inputs/nurates_CCSN/nurates_1.008E+01_dU.txt");
    }

    strcat(filepath, filedir);
    strcat(filepath, outname);

    printf("Input data file: %s\n", filepath);
    printf("\n");

    // Number of entries (must be compatible with # of rows in the file)
    const int n_data = 102;

    // Neutrino energy (the comparison is done at fixed energy)
    double omega; // [MeV]

    // Weak magnetism correction (code output)
    double w, wbar;

    // Data arrays (rates are not corrected for the weak magnetism, WM
    // corrections are provided separately)
    Kokkos::View<int*, LayoutWrapper, HostMemSpace> h_zone(
        "h_zone", n_data); // zone (counter)
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_r(
        "h_r", n_data); // radius [cm]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_rho(
        "h_rho", n_data); // mass density [g/cm3]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_temp(
        "h_temp", n_data); // temperature [MeV]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_ye(
        "h_ye", n_data); // electron fraction [#/baryon]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_e(
        "h_mu_e", n_data); // electron chemical potential (with rest mass) [MeV]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_mu_hat(
        "h_mu_hat",
        n_data); // neutron - proton chemical potential (w/o rest mass) [MeV]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_yh(
        "h_yh", n_data); // heavy abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_ya(
        "h_ya", n_data); // alpha abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_yp(
        "h_yp", n_data); // proton abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_yn(
        "h_yn", n_data); // neuton abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_du(
        "h_du", n_data); // correction to nucleon chemical potential due to
                         // interactions [MeV]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_em_nue(
        "h_em_nue", n_data); // electron neutrino emissivity [1/s]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_ab_nue(
        "h_ab_nue", n_data); // electron neutrino inverse mean free path [1/cm]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_em_anue(
        "h_em_anue", n_data); // electron antineutrino emissivity [1/s]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_ab_anue(
        "h_ab_anue",
        n_data); // electron antineutrino inverse mean free path [1/cm]
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_w_nue("h_w_nue",
                                                               n_data),
        h_w_anue(
            "h_w_anue",
            n_data); // weak magnetism correction to beta rates (electron-type
                     // neutrino and antineutrino, respectvely)
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_r_is(
        "h_r_is", n_data); // source term for quasi-elastic scattering (Eq.(A41)
                           // in Bruenn1985?)
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_w0_n("h_w0_n", n_data),
        h_w0_p("h_w0_p", n_data); // weak magnetism correction to zeroth
                                  // Legendre term of scattering kernel (on
                                  // neutrons and protons, respectively)
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_w1_n("h_w1_n", n_data),
        h_w1_p("h_w1_p", n_data); // weak magnetism correction to first Legendre
                                  // term of scattering kernel (on neutrons and
                                  // protons, respectively)

    // Open input file
    FILE* fptr_in;
    fptr_in = fopen(filepath, "r");
    if (fptr_in == NULL)
    {
        printf("%s: The file %s does not exist!\n", __FILE__, filepath);
        exit(1);
    }

    int i = 0;
    char line[1500];

    // Read in the neutrino energy
    auto fvar = fgets(line, sizeof(line), fptr_in);
    sscanf(line + 14, "%lf\n", &omega);

    printf("Neutrino energy: %.5e MeV\n", omega);
    printf("\n");

    // Define path of output file
    filepath[0] = '\0';
    if (use_dU == false)
    {
        sprintf(outname, "/tests/tests_beta/test_beta_%.3e.txt", omega);
    }
    else
    {
        sprintf(outname, "/tests/tests_beta/test_beta_%.3e_dU.txt", omega);
    }

    strcat(filepath, filedir);
    strcat(filepath, outname);

    // Open output file
    FILE* fptr_out;
    fptr_out = fopen(filepath, "w");
    if (fptr_out == NULL)
    {
        printf("%s: Error in opening in write mode file %s!\n", __FILE__,
               filepath);
        exit(1);
    }

    printf("Output data file: %s\n", filepath);
    printf("\n");

    // Print neutrino energy in output file
    fprintf(fptr_out, "# E_nu = %.5e MeV\n", omega);
    fprintf(fptr_out, "\n");

    // Print legend as output file header
    fprintf(fptr_out, "#  1: zone\n");
    fprintf(fptr_out, "#  2: r [cm]\n");
    fprintf(fptr_out, "#  3: rho [g/cm^3]\n");
    fprintf(fptr_out, "#  4: T [MeV]\n");
    fprintf(fptr_out, "#  5: Ye\n");
    fprintf(fptr_out, "#  6: mu_e [MeV]\n");
    fprintf(fptr_out, "#  7: mu_hat [MeV]\n");
    fprintf(fptr_out, "#  8: Yh\n");
    fprintf(fptr_out, "#  9: Ya\n");
    fprintf(fptr_out, "# 10: Yp\n");
    fprintf(fptr_out, "# 11: Yn\n");
    fprintf(fptr_out, "# 12: dU [MeV]\n");
    fprintf(fptr_out, "# 13-16: Fortran rates with WM (em_nue [s-1], ab_nue "
                      "[cm-1], em_anue [s-1], ab_anue [cm-1])\n");
    fprintf(fptr_out, "# 17-20: C rates with WM (em_nue [s-1], ab_nue [cm-1], "
                      "em_anue [s-1], ab_anue [cm-1])\n");
    fprintf(fptr_out, "\n");

    // Read in data from profile, compute rates and store results in output file
    while (fgets(line, sizeof(line), fptr_in) != NULL)
    {
        if (line[1] == '#')
        {
            continue;
        }

        sscanf(line,
               "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
               "%lf %lf %lf %lf %lf %lf %lf\n",
               &h_zone(i), &h_r(i), &h_rho(i), &h_temp(i), &h_ye(i), &h_mu_e(i),
               &h_mu_hat(i), &h_yh(i), &h_ya(i), &h_yp(i), &h_yn(i), &h_du(i),
               &h_em_nue(i), &h_ab_nue(i), &h_em_anue(i), &h_ab_anue(i),
               &h_w_nue(i), &h_w_anue(i), &h_r_is(i), &h_w0_n(i), &h_w0_p(i),
               &h_w1_n(i), &h_w1_p(i));

        if (i == 0)
        {
            WMAbsEm(omega, &w, &wbar);

            // Print weak magnetism correction in output file
            fprintf(fptr_out, "# Weak magnetism correction (Fortran):\n");
            fprintf(fptr_out, "# W_nue  = %lf\n", h_w_nue(i));
            fprintf(fptr_out, "# W_anue = %lf\n", h_w_anue(i));
            fprintf(fptr_out, "# Weak magnetism correction (C):\n");
            fprintf(fptr_out, "# W_nue  = %lf \n", w);
            fprintf(fptr_out, "# W_anue = %lf\n", wbar);
            fprintf(fptr_out, "\n");
        }
        i++;
    }

    // Data arrays (rates are not corrected for the weak magnetism, WM
    // corrections are provided separately)
    Kokkos::View<int*, LayoutWrapper, DevMemSpace> zone(
        "zone", n_data); // zone (counter)
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> r("r",
                                                        n_data); // radius [cm]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> rho(
        "rho", n_data); // mass density [g/cm3]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> temp(
        "temp", n_data); // temperature [MeV]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> ye(
        "ye", n_data); // electron fraction [#/baryon]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_e(
        "mu_e", n_data); // electron chemical potential (with rest mass) [MeV]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> mu_hat(
        "mu_hat",
        n_data); // neutron - proton chemical potential (w/o rest mass) [MeV]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> yh(
        "yh", n_data); // heavy abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> ya(
        "ya", n_data); // alpha abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> yp(
        "yp", n_data); // proton abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> yn(
        "yn", n_data); // neuton abundance [#/baryon]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> du(
        "du", n_data); // correction to nucleon chemical potential due to
                       // interactions [MeV]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> em_nue(
        "em_nue", n_data); // electron neutrino emissivity [1/s]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> ab_nue(
        "ab_nue", n_data); // electron neutrino inverse mean free path [1/cm]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> em_anue(
        "em_anue", n_data); // electron antineutrino emissivity [1/s]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> ab_anue(
        "ab_anue",
        n_data); // electron antineutrino inverse mean free path [1/cm]
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> w_nue("w_nue", n_data),
        w_anue(
            "w_anue",
            n_data); // weak magnetism correction to beta rates (electron-type
                     // neutrino and antineutrino, respectvely)
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> r_is(
        "r_is", n_data); // source term for quasi-elastic scattering (Eq.(A41)
                         // in Bruenn1985?)
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> w0_n("w0_n", n_data),
        w0_p("w0_p", n_data); // weak magnetism correction to zeroth Legendre
                              // term of scattering kernel (on neutrons and
                              // protons, respectively)
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> w1_n("w1_n", n_data),
        w1_p("w1_p", n_data); // weak magnetism correction to first Legendre
                              // term of scattering kernel (on neutrons and
                              // protons, respectively)

    Kokkos::deep_copy(zone, h_zone);
    Kokkos::deep_copy(r, h_r);
    Kokkos::deep_copy(rho, h_rho);
    Kokkos::deep_copy(temp, h_temp);
    Kokkos::deep_copy(ye, h_ye);
    Kokkos::deep_copy(mu_e, h_mu_e);
    Kokkos::deep_copy(mu_hat, h_mu_hat);
    Kokkos::deep_copy(yh, h_yh);
    Kokkos::deep_copy(ya, h_ya);
    Kokkos::deep_copy(yp, h_yp);
    Kokkos::deep_copy(yn, h_yn);
    Kokkos::deep_copy(du, h_du);
    Kokkos::deep_copy(em_nue, h_em_nue);
    Kokkos::deep_copy(ab_nue, h_ab_nue);
    Kokkos::deep_copy(em_anue, h_em_anue);
    Kokkos::deep_copy(ab_anue, h_ab_anue);
    Kokkos::deep_copy(w_nue, h_w_nue);
    Kokkos::deep_copy(w_anue, h_w_anue);
    Kokkos::deep_copy(r_is, h_r_is);
    Kokkos::deep_copy(w0_n, h_w0_n);
    Kokkos::deep_copy(w0_p, h_w0_p);
    Kokkos::deep_copy(w1_n, h_w1_n);
    Kokkos::deep_copy(w1_p, h_w1_p);

    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_out_em_id_nue(
        "h_out_em_id_nue", n_data);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_out_abs_id_nue(
        "h_out_abs_id_nue", n_data);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_out_em_id_anue(
        "h_out_em_id_anue", n_data);
    Kokkos::View<double*, LayoutWrapper, HostMemSpace> h_out_abs_id_anue(
        "h_out_abs_id_anue", n_data);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> out_em_id_nue(
        "out_em_id_nue", n_data);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> out_abs_id_nue(
        "out_abs_id_nue", n_data);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> out_em_id_anue(
        "out_em_id_anue", n_data);
    Kokkos::View<double*, LayoutWrapper, DevMemSpace> out_abs_id_anue(
        "out_abs_id_anue", n_data);

    Kokkos::parallel_for(
        "loop_over_opacities_test_beta",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, n_data - 1),
        KOKKOS_LAMBDA(const int& i) {
            // Opacity parameters (correction to rates)
            OpacityParams opacity_pars = {.use_dU     = use_dU,
                                          .use_dm_eff = false,
                                          .use_WM_ab  = true,
                                          .use_decay =
                                              false}; // WM correction turned on

            // EOS parameters
            MyEOSParams eos_pars;

            // Populate EOS params structure with data read from the profile
            eos_pars.nb =
                rho(i) / kBS_Mu * 1e-21; // baryon number density [nm-3]
            eos_pars.temp = temp(i);
            eos_pars.yp   = yp(i);
            eos_pars.yn   = yn(i);
            eos_pars.mu_p = 0.;
            eos_pars.mu_n =
                mu_hat(i) +
                kBS_Q; // N.B.: corrected for the nucleon mass difference
            eos_pars.mu_e = mu_e(i);
            eos_pars.dU   = du(i);

            // Compute emissivity and inverse mean free path for electron-type
            // (anti)neutrinos
            MyOpacity out = AbsOpacity(omega, &opacity_pars, &eos_pars);

            out_em_id_nue(i)   = out.em[id_nue];
            out_abs_id_nue(i)  = out.abs[id_nue];
            out_em_id_anue(i)  = out.em[id_anue];
            out_abs_id_anue(i) = out.abs[id_anue];
        });

    Kokkos::deep_copy(h_out_em_id_nue, out_em_id_nue);
    Kokkos::deep_copy(h_out_abs_id_nue, out_abs_id_nue);
    Kokkos::deep_copy(h_out_em_id_anue, out_em_id_anue);
    Kokkos::deep_copy(h_out_abs_id_anue, out_abs_id_anue);

    for (int i = 0; i < n_data; i++)
    {
        fprintf(fptr_out,
                "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
                "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                h_zone(i), h_r(i), h_rho(i), h_temp(i), h_ye(i), h_mu_e(i),
                h_mu_hat(i), h_yh(i), h_ya(i), h_yp(i), h_yn(i), h_du(i),
                h_em_nue(i), h_ab_nue(i), h_em_anue(i), h_ab_anue(i),
                h_out_em_id_nue(i), h_out_abs_id_nue(i) / kBS_Clight * 1e7,
                h_out_em_id_anue(i), h_out_abs_id_anue(i) / kBS_Clight * 1e7);
    }

    // Close files
    fclose(fptr_in);
    fclose(fptr_out);

    // Check if file was entirely read
    if (n_data != h_zone[i - 1])
    {
        printf("Warning: n_data (%d) different for number of lines in the "
               "input file (%d)",
               n_data, h_zone[i - 1]);
    }

    return;
}


int main()
{

    Kokkos::initialize();

    printf("Spectral rates are compared with the results of a tested Fortran "
           "code\n");
    printf("\n");

    bool use_dU = false;

    printf("============================================================== \n");
    printf("Testing opacities for beta processes without dU correction ... \n");
    printf("============================================================== \n");
    printf("\n");

    generate_comparison_data(use_dU); // dU correction is off

    use_dU = true;

    printf("=========================================================== \n");
    printf("Testing opacities for beta processes with dU correction ... \n");
    printf("=========================================================== \n");
    printf("\n");

    generate_comparison_data(use_dU); // dU correction is on

    Kokkos::finalize();

    return 0;
}
