//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  test_opacities_slice_timing.cpp
//  \brief cost of each reaction in the computation of the M1 emissivities and
//         opacities (m1_opacities.hpp) on the points of a 2D slice of a BNS
//         simulation
//
// The M1 opacities are computed on all the points of the slice with:
//   - no reaction (common overhead: distribution function, M1 densities, ...)
//   - each reaction alone (beta processes, isoenergetic scattering, pair
//     processes, NN bremsstrahlung (HR98 + Fischer16 correction), inelastic
//     scattering)
//   - all the reactions together
// The reactions are selected at run time through OpacityFlags, so that the
// same compiled kernel is used in all cases, as in the production code.
// The timings and the fractions of the total time are printed and written to
// a text file.
//
// The input file is produced by inputs/BNS/2d/slice.py,
// the output file is read by tests/tests_cost/plot_opacities_cost.py
//
// Usage: test_opacities_cost <input_file> <output_file> [n_quad] [n_rep]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <Kokkos_Core.hpp>

#include "bns_nurates.hpp"
#include "m1_opacities.hpp"
#include "integration.hpp"
#include "distribution.hpp"
#include "constants.hpp"
#include "functions.hpp"

using DevExeSpace   = Kokkos::DefaultExecutionSpace;
using DevMemSpace   = Kokkos::DefaultExecutionSpace::memory_space;
using HostMemSpace  = Kokkos::HostSpace;
using LayoutWrapper = Kokkos::LayoutRight;

#define n_eos_vars 7 // nb, T, Yn, Yp, mu_e, mu_n, mu_p
#define n_configs 7  // none, 5 single reactions, all

/* Compute the M1 opacities with the given reactions on all the points */
void ComputeOpacitiesSlice(
    const int npts, const MyQuadrature& quad, const OpacityFlags opacity_flags,
    const OpacityParams opacity_pars,
    const Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace>& eos,
    const Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace>& out)
{
    Kokkos::parallel_for(
        "test_opacities_slice_timing",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, npts),
        KOKKOS_LAMBDA(const int& i) {
            GreyOpacityParams grey_pars = {0};

            grey_pars.opacity_flags = opacity_flags;
            grey_pars.opacity_pars  = opacity_pars;

            // Populate EOS parameters
            grey_pars.eos_pars.nb     = eos(i, 0);
            grey_pars.eos_pars.temp   = eos(i, 1);
            grey_pars.eos_pars.yn     = eos(i, 2);
            grey_pars.eos_pars.yp     = eos(i, 3);
            grey_pars.eos_pars.mu_e   = eos(i, 4);
            grey_pars.eos_pars.mu_n   = eos(i, 5);
            grey_pars.eos_pars.mu_p   = eos(i, 6);
            grey_pars.eos_pars.dU     = 0.;
            grey_pars.eos_pars.dm_eff = 1.29333251; // as in AthenaK

            // Neutrino distribution function at equilibrium
            grey_pars.distr_pars = NuEquilibriumParams(&grey_pars.eos_pars);
            ComputeM1DensitiesEq(&grey_pars.eos_pars, &grey_pars.distr_pars,
                                 &grey_pars.m1_pars);
            for (int idx = 0; idx < total_num_species; idx++)
            {
                grey_pars.m1_pars.chi[idx] = 1. / 3.;
            }

            MyQuadrature my_quad = quad;

            M1Opacities rates =
                ComputeM1Opacities(&my_quad, &my_quad, &grey_pars);

            // Store the results, so that the computation cannot be optimized
            // away
            for (int idx = 0; idx < total_num_species; idx++)
            {
                out(i, idx) = rates.eta_0[idx] + rates.eta[idx] +
                              rates.kappa_0_a[idx] + rates.kappa_a[idx] +
                              rates.kappa_s[idx];
            }
        });
}

int main(int argc, char* argv[])
{
    if (argc < 3 || argc > 5)
    {
        printf("Usage: %s <input_file> <output_file> [n_quad] [n_rep]\n",
               argv[0]);
        return 1;
    }

    const char* input_file  = argv[1];
    const char* output_file = argv[2];
    const int n_quad        = (argc > 3) ? atoi(argv[3]) : 6;
    const int n_rep         = (argc > 4) ? atoi(argv[4]) : 10;

    if (n_quad <= 0 || 2 * n_quad > BS_N_MAX || n_rep <= 0)
    {
        printf("Invalid n_quad (= %d, 2*n_quad <= %d) or n_rep (= %d)\n",
               n_quad, BS_N_MAX, n_rep);
        return 1;
    }

    // Read the EOS quantities on the points of the slice
    FILE* fptr = fopen(input_file, "r");
    if (fptr == NULL)
    {
        printf("%s: The file %s does not exist!\n", __FILE__, input_file);
        return 1;
    }

    std::vector<double> eos_data;
    char line[1000];
    while (fgets(line, sizeof(line), fptr) != NULL)
    {
        if (line[0] == '#')
        {
            continue;
        }
        double v[n_eos_vars];
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2],
                   &v[3], &v[4], &v[5], &v[6]) != n_eos_vars)
        {
            printf("Error while reading line: %s\n", line);
            return 1;
        }
        eos_data.insert(eos_data.end(), v, v + n_eos_vars);
    }
    fclose(fptr);

    const int npts = eos_data.size() / n_eos_vars;

    Kokkos::initialize(argc, argv);
    {
        printf("=================================================== \n");
        printf("Cost of each reaction in m1_opacities on a slice\n");
        printf("=================================================== \n");
        printf("Execution space         = %s\n", DevExeSpace::name());
        printf("Input file              = %s\n", input_file);
        printf("Number of points        = %d\n", npts);
        printf("Quadrature points       = %d\n", n_quad);
        printf("Timed repetitions       = %d\n", n_rep);
        printf("\n");

        Kokkos::View<BS_REAL**, LayoutWrapper, HostMemSpace> h_eos(
            "h_eos", npts, n_eos_vars);
        for (int i = 0; i < npts; i++)
        {
            for (int j = 0; j < n_eos_vars; j++)
            {
                h_eos(i, j) = eos_data[i * n_eos_vars + j];
            }
        }

        Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_eos(
            "d_eos", npts, n_eos_vars);
        Kokkos::deep_copy(d_eos, h_eos);

        Kokkos::View<BS_REAL**, LayoutWrapper, DevMemSpace> d_out(
            "d_out", npts, total_num_species);

        // Quadrature in [0,1], as in AthenaK
        MyQuadrature quad = {.type   = kGauleg,
                             .alpha  = -42.,
                             .dim    = 1,
                             .nx     = n_quad,
                             .ny     = 1,
                             .nz     = 1,
                             .x1     = 0.,
                             .x2     = 1.,
                             .y1     = -42.,
                             .y2     = -42.,
                             .z1     = -42.,
                             .z2     = -42.,
                             .points = {0},
                             .w      = {0}};
        GaussLegendre(&quad);

        // Opacity parameters as in the AthenaK run
        OpacityParams opacity_pars       = opacity_params_default_none;
        opacity_pars.use_decay           = true;
        opacity_pars.brem_implementation = BREM_HR98;
        opacity_pars.use_NN_medium_corr  = true;

        // Reactions in each configuration
        const char* names[n_configs] = {"none",  "beta", "iso",      "pair",
                                        "brem",  "inelastic", "all"};
        OpacityFlags flags[n_configs];
        for (int c = 0; c < n_configs; c++)
        {
            flags[c] = opacity_flags_default_none;
        }
        flags[1].use_abs_em          = 1;
        flags[2].use_iso             = 1;
        flags[3].use_pair            = 1;
        flags[4].use_brem            = 1;
        flags[5].use_inelastic_scatt = 1;
        flags[6]                     = opacity_flags_default_all;

        // Warm-up (not timed)
        for (int c = 0; c < n_configs; c++)
        {
            ComputeOpacitiesSlice(npts, quad, flags[c], opacity_pars, d_eos,
                                  d_out);
        }
        Kokkos::fence();

        // Timed runs, cycling over the configurations
        double t_sum[n_configs], t_sum2[n_configs], t_min[n_configs];
        for (int c = 0; c < n_configs; c++)
        {
            t_sum[c]  = 0.;
            t_sum2[c] = 0.;
            t_min[c]  = 1e300;
        }

        Kokkos::Timer timer;
        for (int r = 0; r < n_rep; r++)
        {
            for (int c = 0; c < n_configs; c++)
            {
                timer.reset();
                ComputeOpacitiesSlice(npts, quad, flags[c], opacity_pars, d_eos,
                                      d_out);
                Kokkos::fence();
                const double t = timer.seconds();

                t_sum[c] += t;
                t_sum2[c] += t * t;
                t_min[c] = (t < t_min[c]) ? t : t_min[c];
            }
        }

        double t_mean[n_configs], t_std[n_configs];
        for (int c = 0; c < n_configs; c++)
        {
            t_mean[c] = t_sum[c] / n_rep;
            t_std[c] =
                sqrt(fmax(t_sum2[c] / n_rep - t_mean[c] * t_mean[c], 0.));
        }

        // fraction: time of the configuration / time with all the reactions
        // net fraction: (time - time without reactions) / time with all the
        // reactions, i.e. the cost of the reaction alone
        const double t_all  = t_mean[n_configs - 1];
        const double t_none = t_mean[0];

        printf("%-10s %13s %10s %13s %10s %10s\n", "reactions", "mean [s]",
               "std [s]", "min [s]", "fraction", "net frac.");
        for (int c = 0; c < n_configs; c++)
        {
            const double net = (c == 0) ? t_none : t_mean[c] - t_none;
            printf("%-10s %13.6e %10.2e %13.6e %9.2f%% %9.2f%%\n", names[c],
                   t_mean[c], t_std[c], t_min[c], 100. * t_mean[c] / t_all,
                   100. * net / t_all);
        }

        double sum_net = t_none;
        for (int c = 1; c < n_configs - 1; c++)
        {
            sum_net += t_mean[c] - t_none;
        }
        printf("\n");
        printf("none + sum of the net single-reaction times = %.6e s "
               "(%.2f%% of 'all')\n",
               sum_net, 100. * sum_net / t_all);
        printf("\n");

        // Write the timings
        FILE* fout = fopen(output_file, "w");
        if (fout == NULL)
        {
            printf("Cannot open the output file %s\n", output_file);
        }
        else
        {
            fprintf(fout, "# execution space = %s, npts = %d, n_quad = %d, "
                          "n_rep = %d\n",
                    DevExeSpace::name(), npts, n_quad, n_rep);
            fprintf(fout, "# reactions, mean [s], std [s], min [s]\n");
            for (int c = 0; c < n_configs; c++)
            {
                fprintf(fout, "%s %.10e %.10e %.10e\n", names[c], t_mean[c],
                        t_std[c], t_min[c]);
            }
            fclose(fout);
            printf("Timings written to %s\n", output_file);
        }
    }
    Kokkos::finalize();

    return 0;
}
